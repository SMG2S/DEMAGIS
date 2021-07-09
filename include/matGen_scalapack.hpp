#include <random>
#include <mpi.h>

#include "io.hpp"
#include "myDist.hpp"
#include "blas_templates.hpp"
#include "scalapack_utils.hpp"

const int i_zero = 0, i_one = 1;
const std::size_t sze_one = 1;

typedef std::size_t DESC[ 9 ];

template<typename T, class Fn, typename... Ts>
T* matGen_scalapack(int ctxt, std::size_t N, std::size_t mbsize, std::size_t nbsize,
			Base<T> mean, Base<T> stddev, Fn fn, Ts... args) {

    std::size_t b0 = N / mbsize;
    std::size_t b1 = N / nbsize;

    //BLACS part to setup the grid  
    int irsrc = i_zero;
    int icsrc = i_zero;
    int myproc, nprocs;

    blacs_pinfo( &myproc, &nprocs );

    int dim0, dim1;
    int myrow, mycol;

    blacs_gridinfo( &ctxt, &dim0, &dim1, &myrow, &mycol);

    try{
        if(N % mbsize + N % nbsize + b0 % dim0 + b1 % dim1 != 0 ){
            throw std::logic_error(std::string("Attention:\ndue to the simplication of implementation of parallel IO,\n")
                          + std::string("N should be divisible by mbsize and nbsize")
                          + std::string("N/mbsize and N / nbsize should be divisible by dim0 and dim1, respectively\n"));
        }
    }
    catch(std::exception &e)
    {
        std::cerr << "Caught " << typeid( e ).name( ) << " : "<< e.what( ) << std::endl;

        throw;
    }

    //get local size of matrix = N_loc_r x N_loc_c
    std::size_t N_loc_r, N_loc_c, blocknb_r, blocknb_c;
    std::tie(N_loc_r, blocknb_r) = numroc( N, mbsize, myrow, irsrc, dim0 );
    std::tie(N_loc_c, blocknb_c) = numroc( N, nbsize, mycol, icsrc, dim1 );

    //for column major matrix, the leading dimension
    std::size_t lld_loc = std::max(N_loc_r, (std::size_t)1);

    std::size_t *r_offs = new std::size_t[blocknb_r];
    std::size_t *r_lens = new std::size_t[blocknb_r];
    std::size_t *r_offs_l = new std::size_t[blocknb_r];

    std::size_t *c_offs = new std::size_t[blocknb_c];
    std::size_t *c_lens = new std::size_t[blocknb_c];
    std::size_t *c_offs_l = new std::size_t[blocknb_c];

    get_offs_lens(N, mbsize, nbsize, dim0, dim1, myrow, mycol, blocknb_r, 
  		  blocknb_c, irsrc, icsrc, r_offs, r_lens, r_offs_l, c_offs, 
		  c_lens, c_offs_l);

    //construct scalapack matrix descriptor 
    DESC   desc;
    int    info;

    descinit( desc, &N, &N, &mbsize, &nbsize, &i_zero, &i_zero, &ctxt, &lld_loc, &info );

    //generating ...
    if(myproc == 0) std::cout << "]> start generating ..." << std::endl;

    //initialisation of a random matrix M
    auto M_loc_ptr = std::unique_ptr<T[]>(new T[N_loc_r * N_loc_c]);
    T *M_loc = M_loc_ptr.get();

    std::mt19937 generator(131421);
    std::normal_distribution<T> distribution(mean,stddev);
    
    std::vector<std::size_t> g_offs_r;
    std::vector<std::size_t> g_offs_c;  

    g_offs_r.push_back(0);
    g_offs_c.push_back(0);

    std::size_t r = 0, c = 0;
    for(auto i = 0; i < dim0 - 1; i++){
        r += numroc( &N, &mbsize, &i, &irsrc, &dim0 );
        g_offs_r.push_back(r);
    }

    for(auto j = 0; j < dim1 - 1; j++){
        c += numroc( &N, &nbsize, &j, &icsrc, &dim1 );
        g_offs_c.push_back(c);
    }

    int cnt = 0;
    for(auto j = 0; j < N; j++){
        for(auto i = 0; i < N; i++){
            auto rnd = distribution(generator);
            if((i >= g_offs_r[myrow]) && (i < (g_offs_r[myrow] + N_loc_r)) && (j >= g_offs_c[mycol]) && (j < (g_offs_c[mycol] + N_loc_c))){
	        M_loc[cnt] = rnd;
                cnt++;	
            }
        }
    }  


    //QR factorization of M
    std::vector<T> tau(std::min(N_loc_r, N_loc_c));
    t_pgeqrf<T>(N, N, M_loc, sze_one, sze_one, desc, tau.data());

    //initial matrix definition
    T *A_loc = new T[N_loc_r * N_loc_c];	

    T *d = fn(args...);

    for(std::size_t j = 0; j < blocknb_c; j++){
        for(std::size_t i = 0; i < blocknb_r; i++){
            for(std::size_t q = 0; q < c_lens[j]; q++){
                for(std::size_t p = 0; p < r_lens[i]; p++){
	            if(q + c_offs[j] == p + r_offs[i]){
	    	        A_loc[(q + c_offs_l[j]) * N_loc_r + p + r_offs_l[i]] = d[q + c_offs[j]];
	  	    }else{
	   		A_loc[(q + c_offs_l[j]) * N_loc_r + p + r_offs_l[i]] = 0.0;
	  	    }
		}
      	    }
    	}
    }

    //A = Q*A*Q^T
    t_pmqr<T>('L','N', N, N, N, M_loc, sze_one, sze_one, desc, tau.data(), A_loc, sze_one, sze_one, desc);
    t_pmqr<T>('R','T', N, N, N, M_loc, sze_one, sze_one, desc, tau.data(), A_loc, sze_one, sze_one, desc);
    
    delete[] d;

    return A_loc;

}
