#include <chrono>
#include <random>
#include <mpi.h>
#include <boost/program_options.hpp>
#include "./include/io.hpp"
#include "./include/myDist.hpp"
#include "./wrappers/blas_templates.hpp"
#include "./include/scalapack_utils.hpp"

const int i_zero = 0, i_one = 1;
const std::size_t sze_one = 1;

typedef std::size_t DESC[ 9 ];

template<typename T>
void showMatrix(T *mat, std::size_t m, std::size_t n, std::size_t lda){

        for(std::size_t i = 0; i < m; i ++){

            std::cout << "(" << i << ") " << " : ";

            for(std::size_t j = 0; j < n; j++){
                std::cout << std::scientific;
                std::cout << std::right << std::setw(16) << mat[i + j * lda] << "    ";

            }
            std::cout << std::endl;
        }
}

//gather block-cyclic matrix to rank 0
void GatherMatrix(int ctxt, int M, int N, int Mb, int Nb, int nrows, int ncols, double *A_loc, double *A_glob){

    int iam, nprocs;
    int procrows, proccols;
    int myrow, mycol;
    MPI_Datatype type1, type2;

    int TAG = 3;
    MPI_Status status;

    blacs_pinfo( &iam, &nprocs );

    int sendr = 0, sendc = 0, recvr = 0, recvc = 0;

    blacs_gridinfo( &ctxt, &procrows, &proccols, &myrow, &mycol);

    for (int r = 0; r < M; r += Mb, sendr = (sendr + 1) % procrows) {

        sendc = 0;

        int nr = Mb;

        if (M - r < Mb){
        
            nr = M - r;
	
	}
        
        for (int c = 0; c < N; c += Nb, sendc = (sendc + 1) % proccols){

	    int nc = Nb;

            if (N-c < Nb){

                nc = N - c;

	    }

            MPI_Type_vector(nc, nr, M, MPI_DOUBLE,&type1);
            MPI_Type_commit(&type1);
            MPI_Type_vector(nc, nr, nrows, MPI_DOUBLE,&type2);
            MPI_Type_commit(&type2);

            if (myrow == sendr && mycol == sendc) {

                MPI_Send(A_loc + nrows * recvc + recvr, 1, type2, 0, TAG, MPI_COMM_WORLD);

                recvc = (recvc+nc)%ncols;
            }

            if (myrow == 0 && mycol == 0) {

                MPI_Recv(A_glob + M * c + r, 1, type1, sendr + sendc * procrows, TAG, MPI_COMM_WORLD, &status);
            }

        }

	if (myrow == sendr){
         
	   recvr = (recvr+nr)%nrows;

	}
    }
}

namespace po = boost::program_options;

int main(int argc, char* argv[]){

  MPI_Init(&argc, &argv);

  po::options_description descp("Artificial Matrices with ScaLAPACK: Options");

  descp.add_options()
       ("help,h","show the help"
		 "Attention, for the current implementation of parallel IO, please make sure N/mbsize/dim0 == 0 and N/bbsize/dim1 == 0")
       ("N", po::value<std::size_t>()->default_value(10), "number of row and column of matrices to be generated.")
       ("dim0", po::value<int>()->default_value(1), "first dimension of 2D MPI cartesian grid.")    
       ("dim1", po::value<int>()->default_value(1), "second dimension of 2D MPI cartesian grid.")       
       ("mbsize", po::value<std::size_t>()->default_value(5), "ScaLAPACK block size in the first dimension of 2D MPI cartesian grid.")
       ("nbsize", po::value<std::size_t>()->default_value(5), "ScaLAPACK block size in the second dimension of 2D MPI cartesian grid.")
       ("dmax", po::value<double>()->default_value(1), "A scalar which scales the generated eigenvalues, this makes"
							 " the maximum absolute eigenvalue is abs(dmax)." )
       ("epsilon", po::value<double>()->default_value(0.1), "This value is epsilon." ) 
       ("myDist", po::value<std::size_t>()->default_value(0), "Specifies my externel setup distribution for generating eigenvalues:\n "
							      "0: Uniform eigenspectrum lambda_k = dmax * (epsilon + k * (1 - epsilon) / n for k = 0, ..., n-1)\n "
                                                              "1: Geometric eigenspectrum lambda_k =lambda_k = epsilon^[(n - k) / n] for k = 0, ..., n-1) \n")
       ("mean", po::value<double>()->default_value(0.5), "Mean value of Normal distribution for the randomness." )
       ("stddev", po::value<double>()->default_value(1.0), "Standard deviation value of Normal distribution for the randomness." )
  ;  

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, descp), vm);

  if (vm.count("help")) {
    std::cout << descp << std::endl;
    return 1;
  }
  
  //will be parsered by boost
  //2D dimmension of MPI grid
  int dim0 = vm["dim0"].as<int>();
  int dim1 = vm["dim1"].as<int>();  
  //size of matrix to be generated
  std::size_t N = vm["N"].as<std::size_t>();
  //block size of scalapack
  std::size_t mbsize = vm["mbsize"].as<std::size_t>();
  std::size_t nbsize = vm["nbsize"].as<std::size_t>();
  //for the randomness with normal distribution
  double mean = vm["mean"].as<double>();
  double stddev = vm["stddev"].as<double>();
  double dmax = vm["dmax"].as<double>();
  double eps = vm["epsilon"].as<double>();
  std::size_t dist = vm["myDist"].as<std::size_t>();

  int b0 = N / mbsize; 
  int b1 = N / nbsize;

  if(N % mbsize + N % nbsize + b0 % dim0 + b1 % dim1 != 0 ){
    std::cout << "Attention:\ndue to the simplication of implementation of parallel IO,\n"
		<< "N should be divisible by mbsize and nbsize\n"
		<<"N/mbsize and N / nbsize should be divisible by dim0 and dim1, respectively" << std::endl;
    return 1;
  }

  //BLACS part to setup the grid  
  int irsrc = i_zero;
  int icsrc = i_zero;
  int myproc, nprocs;
  blacs_pinfo( &myproc, &nprocs );
  //std::cout << "mprocs: " << myproc << ", nprocs: " << nprocs << std::endl; 
  int ictxt;
  int val;
  blacs_get( &ictxt, &i_zero, &val );
  blacs_gridinit( &ictxt, 'C', &dim0, &dim1 );
  int myrow, mycol;
  blacs_gridinfo( &ictxt, &dim0, &dim1, &myrow, &mycol);

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
/*
  for(std::size_t i = 0; i < blocknb_r; i++){
      for(std::size_t j = 0; j < blocknb_c; j++){
          std::cout << "[" << myrow << "," 
	  << mycol << "]: ("
	  << r_offs[i] << ":"
          << r_offs_l[i] << ":"		      
          << r_lens[i] << ","
          << c_offs[j] << ":"
          << c_offs_l[j] << ":"		     
	  << c_lens[j] << "),"
	  << std::endl;
    }
  }
*/
  //construct scalapack matrix descriptor 
  DESC   desc;
  int    info;

  descinit( desc, &N, &N, &mbsize, &nbsize, &i_zero, &i_zero, &ictxt, &lld_loc, &info );

  //generating ...
  if(myproc == 0) std::cout << "]> start generating ..." << std::endl;

  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::duration<double> elapsed;

  start = std::chrono::high_resolution_clock::now();

  //initialisation of a random matrix M
  auto M_loc_ptr = std::unique_ptr<double[]>(new double[N_loc_r * N_loc_c]);
  double *M_loc = M_loc_ptr.get();

  std::mt19937 generator(131421);
  std::normal_distribution<double> distribution(mean,stddev);
 
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
  std::vector<double> tau(std::min(N_loc_r, N_loc_c));

  t_pgeqrf<double>(N, N, M_loc, sze_one, sze_one, desc, tau.data());

  //initial matrix definition
  auto A_loc_ptr = std::unique_ptr<double[]>(new double[N_loc_r * N_loc_c]);  
  double *A_loc = A_loc_ptr.get();

  std::string mode;

  if(dist == 0){
    mode= "Uniform";
  }else if(dist == 1){
    mode = "Geometric";
  }
  /*else if(myDist == 2){
    mode = "1-2-1";
  }else if(myDist == 3){
    mode = "Wilkinson";
  }
  */

  //set the diagonal of A by given eigenvalues
  //Array of size m to store generated eigenvalues
  double *d = new double[N];

  if(dist == 0){
    myUniformDist<double>(d, N, eps, dmax);
  }else{
    myGeometricDist<double>(d, N, eps, dmax);
  }

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
  t_pmqr<double>('L','N', N, N, N, M_loc, sze_one, sze_one, desc, tau.data(), A_loc, sze_one, sze_one, desc);

  t_pmqr<double>('R','T', N, N, N, M_loc, sze_one, sze_one, desc, tau.data(), A_loc, sze_one, sze_one, desc);

  end = std::chrono::high_resolution_clock::now();

  elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
  
  if(myproc == 0) std::cout << "]> matrix generated in " << elapsed.count() << " seconds" << std::endl;
  
  //Write generated matrix A into binary file.
  MPI_File fh;
  MPI_Status status;
 
  std::ostringstream out_str;

  out_str << std::scientific << "matgen_m_" << N << "_" << mode << "_eps_"<< eps
            << ".bin";
 
  int handle;

  handle = MPI_File_open(MPI_COMM_WORLD, out_str.str().c_str(),
                  MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &fh);

  MPI_Datatype file_type;
  MPI_Datatype memory_type;

  start = std::chrono::high_resolution_clock::now();

  for(std::size_t j = 0; j < blocknb_c; j++){
    for(std::size_t i = 0; i < blocknb_r; i++){

      MPI_Type_vector(c_lens[j], r_lens[i], N, MPI_DOUBLE, &file_type);
      MPI_Type_commit(&file_type);

      MPI_Type_vector(c_lens[j], r_lens[i], N_loc_r, MPI_DOUBLE, &memory_type);
      MPI_Type_commit(&memory_type);

      MPI_Offset offset;
      offset = (r_offs[i] + c_offs[j] * N) * sizeof(double);

      MPI_File_set_view(fh, offset, MPI_DOUBLE, file_type, "native", MPI_INFO_NULL);
      MPI_File_write_all(fh, A_loc + r_offs_l[i] + c_offs_l[j] * N_loc_r, 1, memory_type, &status);

      MPI_Type_free(&memory_type);
      MPI_Type_free(&file_type);

    }
  }

  MPI_File_close( &fh );

  end = std::chrono::high_resolution_clock::now();

  elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);

  double gb = (double)N * (double)N * sizeof(double)/1024.0/1024.0/1024.0;
  if(myproc == 0) 
    std::cout << "]> matrix wrote in binary " << out_str.str() << " of " << gb << "GB in " 
	<< elapsed.count() << " seconds" << std::endl;

  //verification of writing
  if(N <= 10){
    if(myproc == 0){
      std::cout << "]> Verification of IO for block-cyclic matrix" << std::endl;
    }
    //1. gather matrix to rank 0
    std::vector<double> A_glob(N*N);
  
    GatherMatrix(ictxt, N, N, mbsize, nbsize, N_loc_r, N_loc_c, A_loc, A_glob.data());
 
    if(myproc == 0){
      showMatrix(A_glob.data(), N, N, N);
    }

    //2. read from the file
    std::vector<double> A_load(N*N);
    readMatFromBinary(A_load.data(), out_str.str(), N * N);
    if(myproc == 0){
      showMatrix(A_load.data(), N, N, N);
    }
    
    std::vector<double> diff;

    std::transform(A_glob.begin(), A_glob.end(), A_load.begin(), std::back_inserter(diff), [&](double l, double r)
    {
      return std::abs(l - r);
    });

    double max_abs = *std::max_element(diff.begin(), diff.end());

    if(myproc == 0){
      std::cout << "]> max difference between the gathered matrix and written matrix is: " << max_abs << std::endl;
    }

  }
  MPI_Finalize();	

}
