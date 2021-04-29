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

namespace po = boost::program_options;

int main(int argc, char* argv[]){

  MPI_Init(&argc, &argv);

  po::options_description descp("Artificial Matrices with ScaLAPACK: Options");

  descp.add_options()
       ("help,h","show the help")
       ("N", po::value<std::size_t>()->default_value(10), "number of row and column of matrices to be generated.")
       ("dim0", po::value<int>()->default_value(1), "first dimension of 2D MPI cartesian grid.")    
       ("dim1", po::value<int>()->default_value(1), "second dimension of 2D MPI cartesian grid.")       
       ("mbsize", po::value<std::size_t>()->default_value(5), "ScaLAPACK block size in the first dimension of 2D MPI cartesian grid.")
       ("nbsize", po::value<std::size_t>()->default_value(5), "ScaLAPACK block size in the second dimension of 2D MPI cartesian grid.")
       ("dmax", po::value<double>()->default_value(21), "A scalar which scales the generated eigenvalues, this makes"
							 " the maximum absolute eigenvalue is abs(dmax)." )
       ("epsilon", po::value<double>()->default_value(0.1), "This value is epsilon." ) 
       ("dist", po::value<std::size_t>()->default_value(0), "Specifies my externel setup distribution for generating eigenvalues:\n "
							      "0: Uniform eigenspectrum lambda_k = dmax * (epsilon + k * (1 - epsilon) / n for k = 0, ..., n-1)\n "
                                                              "1: Geometric eigenspectrum lambda_k =lambda_k = epsilon^[(n - k) / n] for k = 0, ..., n-1) \n")
       ("mean", po::value<double>()->default_value(5.0), "Mean value of Normal distribution for the randomness." )
       ("stddev", po::value<double>()->default_value(2.0), "Standard deviation value of Normal distribution for the randomness." )
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
  std::size_t dist = vm["dist"].as<std::size_t>();

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

  //set the diagonal of A by given eigenvalues
  //Array of size m to store generated eigenvalues
  double *d = new double[N];

  if(dist == 0){
    myUniformDist<double>(d, N, eps, dmax);
  }else{
    myGeometricDist<double>(d, N, eps);
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
  
  MPI_Finalize();	

}
