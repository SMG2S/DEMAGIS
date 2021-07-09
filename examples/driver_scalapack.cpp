#include <chrono>
#include <random>
#include <mpi.h>
#include <boost/program_options.hpp>
#include "io.hpp"
#include "myDist.hpp"
#include "matGen_scalapack.hpp"
#include "io_mpi.hpp"
#include "blas_templates.hpp"
#include "scalapack_utils.hpp"


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

  std::string mode;

  if(dist == 0){
    mode= "Uniform";
  }else if(dist == 1){
    mode = "Geometric";
  }
  
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int myproc, nprocs;
  blacs_pinfo( &myproc, &nprocs );

  int ictxt;
  int val;
  blacs_get( &ictxt, &i_zero, &val );
  blacs_gridinit( &ictxt, 'C', &dim0, &dim1 );


  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::duration<double> elapsed;

  start = std::chrono::high_resolution_clock::now();
 
  double *A;
  
  if(dist == 0){
      A = matGen_scalapack<double>(ictxt, N, mbsize, nbsize,
                        	       mean, stddev, myUniformDist<double>, N, eps, dmax);
  }else if(dist == 1){
      A = matGen_scalapack<double>(ictxt, N, mbsize, nbsize,
                                       mean, stddev, myGeometricDist<double>, N, eps, dmax);  
  }else{
      A = matGen_scalapack<double>(ictxt, N, mbsize, nbsize, mean, stddev, 
		      		[](std::size_t n, int x){
                		    double *eigenv = new double[n];
                		    for(auto k = 0; k < n; k++){
                    			eigenv[k] = k * (x + 1);
                		    }
                		    return eigenv;
            		        }, N, 1);
  
  }

  end = std::chrono::high_resolution_clock::now();

  elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
  
  if(rank == 0) std::cout << "]> matrix generated in " << elapsed.count() << " seconds" << std::endl;
 
  std::ostringstream out_str;

  out_str << std::scientific << "matgen_m_" << N << "_" << mode << "_eps_"<< eps << "_dmax_" << dmax
            << ".bin";

  start = std::chrono::high_resolution_clock::now();

  wrtMatIntoBinaryMPI<double>(ictxt, A, out_str.str(), N, mbsize, nbsize);

  end = std::chrono::high_resolution_clock::now();

  elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);

  double gb = (double)N * (double)N * sizeof(double)/1024.0/1024.0/1024.0;

  if(rank == 0) 
    std::cout << "]> matrix wrote in binary " << out_str.str() << " of " << gb << "GB in " 
	<< elapsed.count() << " seconds" << std::endl;

  if(N <= 10){
      double matdiff = matDiff(ictxt, A, N, mbsize, nbsize, out_str.str());
  } 

  delete[] A;

  MPI_Finalize();	

}
