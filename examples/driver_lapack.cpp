#include <chrono>
#include <random>
#include <boost/program_options.hpp>
#include "myDist.hpp"
#include "matGen_lapack.hpp"

namespace po = boost::program_options;

int main(int argc, char* argv[]){

  po::options_description desc("Artificial Matrices MT: Options");
  desc.add_options()
       ("help,h","show the help")
       ("N", po::value<std::size_t>()->default_value(10), "number of row and column of matrices to be generated.")
       ("dmax", po::value<double>()->default_value(21), "A scalar which scales the generated eigenvalues, this makes"
							 " the maximum absolute eigenvalue is abs(dmax)." )
       ("epsilon", po::value<double>()->default_value(0.1), "This value is epsilon." ) 

       ("myDist", po::value<std::size_t>()->default_value(0), "Specifies my externel setup distribution for generating eigenvalues:\n "
							      "0: Uniform eigenspectrum lambda_k = dmax * (epsilon + k * (1 - epsilon) / n for k = 0, ..., n-1)\n "
                                                              "1: Geometric eigenspectrum lambda_k =lambda_k = epsilon^[(n - k) / n] for k = 0, ..., n-1) \n"
							      "2: 1-2-1 matrix\n"
							      "3: Wilkinson matrix\n")
       ("mean", po::value<double>()->default_value(0.5), "Mean value of Normal distribution for the randomness." )
       ("stddev", po::value<double>()->default_value(1.0), "Standard deviation value of Normal distribution for the randomness." );

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }
  
  //number of row and column of matrix to be generated
  std::size_t n = vm["N"].as<std::size_t>();
  double dmax = vm["dmax"].as<double>();
  double eps = vm["epsilon"].as<double>();
  std::size_t myDist = vm["myDist"].as<std::size_t>();
  //for the randomness with normal distribution
  double mean = vm["mean"].as<double>();
  double stddev = vm["stddev"].as<double>();

  //generating ...
  std::cout << "]> start generating ..." << std::endl;

  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::duration<double> elapsed;
  
  start = std::chrono::high_resolution_clock::now();

  double *A;
  
  //myDist = 0: Uniform eigenspectrum
  //myDist = 1: Geometric eigenspectrum
  //myDist = 2: 1-2-1 tridiagonal matrix
  //myDist = 3: Wilkinson tridiagonal matrix
  std::string mode;

  if(myDist == 0){
    mode= "Uniform";
  }else if(myDist == 1){
    mode = "Geometric";
  }else if(myDist == 2){
    mode = "1-2-1";
  }else if(myDist == 3){
    mode = "Wilkinson";
  }else{
    mode = "myLambda";
  }

  if(myDist == 0){
      A = matGen_lapack<double>(n, mean, stddev, myUniformDist<double>, n, eps, dmax);
  }else if(myDist == 1){
      A = matGen_lapack<double>(n, mean, stddev, myGeometricDist<double>, n, eps, dmax);
  }else if(myDist == 2){
      A = matGen_121<double>(n);
  }else if(myDist == 3){
      A = matGen_WilkinsonPlus<double>(n); 
  }else{
    A = matGen_lapack<double>(n, mean, stddev, [](std::size_t n, int x){
                double *eigenv = new double[n];
                for(auto k = 0; k < n; k++){
                    eigenv[k] = k * (x + 1);
                }
                return eigenv;
            }, n, 1);
  }

  end = std::chrono::high_resolution_clock::now();

  elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
  
  std::cout << "]> matrix generated in " << elapsed.count() << " seconds" << std::endl;

  std::ostringstream out_str;

  if(myDist == 0 || myDist == 1){
    out_str << std::scientific << "matgen_m_" << n << "_" << mode << "_eps_"<< eps
            << "_dmax_" << dmax << ".bin";
  }else{
    out_str << std::scientific << "matgen_m_" << n << "_"<< mode << ".bin";
  }

  //save matrix into binary file
  wrtMatIntoBinary<double>(A, out_str.str(), n * n);
 
  delete [] A; 
  return 0;
}

