#include <chrono>
#include <random>
#include <boost/program_options.hpp>
#include "./include/io.hpp"
#include "./include/myDist.hpp"
#include "./wrappers/blas_templates.hpp"

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


namespace po = boost::program_options;

int main(int argc, char* argv[]){

  po::options_description desc("Artificial Matrices MT: Options");
  desc.add_options()
       ("help,h","show the help")
       ("size", po::value<std::size_t>()->default_value(10), "number of row and column of matrices to be generated.")
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
  std::size_t n = vm["size"].as<std::size_t>();
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

  //initialisation of a random matrix M
  auto M_ptr = std::unique_ptr<double[]>(new double[n * n]);
  auto *M = M_ptr.get();

  std::mt19937 generator(131421);
  std::normal_distribution<double> distribution(mean,stddev);

  for(auto i = 0; i < n * n; i++){
    auto rnd = distribution(generator);
    M[i] = rnd;    
  }
 
  if(n <= 10){
    std::cout << "Show Matrix M: "<< std::endl;
    showMatrix<double>(M, n, n, n);
  }

  //QR factorization of M
  std::vector<double> tau(n);
  t_geqrf<double>(n, n, M, n, tau.data());

  //initial matrix definition
  auto A_ptr = std::unique_ptr<double[]>(new double[n * n]);  
  auto *A = A_ptr.get();
  for(auto i = 0; i < n * n; i++){
    A[i] = 0.0;
  }

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
  }

  double *d;

  if(myDist == 0 || myDist == 1){
    d = (double *)malloc(n * sizeof(double));
    if(myDist == 0){
      myUniformDist<double>(d, n, eps, dmax);
    }else{
      myGeometricDist<double>(d, n, eps, dmax);
    }
    for (auto i = 0; i < n; ++i) {
      A[i + n * i] = d[i];
    }
  }else if(myDist == 2){
    for (auto i = 0; i < n; ++i) {
      A[i + n * i] = 2.0;
      if (i != n - 1) {
 	A[i + 1 + n * i] = 1.0;
        A[i + n * (1 + i)] = 1.0;	
      }

    } 
  }else if(myDist == 3){
    double m = ((double)n - 1.0) * 0.5;
    for (auto i = 0; i < n; ++i) {
      double a1 = double(n - 1 - i); double a2 = 0. - i;
      int v = std::min(abs(a1), abs(a2));
      A[i + n * i] = m - v;
      if (i != n - 1) {
        A[i + 1 + n * i] = 1.0;
        A[i + n * (1 + i)] = 1.0;
      }
    }
  }else{
    std::cout << "myDist's value should be one of 0,1,2,3\n";
    return 1;
  }
 
  if(n <= 10){
    std::cout << "Show Initial Matrix A: "<< std::endl;
    showMatrix<double>(A, n, n, n);
  }

  //  A := Q*A*Q^T
  t_mqr<double>('L','N', n, n, n, M, n, tau.data(), A, n);
  t_mqr<double>('R','T', n, n, n, M, n, tau.data(), A, n);

  if(n <= 10){
    std::cout << "Show Final Matrix A: "<< std::endl;
    showMatrix<double>(A, n, n, n);
  }

  end = std::chrono::high_resolution_clock::now();

  elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
  
  std::cout << "]> matrix generated in " << elapsed.count() << " seconds" << std::endl;

  std::ostringstream out_str;

  if(myDist == 0 || myDist == 1){
    out_str << std::scientific << "matgen_m_" << n << "_" << mode << "_eps_"<< eps
            << ".bin";
  }else{
    out_str << std::scientific << "matgen_m_" << n << "_"<< mode << ".bin";
  }

  //save matrix into binary file
  wrtMatIntoBinary<double>(A, out_str.str(), n * n);
  
  if(myDist == 0 || myDist == 1){
    delete[] d;
  }
  
  return 0;
}

