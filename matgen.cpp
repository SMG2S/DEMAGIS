#include <chrono>
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
/*
template<typename T>
void myUniformDist(Base<T>* &eigenv, std::size_t n, Base<T> epsilon, Base<T> dmax){
  for(auto k = 0; k < n; k++){
    eigenv[k] = dmax * (epsilon + (Base<T>)k * (1.0 - epsilon) / (Base<T>)n);
  }
}

template<typename T>
void myGeometricDist(Base<T>* &eigenv, std::size_t n, Base<T> epsilon){
  for(auto k = 0; k < n; k++){
    Base<T> exp = (Base<T>)(n - k) / (Base<T>)n; 
    eigenv[k] = std::pow(epsilon, exp);
  }
}
*/

namespace po = boost::program_options;

int main(int argc, char* argv[]){

  po::options_description desc("Artificial Matrices: Options");

  desc.add_options()
       ("help,h","show the help")
       ("size", po::value<std::size_t>()->default_value(10), "number of row and column of matrices to be generated.")
       ("mode", po::value<std::size_t>()->default_value(1), "Describes how the singular/eigenvalues are specified.")
       ("cond", po::value<double>()->default_value(100), "Condition number of generated matrices, this parameter affect only" 
							  " when using internally provided distribution to generate eigenvalues.")
       ("dmax", po::value<double>()->default_value(21), "A scalar which scales the generated eigenvalues, this makes"
							 " the maximum absolute eigenvalue is abs(dmax)" )
       ("klu", po::value<std::size_t>()->default_value(9), "Specifies the upper/lower bandwidth of the matrix.")
       ("epsilon", po::value<double>()->default_value(0.1), "This value affects only if mode == 0, in which eigenvalues are " 
							    "generated explicitly by users" )       
       ("dist", po::value<std::string>()->default_value("U"), "Specifies the internal distribution for generating eigenvalues: "
							      "'U': uniform distribution (0, 1); 'S': symmetric uniform distribution (-1, 1)"
							      "'N': normal distribution (0, 1)")
       ("myDist", po::value<std::size_t>()->default_value(0), "Specifies my externel setup distribution for generating eigenvalues:\n "
							      "0: Uniform eigenspectrum lambda_k = dmax * (epsilon + k * (1 - epsilon) / n for k = 0, ..., n-1)\n "
                                                              "1: Geometric eigenspectrum lambda_k =lambda_k = epsilon^[(n - k) / n] for k = 0, ..., n-1) \n");
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }

  //number of row and column of matrix to be generated
  std::size_t m = vm["size"].as<std::size_t>();
  //Describes how the singular/eigenvalues are specified.
  std::size_t mode = vm["mode"].as<std::size_t>();
  //condition number for the generated matrices;
  double cond = vm["cond"].as<double>();
  double dmax = vm["dmax"].as<double>();
  std::size_t klu = vm["klu"].as<std::size_t>();
  double eps = vm["epsilon"].as<double>();
  //Internal distribution for generating eigenvalues
  char dist = vm["dist"].as<std::string>().c_str()[0];
  std::size_t myDist = vm["myDist"].as<std::size_t>();

  if(klu > m - 1 || klu < 0){
    std::cout << "klu value should in [0,m-1]." << std::endl;
    return 1;
  }

  if(mode < -6 && mode > 6){
    std::cout << "mode value must be one of (-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6)." << std::endl;
    return 1;
  }

  //Specifies packing of matrix
  char pack = 'N';
  //Specifies the seed of the random number generator
  std::size_t iseed[4];
  //iseed needs to be odd
  iseed[0] = 92; iseed[1] = 203; iseed[2] = 3921; iseed[3] = 1393;
  //Generate symmetric matrices, whose eigenvalues can be postive/negative.  
  char sym = 'S';

  //Array of size m to store generated eigenvalues
  double *d = new double[m];

  //if mode = 0, using myDist
  if(mode == 0){
    if(myDist == 0){
      myUniformDist<double>(d, m, eps, dmax);
    }else{
      myGeometricDist<double>(d, m, eps);
    }
  }

  //Matrix generated
  double *a = new double[m * m];
  std::size_t lda = m;

  //generating ...
  std::cout << "]> start generating ..." << std::endl;

  std::chrono::high_resolution_clock::time_point start, end;
  std::chrono::duration<double> elapsed;
  
  start = std::chrono::high_resolution_clock::now();

  t_latms<double>(m, m, dist, iseed, sym, d, mode, cond, dmax, klu, klu, pack, a, lda);

  end = std::chrono::high_resolution_clock::now();

  elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
  
  std::cout << "]> matrix generated in " << elapsed.count() << " seconds" << std::endl;

  //if matrix is small, show it 
  if(m <= 10)
    showMatrix<double>(a, m, m, lda);

  //path out to save the matrices into binary
  std::ostringstream out_str;
  if(mode == 0){
    out_str << std::scientific << "matgen_m_" << m << "_mydist_" << myDist << "_eps_"<< eps << "_mode_"
            << mode << "_dmax_" << dmax << "_klu_" << klu
            << ".bin";
  }else{
    out_str << std::scientific << "matgen_m_" << m << "_dist_" << dist << "_mode_"
            << mode << "_cond_"<< cond << "_dmax_" << dmax << "_klu_" << klu
            << ".bin";
  }  
  //save matrix into binary file
  wrtMatIntoBinary<double>(a, out_str.str(), m * m);

  
  delete[] d;
  delete[] a;

  return 0;

}
