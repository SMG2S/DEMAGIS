#include "./include/io.hpp"
#include "./wrappers/blas_templates.hpp"

template<typename T>
void showMatrix(T *mat, std::size_t m, std::size_t n, std::size_t lda){

	for(std::size_t i = 0; i < m; i ++){

	    std::cout << "(" << i << ") " << " : ";

    	    for(std::size_t j = 0; j < n; j++){
	        std::cout << std::setprecision(8);
                std::cout << std::scientific;
                std::cout << std::right << std::setw(12) << mat[i + j * lda] << "    ";

	    }
	    std::cout << std::endl;
	}

	std::cout << "..." << std::endl;
}



int main(int argc, char* argv[]){

  //number of row and column of matrix to be generated
  std::size_t m = 10;
  //Internal distribution for generating eigenvalues
  char dist = 'U';
  //Specifies the seed of the random number generator
  std::size_t iseed[4];
  iseed[0] = 92; iseed[1] = 203; iseed[2] = 3921; iseed[3] = 1393;
  
  char sym = 'S';
  //Array of size m to store generated eigenvalues
  double *d = new double[m];

  //Describes how the singular/eigenvalues are specified.
  std::size_t mode = 1;

  //cond number for the generated matrices;
  double cond = 100;

  double dmax = 21;

  //Specifies the lower bandwidth of the matrix
  std::size_t kl = m - 1;
  //Specifies the upper bandwidth of the matrix
  std::size_t ku = m - 1;
  //Specifies packing of matrix
  char pack = 'N';

  std::size_t lda = m;

  //Matrix generated
  double *a = new double[m * m];
   
  t_latms<double>(m, m, dist, iseed, sym, d, mode, cond, dmax, kl, ku, pack, a, lda);
 
  showMatrix<double>(a, m, m, lda);

  delete[] d;
  delete[] a;

  return 0;

}
