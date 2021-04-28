#include <chrono>
#include <mpi.h>
#include <boost/program_options.hpp>
#include "./include/io.hpp"
#include "./wrappers/blas_templates.hpp"

const int i_zero = 0, i_one = 1;
const std::size_t sze_one = 1;

typedef std::size_t DESC[ 9 ];

int main(int argc, char* argv[]){

  MPI_Init(&argc, &argv);

  //will be parsered by boost
  //2D dimmension of MPI grid
  int dim0=1, dim1=1;
  //size of matrix to be generated
  std::size_t N = 100;
  //block size of scalapack
  std::size_t mbsize = 10;
  std::size_t nbsize = 10;

  int irsrc = i_zero;
  int icsrc = i_zero;

  //BLACS part to setup the grid
  int myproc, nprocs;
  blacs_pinfo( &myproc, &nprocs );
  std::cout << "mprocs: " << myproc << ", nprocs: " << nprocs << std::endl; 
  int ictxt;
  int val;
  blacs_get( &ictxt, &i_zero, &val );
  blacs_gridinit( &ictxt, 'C', &dim0, &dim1 );
  int myrow, mycol;
  blacs_gridinfo( &ictxt, &dim0, &dim1, &myrow, &mycol);

  //get local size of matrix = N_loc_r x N_loc_c
  std::size_t N_loc_r, N_loc_c;
  N_loc_r = numroc( &N, &mbsize, &myrow, &irsrc, &dim0 );
  N_loc_c = numroc( &N, &nbsize, &mycol, &icsrc, &dim1 );

  //for column major matrix, the leading dimension
  std::size_t lld_loc = std::max(N_loc_r, (std::size_t)1);

  //construct scalapack matrix descriptor 
  DESC   desc;
  int    info;

  descinit( desc, &N, &N, &mbsize, &nbsize, &irsrc, &irsrc, &ictxt, &lld_loc, &info );

  //initialisation of a random matrix M
  auto M_loc_ptr = std::unique_ptr<double[]>(new double[N_loc_r * N_loc_c]);
  double *M_loc = M_loc_ptr.get();

  //this part will be replaced by the random generation later
  for(auto i = 0; i < N_loc_r; i++){
    for(auto j = 0; j < N_loc_c; j++){
      M_loc[i + j * lld_loc] = 5;
    }
  }

  //QR factorization of M
  std::vector<double> tau(std::min(N_loc_r, N_loc_c));

  t_pgeqrf<double>(N, N, M_loc, sze_one, sze_one, desc, tau.data());

  //initial matrix definition
  auto A_loc_ptr = std::unique_ptr<double[]>(new double[N_loc_r * N_loc_c]);  
  double *A_loc = A_loc_ptr.get();

  //set the diagonal of A by given eigenvalues


  //A = Q*A*Q^T
  t_pmqr<double>('L','N', N, N, N, M_loc, sze_one, sze_one, desc, tau.data(), A_loc, sze_one, sze_one, desc);

  t_pmqr<double>('R','T', N, N, N, M_loc, sze_one, sze_one, desc, tau.data(), A_loc, sze_one, sze_one, desc);
  

  //Write generated matrix A into binary file.
  
  MPI_Finalize();	

}
