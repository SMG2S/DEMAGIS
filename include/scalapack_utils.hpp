#pragma once

#include <mpi.h>

/*
 * Logic for an MPI parallel HEMM
 */

template <typename T>
MPI_Datatype getMPI_Type();

template <>
MPI_Datatype getMPI_Type<float>() {
  return MPI_FLOAT;
}

template <>
MPI_Datatype getMPI_Type<double>() {
  return MPI_DOUBLE;
}

template <>
MPI_Datatype getMPI_Type<std::complex<float> >() {
  return MPI_COMPLEX;
}

template <>
MPI_Datatype getMPI_Type<std::complex<double> >() {
  return MPI_DOUBLE_COMPLEX;
}


std::pair<std::size_t, std::size_t> numroc(std::size_t n, std::size_t nb, int iproc, int isrcproc, int nprocs){

    std::size_t numroc;
    std::size_t extrablks, mydist, nblocks;
    mydist = (nprocs + iproc - isrcproc) % nprocs;
    nblocks = n / nb;
    numroc = (nblocks / nprocs) * nb;
    extrablks = nblocks % nprocs;

    if(mydist < extrablks)
        numroc = numroc + nb;
    else if(mydist == extrablks)
        numroc = numroc + n % nb;

    std::size_t nb_loc = numroc / nb;

    if(numroc % nb != 0){
        nb_loc += 1;
    }

    return std::make_pair(numroc, nb_loc);
}


void get_offs_lens(std::size_t N_, std::size_t mb_, std::size_t nb_, int dim0, int dim1, int myrow, int mycol,
                   std::size_t mblocks_, std::size_t nblocks_, std::size_t irsrc_,
                   std::size_t icsrc_, std::size_t* &r_offs_, std::size_t* &r_lens_,
                   std::size_t* &r_offs_l_, std::size_t* &c_offs_, std::size_t* &c_lens_, std::size_t* &c_offs_l_)
{
/*
  std::unique_ptr<std::size_t[]> r_offs_;
  std::unique_ptr<std::size_t[]> r_lens_;
  std::unique_ptr<std::size_t[]> r_offs_l_;
  std::unique_ptr<std::size_t[]> c_offs_;
  std::unique_ptr<std::size_t[]> c_lens_;
  std::unique_ptr<std::size_t[]> c_offs_l_;

  r_offs_.reset(new std::size_t[mblocks_]());
  r_lens_.reset(new std::size_t[mblocks_]());
  r_offs_l_.reset(new std::size_t[mblocks_]());
  c_offs_.reset(new std::size_t[nblocks_]());
  c_lens_.reset(new std::size_t[nblocks_]());
  c_offs_l_.reset(new std::size_t[nblocks_]());
*/

  std::size_t sendr = irsrc_;
  std::size_t sendc = icsrc_;

  int cnt = 0;
  std::size_t nr, nc;

  for(std::size_t r = 0; r < N_; r += mb_, sendr = (sendr + 1) % dim0){
    nr = mb_;
    if (N_ - r < mb_){
      nr = N_ - r;
    }

    if(myrow == sendr){
      r_offs_[cnt] = r;
      r_lens_[cnt] = nr;
      cnt++;
    }
  }

  cnt = 0;

  for(std::size_t c = 0; c < N_; c += nb_, sendc = (sendc + 1) % dim1){
    nc = nb_;
    if(N_ - c < nb_){
      nc = N_ - c;
    }
    if(mycol == sendc){
      c_offs_[cnt] = c;
      c_lens_[cnt] = nc;
      cnt++;
    }
  }

  r_offs_l_[0] = 0;
  c_offs_l_[0] = 0;

  for(std::size_t i = 1; i < mblocks_; i++){
    r_offs_l_[i] = r_offs_l_[i - 1] + r_lens_[i - 1];
  }

  for(std::size_t j = 1; j < nblocks_; j++){
    c_offs_l_[j] = c_offs_l_[j - 1] + c_lens_[j - 1];
  }
/*
  r_offs = r_offs_.get();
  r_lens = r_lens_.get();
  r_offs_l = r_offs_l_.get();
  c_offs = c_offs_.get();
  c_lens = c_lens_.get();
  c_offs_l = c_offs_l_.get();
*/
}

//gather block-cyclic matrix to rank 0
template<typename T>
void GatherMatrix(int ctxt, std::size_t M, std::size_t N, std::size_t Mb, std::size_t Nb, 
			int nrows, int ncols, T *A_loc, T *A_glob){

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

            MPI_Type_vector(nc, nr, M, getMPI_Type<T>(),&type1);
            MPI_Type_commit(&type1);
            MPI_Type_vector(nc, nr, nrows, getMPI_Type<T>(),&type2);
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

//verification of writing
template<typename T>
Base<T> matDiff(int ctxt, T *A_loc, std::size_t N, std::size_t mbsize, std::size_t nbsize, std::string out_str){
    
    //BLACS part to setup the grid
    int myproc, nprocs;

    blacs_pinfo( &myproc, &nprocs );

    int dim0, dim1;
    int myrow, mycol;

    blacs_gridinfo( &ctxt, &dim0, &dim1, &myrow, &mycol);

    //get local size of matrix = N_loc_r x N_loc_c
    std::size_t N_loc_r, N_loc_c, blocknb_r, blocknb_c;
    std::tie(N_loc_r, blocknb_r) = numroc( N, mbsize, myrow, 0, dim0 );
    std::tie(N_loc_c, blocknb_c) = numroc( N, nbsize, mycol, 0, dim1 );

    std::vector<T> A_glob(N*N);

    GatherMatrix<T>(ctxt, N, N, mbsize, nbsize, N_loc_r, N_loc_c, A_loc, A_glob.data());

    if(myproc == 0){
      showMatrix<T>(A_glob.data(), N, N, N);
    }

    //2. read from the file
    std::vector<T> A_load(N*N);
    readMatFromBinary<T>(A_load.data(), out_str, N * N);
    if(myproc == 0){
      showMatrix(A_load.data(), N, N, N);
    }
    
    std::vector<Base<T>> diff;

    std::transform(A_glob.begin(), A_glob.end(), A_load.begin(), std::back_inserter(diff), [&](T l, T r)
    {
      return std::abs(l - r);
    });

    Base<T> max_abs = *std::max_element(diff.begin(), diff.end());

    if(myproc == 0){
      std::cout << "]> max difference between the gathered matrix and written matrix is: " << max_abs << std::endl;
    }
    
    return max_abs;
    
}

