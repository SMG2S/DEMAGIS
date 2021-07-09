#pragma once

#include "scalapack_utils.hpp"

template <typename T>
void wrtMatIntoBinaryMPI(int ctxt, T *H, std::string path_out, std::size_t N, 
		std::size_t mbsize, std::size_t nbsize){
  
    std::ostringstream problem(std::ostringstream::ate);
    problem << path_out;

    std::size_t b0 = N / mbsize;
    std::size_t b1 = N / nbsize;

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

    //for column major matrix, the leading dimension
    std::size_t lld_loc = std::max(N_loc_r, (std::size_t)1);

    std::size_t *r_offs = new std::size_t[blocknb_r];
    std::size_t *r_lens = new std::size_t[blocknb_r];
    std::size_t *r_offs_l = new std::size_t[blocknb_r];

    std::size_t *c_offs = new std::size_t[blocknb_c];
    std::size_t *c_lens = new std::size_t[blocknb_c];
    std::size_t *c_offs_l = new std::size_t[blocknb_c];

    get_offs_lens(N, mbsize, nbsize, dim0, dim1, myrow, mycol, blocknb_r,
                  blocknb_c, 0, 0, r_offs, r_lens, r_offs_l, c_offs,
                  c_lens, c_offs_l);

    //Write generated matrix A into binary file.
    MPI_File fh;
    MPI_Status status;

    int handle;

    handle = MPI_File_open(MPI_COMM_WORLD, path_out.c_str(),
                  MPI_MODE_CREATE | MPI_MODE_WRONLY,
                  MPI_INFO_NULL, &fh);

    MPI_Datatype file_type;
    MPI_Datatype memory_type;

    for(std::size_t j = 0; j < blocknb_c; j++){

    	for(std::size_t i = 0; i < blocknb_r; i++){

            MPI_Type_vector(c_lens[j], r_lens[i], N, getMPI_Type<T>(), &file_type);
            MPI_Type_commit(&file_type);

            MPI_Type_vector(c_lens[j], r_lens[i], N_loc_r, getMPI_Type<T>(), &memory_type);
            MPI_Type_commit(&memory_type);

	    MPI_Offset offset;
            offset = (r_offs[i] + c_offs[j] * N) * sizeof(T);

            MPI_File_set_view(fh, offset, getMPI_Type<T>(), file_type, "native", MPI_INFO_NULL);
            MPI_File_write_all(fh, H + r_offs_l[i] + c_offs_l[j] * N_loc_r, 1, memory_type, &status);

            MPI_Type_free(&memory_type);
            MPI_Type_free(&file_type);

        }
    }  
    
    MPI_File_close( &fh );


}

