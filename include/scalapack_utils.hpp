#pragma once

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
                   std::size_t icsrc_, std::size_t* &r_offs, std::size_t* &r_lens,
                   std::size_t* &r_offs_l, std::size_t* &c_offs, std::size_t* &c_lens, std::size_t* &c_offs_l)
{

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

  r_offs = r_offs_.get();
  r_lens = r_lens_.get();
  r_offs_l = r_offs_l_.get();
  c_offs = c_offs_.get();
  c_lens = c_lens_.get();
  c_offs_l = c_offs_l_.get();

}


