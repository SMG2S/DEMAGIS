#pragma once

#include <complex>
#include "./types.hpp"

template <typename T>
void t_latms(const std::size_t m, const std::size_t n,
	     const char dist, std::size_t* iseed, char sym, Base<T> *d, const std::size_t mode,
	     const Base<T> cond, const Base<T> dmax, const std::size_t kl, const std::size_t ku,
	     char pack, T* a, const std::size_t lda);

#include "blas_templates.inc"


