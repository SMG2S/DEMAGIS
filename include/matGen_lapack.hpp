#pragma once

#include <chrono>
#include <random>
#include "./io.hpp"
#include "../wrappers/blas_templates.hpp"

template<typename T, class Fn, typename... Ts>
T* matGen_lapack(std::size_t n, Base<T> mean, Base<T> stddev, Fn fn, Ts... args) {	

    auto M_ptr = std::unique_ptr<T[]>(new T[n * n]);

    auto *M = M_ptr.get();

    std::mt19937 generator(131421);
    std::normal_distribution<double> distribution(mean,stddev);

    for(auto i = 0; i < n * n; i++){
        auto rnd = distribution(generator);
        M[i] = rnd;    
    }

    if(n <= 10){
        std::cout << "Show Matrix M: "<< std::endl;
        showMatrix<T>(M, n, n, n);
    }

    //QR factorization of M
    std::vector<T> tau(n);
    t_geqrf<T>(n, n, M, n, tau.data());

    //initial matrix definition
    auto A_ptr = std::unique_ptr<T[]>(new T[n * n]);  
    auto *A = A_ptr.get();
    for(auto i = 0; i < n * n; i++){
      A[i] = 0.0;
    }

    //generating spectrum
    T *d;
    d = fn( args...);

    std::cout << "n = " << n << std::endl;

    for(auto i = 0; i < n; i++){
            std::cout << d[i] << ",";
    }
    std::cout << std::endl;

    for (auto i = 0; i < n; ++i) {
      A[i + n * i] = d[i];
    }
    
    if(n <= 10){
        std::cout << "Show Initial Matrix A: "<< std::endl;
        showMatrix<T>(A, n, n, n);
    }

    // A := Q*A*Q^T
    t_mqr<T>('L','N', n, n, n, M, n, tau.data(), A, n);
    t_mqr<T>('R','T', n, n, n, M, n, tau.data(), A, n);

    if(n <= 10){
        std::cout << "Show Final Matrix A: "<< std::endl;
        showMatrix<T>(A, n, n, n);
    }

    delete[] d;

    return A;
}


template<typename T>
T *matGen_121(std::size_t n){

    auto A_ptr = std::unique_ptr<T[]>(new T[n * n]);
    auto *A = A_ptr.get();
    for(auto i = 0; i < n * n; i++){
      A[i] = 0.0;
    }

    for (auto i = 0; i < n; ++i) {
        A[i + n * i] = 2.0;
        if (i != n - 1) {
 	    A[i + 1 + n * i] = 1.0;
            A[i + n * (1 + i)] = 1.0;	
        }
    }

    if(n <= 10){
        std::cout << "Show Final Matrix A: "<< std::endl;
        showMatrix<T>(A, n, n, n);
    }

    return A;
}

template<typename T>
T *matGen_WilkinsonPlus(std::size_t n){

    auto A_ptr = std::unique_ptr<T[]>(new T[n * n]);
    auto *A = A_ptr.get();
    for(auto i = 0; i < n * n; i++){
      A[i] = 0.0;
    }

    T m = ((T)n - 1.0) * 0.5;
    for (auto i = 0; i < n; ++i) {
      T a1 = T(n - 1 - i); T a2 = 0. - i;
      int v = std::min(abs(a1), abs(a2));
      A[i + n * i] = m - v;
      if (i != n - 1) {
        A[i + 1 + n * i] = 1.0;
        A[i + n * (1 + i)] = 1.0;
      }
    }

    if(n <= 10){
        std::cout << "Show Final Matrix A: "<< std::endl;
        showMatrix<T>(A, n, n, n);
    }

    return A;
}

