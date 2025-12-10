#include <benchmark/benchmark.h>
#include "matrix.h"
#include <armadillo>

static void Multiplication_MatrixClass(benchmark::State& state) {
    int n = state.range(0);
    Matrix A = Matrix::Random(n, n);
    Matrix B = Matrix::Random(n, n);

    for (auto _ : state) {
        benchmark::DoNotOptimize(A * B);
    }

    state.SetItemsProcessed(state.iterations() * n * n);
}

static void Multiplication_Armadillo(benchmark::State& state) {
    int n = state.range(0);
    arma::mat A = arma::randu<arma::mat>(n, n);
    arma::mat B = arma::randu<arma::mat>(n, n);

    for (auto _ : state) {
      arma::mat C = A * B;
      benchmark::DoNotOptimize(C.memptr()); // needs to be pointer to memory NOT arma object
    }

    state.SetItemsProcessed(state.iterations() * n * n);
}

// Register benchmarks for different sizes
BENCHMARK(Multiplication_MatrixClass)
  ->Arg(10)
  ->Arg(100)
  ->Arg(200)
  ->Arg(400);

BENCHMARK(Multiplication_Armadillo)
  ->Arg(10)
  ->Arg(100)
  ->Arg(200)
  ->Arg(400);
