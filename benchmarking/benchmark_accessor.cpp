#include <benchmark/benchmark.h>
#include "matrix.h"
#include <armadillo>

static void Accessor_MatrixClass(benchmark::State& state) {
  int n = state.range(0);

  Matrix A = Matrix::Zeros(n, n);

  for (auto _ : state) {
    double sum = 0.0;
    for (int i = 0; i < n; i++ ) {
      for ( int j = 0; j < n; j++ ) {
        sum += A(i, j);
      }
    }
    benchmark::DoNotOptimize(sum);
  }

  state.SetItemsProcessed(state.iterations() * n * n);
}

static void Accessor_Armadillo(benchmark::State& state) {
  int n = state.range(0);
  arma::mat A = arma::zeros(n, n);

  for (auto _ : state) {
    double sum = 0.0;
    for (int i = 0; i < n; i++ ) {
      for ( int j = 0; j < n; j++ ) {
        sum += A(i, j);
      }
    }
    benchmark::DoNotOptimize(sum);
  }

  state.SetItemsProcessed(state.iterations() * n * n);
}


// benchmark for differenct sizes
BENCHMARK(Accessor_MatrixClass)
  ->Arg(10)
  ->Arg(100)
  ->Arg(200)
  ->Arg(400);

BENCHMARK(Accessor_Armadillo)
  ->Arg(10)
  ->Arg(100)
  ->Arg(200)
  ->Arg(400);
