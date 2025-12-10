#include <benchmark/benchmark.h>
#include "matrix.h"
#include <armadillo>

static void EigSym_MatrixClass(benchmark::State& state) {
  int n = state.range(0);
  Matrix B = Matrix::Random(n, n);
  Matrix A = B + B.transpose(); // needs to be symmetric

  for (auto _ : state) {
    auto result = A.eigsym();
    benchmark::DoNotOptimize(result);
  }

  state.SetItemsProcessed(state.iterations() * n * n);
}

static void EigSym_Armadillo(benchmark::State& state) {
  int n = state.range(0);
  arma::mat B = arma::randu<arma::mat>(n, n);
  arma::mat A = B + B.t(); // needs to be symmetric

  for (auto _ : state) {
    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, A);
    benchmark::DoNotOptimize(eigval.memptr());
    benchmark::DoNotOptimize(eigvec.memptr());
  }

  state.SetItemsProcessed(state.iterations() * n * n);
}

// benchmark for differenct sizes
BENCHMARK(EigSym_MatrixClass)
  ->Arg(10)
  ->Arg(100)
  ->Arg(200)
  ->Arg(400);

BENCHMARK(EigSym_Armadillo)
  ->Arg(10)
  ->Arg(100)
  ->Arg(200)
  ->Arg(400);
