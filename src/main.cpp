#include <armadillo>
#include "benchmark.cpp"
#include "matrix.cpp"

int main(int argc, char* argv[]) {

  // initialize zero matrix
  Matrix mat1(2, 2);
  mat1.print();

  // initialize matrix with values
  Matrix mat2({1, 2, 3, 4}, 2, 2);
  mat2.print();

  // initialize matrix with values
  Matrix mat3({5, 6, 7, 8}, 2, 2);
  mat3.print();

  // add matrixes
  Matrix mat4 = mat2 + mat3;
  mat4.print();

  // benchmark the time
  auto matrix_sum_func = [&]() {
    Matrix temp = mat2 + mat3;
  };

  // initialize matrixes with same values in armadillo
  arma::mat A = { {1, 2},
                  {3, 4} };

  arma::mat B = { {5, 6},
                  {7, 8} };

  auto arma_sum_func = [&]() {
    arma::mat temp = A + B;
  };

  double matrix_time = benchmark_time(matrix_sum_func);
  double arma_time = benchmark_time(arma_sum_func);

  std::cout << "Matrix Sum Time: " << matrix_time << '\n';
  std::cout << "Armadillo Sum Time: " << arma_time << '\n';

  return 0;
}