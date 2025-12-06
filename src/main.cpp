#include <armadillo>
#include "benchmark.hpp"
#include "Matrix.hpp"
#include <vector>
#include "naive_matrix.hpp"

#include <Eigen/Dense>


int main(int argc, char* argv[]) {

  // initialize zero matrix
  Matrix mat1(2, 2);
  std::cout << "Matrix 1: " << mat1 << std::endl;

  // initialize matrix with values
  Matrix mat2({1, 2, 3, 4}, 2, 2);
  std::cout << "Matrix 2: " << mat2 << std::endl;

  // initialize matrix with values
  Matrix mat3({5, 6, 7, 8}, 2, 2);
  std::cout << "Matrix 3: " << mat3 << std::endl;

  // add matrixes
  Matrix mat4 = mat2 + mat3;
  std::cout << "Matrix 4: " << mat4 << std::endl;

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

  Eigen::Matrix<double, 2, 2> E1;
    E1 << 1, 2,
          3, 4;

    Eigen::Matrix<double, 2, 2> E2;
    E2 << 5, 6,
          7, 8;

    auto eigen_sum_func = [&]() {
        auto temp = E1 + E2;
    };

  naive_matrix N1 = make_naive_matrix(2, 2, {1,2,3,4});
  naive_matrix N2 = make_naive_matrix(2, 2, {5,6,7,8});

  auto naive_sum_func = [&]() {
      auto temp = N1 + N2;
  };
  double naive_time = benchmark_time(naive_sum_func);
  double eigen_time = benchmark_time(eigen_sum_func);
  double matrix_time = benchmark_time(matrix_sum_func);
  double arma_time = benchmark_time(arma_sum_func);

  std::cout << "Naive C++ Implementation Sum Time: " << naive_time << '\n';
  std::cout << "Eigen Sum Time: " << eigen_time << '\n';
  std::cout << "Matrix Sum Time: " << matrix_time << '\n';
  std::cout << "Armadillo Sum Time: " << arma_time << '\n';


  return 0;
}