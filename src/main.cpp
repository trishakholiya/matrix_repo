#include <armadillo>
#include "benchmark.hpp"
#include "Matrix.h"
#include <vector>
#include "naive_matrix.hpp"

#include <Eigen/Dense>

void print_vec(const std::vector<double>& v, const std::string& name) {
    std::cout << name << ": ";
    for (double x : v) std::cout << x << " ";
    std::cout << "\n";
}

Matrix identity_matrix(int n) {
    std::vector<double> flat(n * n, 0.0);
    for (int i = 0; i < n; ++i) flat[i*n + i] = 1.0;
    return Matrix(flat, n, n);
}

Matrix build_diagonal(const std::vector<double>& vals) {
    int n = static_cast<int>(vals.size());
    std::vector<double> flat(n * n, 0.0);
    for (int i = 0; i < n; ++i) {
        flat[i*n + i] = vals[i];
    }
    return Matrix(flat, n, n);
}

double frob_norm(const Matrix& A, const Matrix& B) {
    if (A.get_num_rows() != B.get_num_rows() || A.get_num_cols() != B.get_num_cols()) {
        throw std::runtime_error("Size mismatch in frob_norm");
    }
    int n = A.get_num_rows();
    int m = A.get_num_cols();
    double sumsq = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            double diff = A(i,j) - B(i,j);
            sumsq += diff * diff;
        }
    }
    return std::sqrt(sumsq);
}

Matrix make_symmetric(int n) {
    std::vector<double> flat(n * n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            // Any deterministic value works, just make it symmetric
            double val = ((i + 1) * (j + 1)) / double(n + 1);
            flat[i*n + j] = val;
            flat[j*n + i] = val;   // enforce symmetry
        }
    }

    return Matrix(flat, n, n);
}


// Diagonal test matrix with distinct entries
Matrix make_diagonal(int n) {
    std::vector<double> flat(n * n, 0.0);
    for (int i = 0; i < n; ++i) {
        flat[i*n + i] = i + 1.0;  // 1,2,3,...
    }
    return Matrix(flat, n, n);
}

// ---------- Main eigsym test ----------

void test_eigsym_on(const Matrix& A, const std::string& label) {
    std::cout << "=== Testing " << label << " ===\n";
    int n = A.get_num_rows();
    std::cout << "size n = " << n << "\n";

    EigsymResult eig = A.eigsym();

    Matrix P = eig.eigenvectors;
    Matrix D = build_diagonal(eig.eigenvalues);
    Matrix Pt = P.transpose();

    // Reconstruction: P D P^T
    Matrix A_rec = P * D * Pt;
    double recon_err = frob_norm(A_rec, A);

    // Orthogonality: P^T P ≈ I
    Matrix I = identity_matrix(n);
    double orth_err = frob_norm(Pt * P, I);

    // Direct diagonalization check: P^T A P should ≈ D
    Matrix D_alt = Pt * A * P;
    double diag_err = frob_norm(D_alt, D);

    std::cout << "reconstruction error ||P D P^T - A||_F = " << recon_err << "\n";
    std::cout << "orthogonality error ||P^T P - I||_F = " << orth_err << "\n";
    std::cout << "diag error ||P^T A P - D||_F       = " << diag_err << "\n\n";
}

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

  //multiply matrixes
  Matrix mat5 = mat2 * mat3;
  std::cout << "Matrix 5: " << mat5 << std::endl;

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


  std::vector<int> sizes = {2, 3, 5, 10};

    for (int n : sizes) {
        Matrix I = identity_matrix(n);
        Matrix D = make_diagonal(n);
        Matrix S = make_symmetric(n);

        test_eigsym_on(I, "Identity");
        test_eigsym_on(D, "Diagonal");
        test_eigsym_on(S, "Random-like symmetric");
    }

  

  return 0;
}