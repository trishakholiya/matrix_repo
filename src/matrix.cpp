#pragma once

#include <iostream>
#include <vector>
#include "custom_exception.cpp"

typedef std::vector<double> vec;
typedef std::vector<vec> mat;

// Forward declarations of global result types
struct TridiagonalResult;
struct QREigenResult;

class Matrix {
private:
  int num_rows;
  int num_cols;
  int size;
  mat matrix;
  vec data; // store flattened version of matrix for speed, MAYBE SHOULD CONSIDER ONLY STORING FLATTENED MATRIX

public:
  // initialize matrix of specific size with zeros
  Matrix(int rows, int cols) : num_rows(rows), num_cols(cols), size(0), matrix(rows, vec(cols, 0.0)) {}

  // empty constructor
  Matrix() : num_rows(0), num_cols(0), size(0) {}

  // constructor from FLAT vector
  Matrix(const vec& values, int rows, int cols)
    : num_rows(rows), num_cols(cols), size(num_rows*num_cols)
  {
    // make sure this matrix can be constructed
    if (values.size() != rows * cols) {
        throw InvalidMatrixSize("Flat vector size does not match requested matrix dimensions");
    }

    matrix.resize(rows, vec(cols));
    int idx = 0; // for indexing into vector
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i][j] = values[idx];
            idx++;
        }
    }
  }

  TridiagonalResult householder_tridiagonalize();
  QREigenResult qr_eigen_tridiagonal();


  // equal to a matrix operator
  // overload print << operator

  
  // overload addition operator
  inline Matrix operator+(const Matrix& other) const {
    // check if matrices are able to be added
    if (this->num_cols != other.num_cols || this->num_rows != other.num_rows) {
      throw InvalidMatrixSize("Matrix sizes must be equivalent for summation");
    }

    Matrix result(num_rows, num_cols);

    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            result.matrix[i][j] = matrix[i][j] + other.matrix[i][j];
        }
    }
      
    return result;
  }

  void print() const {
    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            std::cout << matrix[i][j] << "  "; // need to figure out spacing between numbers
        }
        std::cout << "\n"; // line between rows
    }
    std::cout << "\n";
  }

  double L2_norm() {
    // need to implement this for benchmarking
    return 0.0;
  }

  // transpose matrix
  Matrix transpose() const {
    Matrix result(num_cols, num_rows);

    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            result.matrix[j][i] = matrix[i][j];
        }
    }

    return result;
  }
};

struct TridiagonalResult {
    Matrix T; // tridiagonal matrix
    Matrix Q_house; // accumulated Householder transforms
};

struct QREigenResult {
    std::vector<double> eigenvalues;
    Matrix Q_qr; // accumulated QR transforms
};

inline TridiagonalResult Matrix::householder_tridiagonalize() {
    TridiagonalResult result;
    // juliet implement
    return result;
}

inline QREigenResult Matrix::qr_eigen_tridiagonal() {
    QREigenResult result;
    // trisha implement
    return result;
}

  /*
  std::pair<vec, mat> compute eigsym_decomposition() {
    // returns eigenvalues and eigenvectors at the same time 
    // so you don't have to recompute 3 times

    // implement
    int l, k, j, i;
    double scale, hh, h, g, f;
    for (i = n - 1; i > 0; i--) {
      l = i -1;
      h, scale = 0.0;
      if (l > 0) {
        for (k = 0; k < i; k++) {
          scale += abs(z[i][k]);
        }
        if (scale == 0) {
          e[i] = z[i][j];
        } else {
          for (k = 0; k < i; k++) {
            z[i][k] /= scale;
            h += z[i][k] * z[i][k];
          }
        }
      }
    }
  }

  vec get_eigenvalues() const {
    // implement getter
  }

  Matrix get_eigenvectors() const {
    // implement getter
  }

  // diagonalize matrix
  Matrix diagonalize() const {
    Matrix result(num_rows, num_cols);
    // tridiagonal 
    // A = PDP^-1
    // compute P (eigenvectors matrix of A) and check if it has full rank
    // if full rank then it's diagonalizable, if not throw error (need to include a tolerance in the rank check)

    // implement

    return result;
  }
  */

