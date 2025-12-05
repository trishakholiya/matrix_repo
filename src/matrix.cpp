#pragma once

#include <iostream>
#include <vector>
#include "custom_exception.cpp"
#include "helper_func.cpp" // numerical recipes helper functions

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
  QREigenResult QL(std::vector<double> d, std::vector<double> e);

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

inline QREigenResult Matrix::QL(std::vector<double> d, std::vector<double> e) {
  QREigenResult result;
  int n = d.size();
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;
  const double eps=std::numeric_limits<double>::epsilon();


  // initialize Z matrix (to store eigenvectors)
  std::vector<std::vector<double>> z(n, std::vector<double>(n, 0.0));
  for (i = 0; i < n; i++) {
    z[i][i] = 1.0;
  }

  for (i = 1; i < n; i++) {
    e[i-1] = e[i];
  }
  e[n-1] = 0.0;

  for ( l=0; l<n; l++ ) {
    iter = 0;
    do {
      for ( m=l; m<n-1; m++ ) {
        dd = std::abs(d[m]) + std::abs(d[m+1]);
        if ( std::abs(e[m]) <= eps*dd ) break;
      }
      if ( m!=l ) {
        if ( iter++ == 30 ) throw std::runtime_error("Too many iterations in tqli");
        g = (d[l+1] - d[l]) / (2.0 * e[l]);
        r = pythag(g, 1.0);
        g = d[m] - d[l] + e[l]/(g + SIGN(r,g));
        s = 1.0;
        c = 1.0;
        p = 0.0;

        for (i = m-1; i >= l; i--) {
          f = s * e[i];
          b = c * e[i];

          e[i+1] = (r=pythag(f, g));
          if ( r == 0 ) {
            d[i+1] -= p;
            e[m] = 0.0;
            break;
          }
          s = f / r;
          c = g / r;
          g = d[i+1]-p;
          r =  (d[i] - g)*s+2.0*c*b;
          d[i+1]  = g + (p=s*r);
          g = c*r - b;
          
          // FORM EIGENVECTORS
          for ( k=0; k<n; k++ ) {
            f = z[k][i+1];
            z[k][i+1] = s*z[k][i]+c*f;
            z[k][i] = c*z[k][i]-s*f;
          }

        }
        if ( r== 0.0 && i >= l) continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while ( m != l );
  }

  // flatten z for now to construct matrix
  vec flatZ;
  flatZ.reserve(n * n);

  for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
          flatZ.push_back(z[i][j]);
      }
  }
  result.eigenvalues = d;
  result.Q_qr = Matrix(flatZ, n, n);
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

