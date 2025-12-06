#pragma once

#include <iostream>
#include <vector>
#include "custom_exception.hpp"
#include "helper_func.hpp" // numerical recipes helper functions

typedef std::vector<double> vec;
typedef std::vector<vec> mat;

// Forward declarations of global result types
struct TridiagonalResult;
struct QREigenResult;
struct EigsymResult;
struct QLEigenResult;

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

  TridiagonalResult householder_tridiagonalize(bool yesvecs = true) const;
  QREigenResult QL(std::vector<double> d, std::vector<double> e) const;
  EigsymResult eigsym() const;

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
    std::vector<double> d; // diagonal elements
    std::vector<double> e; // off-diagonal
    Matrix Q_house; // accumulated Householder transforms
};

struct QLEigenResult {
    std::vector<double> eigenvalues;
    Matrix Q_ql; // accumulated QR transforms
};

struct EigsymResult {
    std::vector<double> eigenvalues;
    Matrix eigenvectors;   // columns = eigenvectors
};

inline TridiagonalResult Matrix::householder_tridiagonalize(bool yesvecs) const {
    TridiagonalResult result;
    // juliet implement
    int l, k, j, i;
    int n = num_rows;
    mat z = matrix;
    vec e(n), d(n);
    double scale, hh, h, g, f;
    for (i = n - 1; i > 0; i--) {
      l = i - 1;
      h = scale = 0.0;
      if (l > 0) {
        for (k = 0; k < i; k++)
          scale += abs(z[i][k]);
        if (scale == 0.0)
          e[i] = z[i][l];
        else {
          for (k = 0; k < i; k++) {
            z[i][k] /= scale;
            h += z[i][k] * z[i][k];
          }
        }
        f = z[i][l];
        g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
        e[i] = scale * g;
        h -= f * g;
        z[i][l] = f - g;
        f = 0.0;
        for (j = 0; j < i; j ++) {
          if (yesvecs)
            z[j][i] = z[i][j] / h;
          g = 0.0;
          for (k = 0; k < j; k++)
            g += z[j][k] * z[i][k];
          for (k = j + 1; k < i; k++)
            g += z[k][j] * z[i][k];
          e[j] = g / h;
          f += e[j] * z[i][j];
        }
        hh = f / (h + h);
        for (j = 0; j < i; j++) {
          f = z[i][j];
          e[j] = g = e[j] - hh * f;
          for (k = 0; k < j + 1; k++)
            z[j][k] -= (f * e[k] + g * z[i][k]);
        }
      } else
          e[i] = z[i][l];
        d[i] = h;
    }
    if (yesvecs)
      d[0] = 0.0;
    e[0] = 0.0;
    for (i = 0; i < n; i++) {
      if (yesvecs) {
        if (d[i] != 0.0) {
          for (j = 0; j < i; j++) {
            g = 0.0;
            for (k = 0; k < i; k++) 
              g += z[i][k] * z[k][j];
            for (k = 0; k < i; k++)
              z[k][j] -= g * z[k][i];
          }
        }
        
        d[i] = z[i][i];
        z[i][i] = 1.0;
        for (j = 0; j < i; j++)
          z[j][i] = z[i][j] = 0.0;

      } else {
        d[i] = z[i][i];
      }        
    }

    if (yesvecs) {
      // flatten z for now to construct matrix
      vec flatZ;
      flatZ.reserve(n * n);

        for (int r = 0; r < n; r++) {
          for (int c = 0; c < n; c++) {
            flatZ.push_back(z[r][c]);
          }
        }
      result.Q_house = Matrix(flatZ, n, n);
    }

    result.d = d;
    result.e = e;

    return result;
}

inline QLEigenResult Matrix::QL(std::vector<double> d, std::vector<double> e) const {
  QLEigenResult result;
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
  result.Q_ql = Matrix(flatZ, n, n);
  return result;
}

inline EigsymResult Matrix::eigsym() const {
    EigsymResult result;
    
    TridiagonalResult tri = householder_tridiagonalize(true);

    QREigenResult qr = QL(tri.d, tri.e);

    Matrix P = tri.Q_house * qr.Q_qr;

    result.eigenvalues = qr.eigenvalues;
    result.eigenvectors = P;

    return result;
}


