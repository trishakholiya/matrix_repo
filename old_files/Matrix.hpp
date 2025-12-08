/*
#pragma once

#include <iostream>
#include <vector>
#include <limits>
#include "custom_exception.hpp"
#include "helper_func.hpp" // numerical recipes helper functions
#include <cmath>

typedef std::vector<double> vec;
typedef std::vector<vec> mat;

// Forward declarations of global result types
struct TridiagonalResult;
// struct QREigenResult;
struct EigsymResult;
struct QLEigenResult;

class Matrix {
private:
  int num_rows;
  int num_cols;
  int size;
  vec matrix;

public:
  // initialize matrix of specific size
  Matrix(int rows, int cols) : num_rows(rows), num_cols(cols), size(rows*cols), matrix(size) {}

  // empty constructor
  Matrix() : num_rows(0), num_cols(0), size(0) {}

  // constructor from FLAT vector
  Matrix(const vec& values, int rows, int cols) // l-value
    : num_rows(rows), 
      num_cols(cols), 
      size(rows*cols)
  {
    // make sure this matrix can be constructed
    if (values.size() != rows * cols) {
        throw InvalidMatrixSize("Flat vector size does not match requested matrix dimensions");
    }

    matrix = values;
  }

  Matrix(vec&& values, int rows, int cols) // r-value
    : num_rows(rows), 
      num_cols(cols), 
      size(rows*cols),
      matrix(std::move(values)) {
      // make sure this matrix can be constructed
      if (values.size() != rows * cols) {
        throw InvalidMatrixSize("Flat vector size does not match requested matrix dimensions");
      }
    }

  static Matrix Ones(int rows, int cols) {
    Matrix M(rows, cols);
    std::fill(M.matrix.begin(), M.matrix.end(), 1.0);
    return M;
  }

  static Matrix Zeros(int rows, int cols) {
    Matrix M(rows, cols);
    std::fill(M.matrix.begin(), M.matrix.end(), 0.0);
    return M;
  }

  static Matrix Identity(int n) {
    Matrix M(n, n);
    std::fill(M.matrix.begin(), M.matrix.end(), 0.0);
    for (int i = 0; i < n; i++)
        M(i, i) = 1.0;
    return M;
  }

  // Accessors (getters)
  int get_num_rows() const
  {
    return this->num_rows;
  }

  int get_num_cols() const
  {
    return this->num_cols;
  }

  int get_size() const
  {
    return this->size;
  }

  const vec& get_data() const
  {
    return this->matrix;
  }

  // operator to access the matrix so you can call matrix(1, 1)
  inline double& operator()(int x, int y) {
    if (x < 0 || x >= num_rows || y < 0 || y >= num_cols) {
      throw std::out_of_range("Matrix index out of range");
    }
    return matrix[x * num_cols + y];
  }

  // operator to access the matrix so you can call matrix(1, 1)
  inline const double& operator()(int x, int y) const {
    if (x < 0 || x >= num_rows || y < 0 || y >= num_cols) {
      throw std::out_of_range("Matrix index out of range");
    }
    return matrix[x * num_cols + y];
  }

  TridiagonalResult householder_tridiagonalize(bool yesvecs = true) const;
  QLEigenResult QL(std::vector<double> d, std::vector<double> e) const;
  EigsymResult eigsym() const;

  // equal to a matrix operator
  inline bool operator==(const Matrix& other) const {
    if (num_rows != other.num_rows || num_cols != other.num_cols) {
      return false;
    }

    const double* m1 = this->matrix.data();
    const double* m2 = other.matrix.data();

    double epsilon = std::numeric_limits<double>::epsilon(); // NEED TO SEE IF THIS IS TOO STRICT?

    for ( int i = 0; i < this->size ; i++ ) {
      if ( std::abs(m1[i] - m2[i]) > epsilon ) return false;
    }

    return true;
  }

  // overload addition operator
  inline Matrix operator+(const Matrix& other) const {

    if (num_rows != other.num_rows || num_cols != other.num_cols) {
        throw InvalidMatrixSize("Matrix sizes must match for addition");
    }

    vec result(size);

    const double* m1 = this->matrix.data();
    const double* m2 = other.matrix.data();
    double* r = result.data();

    for ( int i = 0; i < this->size ; i++ ) {
      r[i] = m1[i] + m2[i];
    }

    return Matrix(result, num_rows, num_cols);
  }

  // overload subtraction operator
  inline Matrix operator-(const Matrix& other) const {
    vec result(size);

    const double* m1 = this->matrix.data();
    const double* m2 = other.matrix.data();
    double* r = result.data();

    for ( int i = 0; i < this->size ; i++ ) {
      r[i] = m1[i] - m2[i];
    }

    return Matrix(result, num_rows, num_cols);
  }

  inline Matrix operator*(const Matrix& other) const {
    if (this->num_cols != other.num_rows) {
        throw InvalidMatrixSize("Matrix dimensions incompatible for multiplication");
    }

    Matrix result::Zeros(this->num_rows, other.num_cols);

    for (int i = 0; i < this->num_rows; ++i) {
        for (int j = 0; j < other.num_cols; ++j) {
            double sum = 0.0;
            for (int k = 0; k < this->num_cols; ++k) {
                sum += (*this)(i, k) * other(k, j);
            }
            result(i, j) = sum;
        }
    }

    return result;
  }

  Matrix operator*(double s) const {
    vec result_vec(get_size());

    std::transform(matrix.begin(), matrix.end(), result_vec.begin(), 
                   [&s](double val) { return val * s; });
    return Matrix(result_vec, get_num_rows(), get_num_cols());
  }

  bool is_symmetric(double tol = 1e-12) const {
    if (num_rows != num_cols)
        return false;

    for (int i = 0; i < num_rows; i++) {
        for (int j = i + 1; j < num_cols; j++) {  
            if (std::abs((*this)(i, j) - (*this)(j, i)) > tol)
                return false;
        }
    }
    return true;
}

  // transpose matrix
  Matrix transpose() const {
    Matrix result(num_cols, num_rows);

    for (int i = 0; i < num_rows; i++) {
        for (int j = 0; j < num_cols; j++) {
            result(j, i) = (*this)(i, j);
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
    Matrix Q_ql; // accumulated QL transforms
};

struct EigsymResult {
    vec eigenvalues;
    Matrix eigenvectors;   // columns = eigenvectors
};

inline TridiagonalResult Matrix::householder_tridiagonalize(bool yesvecs) const {
    TridiagonalResult result;
    // juliet implement
    int l, k, j, i;
    int n = num_rows;
    Matrix z = *this;
    vec e(n), d(n);
    double scale, hh, h, g, f;
    for (i = n - 1; i > 0; i--) {
      l = i - 1;
      h = scale = 0.0;
      if (l > 0) {
        for (k = 0; k < i; k++)
          scale += std::abs(z(i, k));
        if (scale == 0.0) {
          e[i] = z(i, l);
        } else {
          for (k = 0; k < i; k++) {
            z(i, k) /= scale;
            h += z(i, k) * z(i, k);
          }

          f = z(i, l);
          g = (f >= 0.0 ? -std::sqrt(h) : std::sqrt(h));
          e[i] = scale * g;
          h -= f * g;
          z(i, l) = f - g;
          f = 0.0;
          for (j = 0; j < i; j ++) {
            if (yesvecs)
              z(j, i) = z(i, j) / h;
            g = 0.0;
            for (k = 0; k < j+1; k++)
              g += z(j, k) * z(i, k);
            for (k = j + 1; k < i; k++)
              g += z(k, j) * z(i, k);
            e[j] = g / h;
            f += e[j] * z(i, j);
        }

        hh = f / (h + h);

        for (j = 0; j < i; j++) {
          f = z(i, j);
          e[j] = g = e[j] - hh * f;

          for (k = 0; k < j + 1; k++)
            z(j, k) -= (f * e[k] + g * z(i, k));
        }
      }
      } else
          e[i] = z(i, l);
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
              g += z(i, k) * z(k, j);
            for (k = 0; k < i; k++)
              z(k, j) -= g * z(k, i);
          }
        }
        
        d[i] = z(i, i);
        z(i, i) = 1.0;
        for (j = 0; j < i; j++)
          z(j, i) = z(i, j) = 0.0;

      } else {
        d[i] = z(i, i);
      }        
    }

    if (yesvecs) {
      result.Q_house = z;
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

  Matrix z = Matrix::Identity(n); // to store eigenvectors

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
            f = z(k, i+1);
            z(k, i+1) = s*z(k, i)+c*f;
            z(k, i) = c*z(k, i)-s*f;
          }

        }
        if ( r== 0.0 && i >= l) continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while ( m != l );
  }

  result.eigenvalues = d;
  result.Q_ql = z;
  return result;
}

inline EigsymResult Matrix::eigsym() const {
  if (num_rows != num_cols) {
    throw InvalidMatrixSize("householder_tridiagonalize requires a square matrix");
  }

  if (!is_symmetric()) {
    throw InvalidMatrixSize("Matrix must be symmetric for eigsym()");
  }
  
  EigsymResult result;
  
  TridiagonalResult tri = householder_tridiagonalize(true);

  QLEigenResult ql = QL(tri.d, tri.e);

  Matrix P = tri.Q_house * ql.Q_ql;

  result.eigenvalues = ql.eigenvalues;
  result.eigenvectors = P;

  return result;
}

// overload print << operator
std::ostream& operator<<(std::ostream& out, const Matrix & M) {
  out << std::fixed << std::setprecision(6); // set precision
  out << "\n";
  for (int i = 0; i < M.get_num_rows(); i++) {
      for (int j = 0; j < M.get_num_cols(); j++) {
          out << set::sw(10) << M(i, j) << "  "; // set a constant width
      }
      out << "\n"; // line between rows
  }
  out << "\n";
  return out;
}
*/