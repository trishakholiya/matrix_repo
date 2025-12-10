#include "matrix.h"
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <limits>

// operator to access the matrix so you can call matrix(1, 1)
double& Matrix::operator()(int x, int y) {
  if (x < 0 || x >= num_rows || y < 0 || y >= num_cols) {
    throw std::out_of_range("Matrix index out of range");
  }
  return matrix[x * num_cols + y];
}

// operator to access the matrix so you can call matrix(1, 1)
const double& Matrix::operator()(int x, int y) const {
  if (x < 0 || x >= num_rows || y < 0 || y >= num_cols) {
    throw std::out_of_range("Matrix index out of range");
  }
  return matrix[x * num_cols + y];
}

// equal to a matrix operator
bool Matrix::operator==(const Matrix& other) const {
  if (num_rows != other.num_rows || num_cols != other.num_cols) {
    return false;
  }

  const double* m1 = this->matrix.data();
  const double* m2 = other.matrix.data();

  double epsilon = 1e-10;

  for ( int i = 0; i < this->size ; i++ ) {
    if ( std::abs(m1[i] - m2[i]) > epsilon ) return false;
  }

  return true;
}

// overload addition operator
Matrix Matrix::operator+(const Matrix& other) const {

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
Matrix Matrix::operator-(const Matrix& other) const {
  vec result(size);

  const double* m1 = this->matrix.data();
  const double* m2 = other.matrix.data();
  double* r = result.data();

  for ( int i = 0; i < this->size ; i++ ) {
    r[i] = m1[i] - m2[i];
  }

  return Matrix(result, num_rows, num_cols);
}

Matrix Matrix::operator*(const Matrix& other) const {
  if (this->num_cols != other.num_rows) {
      throw InvalidMatrixSize("Matrix dimensions incompatible for multiplication");
  }
  int rows = this->num_rows;
  int cols = other.num_cols;
  int my_cols = this->num_cols;

  Matrix result(rows, cols);

  const double* m1 = this->matrix.data();
  const double* m2 = other.matrix.data();
  double* r = result.matrix.data();

  for (int i = 0; i < rows; ++i) {
    for (int k = 0; k < my_cols; ++k) {
        double a = m1[i * my_cols + k];
        for (int j = 0; j < cols; ++j) {
            r[i * cols + j] += a * m2[k * cols + j];
        }
    }
  }
  return result;
}

Matrix Matrix::operator*(double s) const {
  vec result_vec(get_size());

  std::transform(matrix.begin(), matrix.end(), result_vec.begin(), 
                  [&s](double val) { return val * s; });
  return Matrix(result_vec, get_num_rows(), get_num_cols());
}

std::ostream& operator<<(std::ostream& out, const Matrix & M) {
  out << std::fixed << std::setprecision(4);
  out << "\n";
  for (int i = 0; i < M.get_num_rows(); i++) {
      for (int j = 0; j < M.get_num_cols(); j++) {
          out << std::setw(8) << M(i, j) << "  "; // set a constant width
      }
      out << "\n"; // line between rows
  }
  out << "\n";
  return out;
}