#include "matrix.h"
#include <algorithm>

Matrix::Matrix(int rows, int cols) : num_rows(rows), num_cols(cols), size(rows*cols), matrix(size) {}

// empty constructor
Matrix::Matrix() : num_rows(0), num_cols(0), size(0) {}

// constructor from FLAT vector
Matrix::Matrix(const vec& values, int rows, int cols) // l-value
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

Matrix::Matrix(vec&& values, int rows, int cols) // r-value
  : num_rows(rows), 
    num_cols(cols), 
    size(rows*cols),
    matrix(std::move(values)) {
    // make sure this matrix can be constructed
    if (matrix.size() != rows * cols) {
      throw InvalidMatrixSize("Flat vector size does not match requested matrix dimensions");
    }
  }

Matrix Matrix::Ones(int rows, int cols) {
  Matrix M(rows, cols);
  std::fill(M.matrix.begin(), M.matrix.end(), 1.0);
  return M;
}

Matrix Matrix::Zeros(int rows, int cols) {
  Matrix M(rows, cols);
  std::fill(M.matrix.begin(), M.matrix.end(), 0.0);
  return M;
}

Matrix Matrix::Identity(int n) {
  Matrix M(n, n);
  std::fill(M.matrix.begin(), M.matrix.end(), 0.0);
  for (int i = 0; i < n; i++)
      M(i, i) = 1.0;
  return M;
}

// Accessors (getters)
int Matrix::get_num_rows() const
{
  return this->num_rows;
}

int Matrix::get_num_cols() const
{
  return this->num_cols;
}

int Matrix::get_size() const
{
  return this->size;
}

const vec& Matrix::get_data() const
{
  return this->matrix;
}