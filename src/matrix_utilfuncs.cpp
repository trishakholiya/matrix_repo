#include "matrix.h"
#include <cmath>

Matrix Matrix::diagmat(vec vector) {
  Matrix result = Matrix::Zeros(vector.size(), vector.size());
  for (int i = 0; i < vector.size(); i++ ) {
    result(i, i) = vector[i];
  }
  return result;
}


Matrix Matrix::diagmat(Matrix mat) {
  Matrix result = Matrix::Zeros(mat.get_num_rows(), mat.get_num_cols());
  for (int i = 0; i < std::min(mat.get_num_rows(), mat.get_num_cols()); i++) {
    result(i, i) = mat(i, i);
  }
  return result;
}

bool Matrix::is_symmetric(double tol = 1e-12) const {
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
Matrix Matrix::transpose() const {
  Matrix result(num_cols, num_rows);

  for (int i = 0; i < num_rows; i++) {
      for (int j = 0; j < num_cols; j++) {
          result(j, i) = (*this)(i, j);
      }
  }

  return result;
}