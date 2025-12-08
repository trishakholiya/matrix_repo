#include "matrix.h"
#include <cmath>

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