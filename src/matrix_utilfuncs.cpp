#include "matrix.h"
#include <cmath>
#include <highfive/H5File.hpp>


Matrix Matrix::diagmat(const vec& vector) {
  Matrix result = Matrix::Zeros(vector.size(), vector.size());
  for (int i = 0; i < vector.size(); i++ ) {
    result(i, i) = vector[i];
  }
  return result;
}


Matrix Matrix::diagmat(const Matrix& mat) {
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

void Matrix::save_hdf5(const Matrix& data, 
                        const std::string& filename,
                       const std::string& dataset_name)
{
  int rows = data.get_num_rows();
  int cols = data.get_num_cols();

  // get file, might need to create it
  HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create);

  // need to use hsize_t type from HDF5
  std::vector<hsize_t> dimensions = { 
      static_cast<hsize_t>(rows),
      static_cast<hsize_t>(cols)
  };

  // Convert flat to 2D
  std::vector<vec> reshaped(rows, vec(cols));
    const vec& flat = data.get_data();
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            reshaped[i][j] = flat[i * cols + j];
        }
    }

  // Create dataset
  HighFive::DataSet dataset = file.createDataSet<double>(
      dataset_name,
      HighFive::DataSpace({rows, cols}));

    // add matrix data to the dataset
  dataset.write(reshaped);
}

void Matrix::save_hdf5(const vec& data, 
                        const std::string& filename,
                        const std::string& dataset_name)
{
  // get file, might need to create it
  HighFive::File file(filename, HighFive::File::ReadWrite | HighFive::File::Create);

  // need to use hsize_t type from HDF5
  std::vector<hsize_t> dimensions = {
    static_cast<hsize_t>(data.size())
    };

  // Create dataset
  HighFive::DataSet dataset = file.createDataSet<double>(
    dataset_name,
    HighFive::DataSpace(dimensions)
  );

  // add vec data to the dataset
  dataset.write(data);
}