#pragma once

#include <iostream>
#include <vector>
#include <limits>
#include "custom_exception.hpp"
#include "helper_func.hpp" // numerical recipes helper functions
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <random>

typedef std::vector<double> vec;

// Forward declarations of global result types
struct TridiagonalResult;
struct EigsymResult;
struct QLEigenResult;

class Matrix {
private:
  int num_rows;
  int num_cols;
  int size;
  vec matrix;

public:
  // constructors
  Matrix(int rows, int cols);
  Matrix();
  Matrix(const vec& values, int rows, int cols);
  Matrix(vec&& values, int rows, int cols);

  static Matrix Ones(int rows, int cols);
  static Matrix Zeros(int rows, int cols);
  static Matrix Random(int rows, int cols);
  static Matrix Identity(int n);

  // Accessors (getters)
  int get_num_rows() const;
  int get_num_cols() const;
  int get_size() const;
  const vec& get_data() const;

  // operators
  double& operator()(int x, int y);
  const double& operator()(int x, int y) const;
  bool operator==(const Matrix& other) const;
  Matrix operator+(const Matrix& other) const;
  Matrix operator-(const Matrix& other) const;
  Matrix operator*(const Matrix& other) const;
  Matrix operator*(double s) const;

  // linear algebra functions
  static Matrix diagmat(const vec& vector);
  static Matrix diagmat(const Matrix& mat);
  bool is_symmetric(double tol) const;
  Matrix transpose() const;
  TridiagonalResult householder_tridiagonalize(bool yesvecs = true) const;
  QLEigenResult QL(vec d, vec e) const;
  EigsymResult eigsym() const;

  // saving
  static void save_hdf5(const Matrix& data, const std::string& filename, const std::string& dataset_name);
  static void save_hdf5(const vec& data, const std::string& filename, const std::string& dataset_name);

};

// structs for eigen decomposition
struct TridiagonalResult { 
    vec d;
    vec e;
    Matrix Q_house;
};

struct QLEigenResult {
    vec eigenvalues;
    Matrix Q_ql;
};

struct EigsymResult {
    vec eigenvalues;
    Matrix eigenvectors;
};

// printing matrix class
std::ostream& operator<<(std::ostream& out, const Matrix & M);