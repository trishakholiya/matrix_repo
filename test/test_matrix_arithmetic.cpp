// test/test_matrix_arithmetic.cpp
#include <gtest/gtest.h>
#include "matrix.h"
#include "test_helpers.hpp"

// file includes tests for add, subract, multiply, multiply by scalar

// creates 2 random 3 x 3 matrices using Matrix lib, adds them.
// converts same matrices to armadillo and adds, then compares the results.
TEST(MatrixArithmetic, AddMatchesArmadillo) {
    Matrix A = Matrix::Random(3,3);
    Matrix B = Matrix::Random(3,3);

    Matrix C_my = A + B;

    arma::mat A_ref = to_arma(A);
    arma::mat B_ref = to_arma(B);
    arma::mat C_ref = A_ref + B_ref;

    EXPECT_TRUE(mats_close(C_my, C_ref));
}

// same but for multiplication of 2 matrices
TEST(MatrixArithmetic, MulMatchesArmadillo) {
    Matrix A = Matrix::Random(2,3);
    Matrix B = Matrix::Random(3,4);

    Matrix C_my = A * B;

    arma::mat A_ref = to_arma(A);
    arma::mat B_ref = to_arma(B);
    arma::mat C_ref = A_ref * B_ref;

    EXPECT_TRUE(mats_close(C_my, C_ref));
}

// same but for multiplying matrix by a scalar
TEST(MatrixArithmetic, ScalarMulMatchesArmadillo) {
    Matrix A = Matrix::Random(3,3);
    double s = -2.5;

    Matrix C_my = A * s;
    arma::mat C_ref = to_arma(A) * s;

    EXPECT_TRUE(mats_close(C_my, C_ref));
}

/*
TEST(MatrixArithmetic, AddMatchesArmadilloMultipleSizes) {
    for (int n = 1; n <= 10; ++n) {
        Matrix A = Matrix::Random(n, n);
        Matrix B = Matrix::Random(n, n);

        Matrix C_my = A + B;

        arma::mat A_ref = to_arma(A);
        arma::mat B_ref = to_arma(B);
        arma::mat C_ref = A_ref + B_ref;

        EXPECT_TRUE(mats_close(C_my, C_ref))
            << "Addition failed for size " << n << "x" << n;
    }
}
*/