#include <iostream>
#include <random>
#include <cmath>
#include "Matrix.hpp"

// Build a diagonal Matrix from eigenvalues
Matrix make_diag(const std::vector<double>& eigvals) {
    int n = static_cast<int>(eigvals.size());
    std::vector<double> flat(n * n, 0.0);
    for (int i = 0; i < n; ++i) {
        flat[i * n + i] = eigvals[i];
    }
    return Matrix(flat, n, n);
}

// Make a random symmetric matrix of size n
Matrix make_random_symmetric(int n, std::mt19937& rng, double scale = 1.0) {
    std::uniform_real_distribution<double> dist(-scale, scale);
    std::vector<double> flat(n * n, 0.0);

    // fill upper triangle and mirror
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double v = dist(rng);
            flat[i * n + j] = v;
            flat[j * n + i] = v;
        }
    }
    return Matrix(flat, n, n);
}

// Frobenius norm of a Matrix
double frobenius_norm(const Matrix& A) {
    int n = A.rows();
    int m = A.cols();
    double sumsq = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            double v = A.at(i, j);
            sumsq += v * v;
        }
    }
    return std::sqrt(sumsq);
}

// Check orthogonality: ||P^T P - I||_F
double orthogonality_error(const Matrix& P) {
    int n = P.rows();
    int m = P.cols(); // should be same as n for full eigensystem

    // build identity I
    std::vector<double> flatI(m * m, 0.0);
    for (int i = 0; i < m; ++i) flatI[i * m + i] = 1.0;
    Matrix I(flatI, m, m);

    Matrix Pt = P.transpose();
    Matrix PtP = Pt * P;
    Matrix diff = PtP - I;
    return frobenius_norm(diff);
}

// Test pipeline on a single matrix A
void test_matrix(const Matrix& A, const std::string& name) {
    std::cout << "=== Testing " << name << " ===\n";
    int n = A.rows();

    // Full eigensystem
    EigsymResult res = A.eigsym();
    const std::vector<double>& eigvals = res.eigenvalues;
    const Matrix& P = res.eigenvectors;

    // Build D and reconstruct A_hat = P D P^T
    Matrix D = make_diag(eigvals);
    Matrix A_hat = P * D * P.transpose();

    Matrix diff = A_hat - A;
    double recon_err = frobenius_norm(diff);
    double ortho_err = orthogonality_error(P);

    std::cout << "size n = " << n << "\n";
    std::cout << "reconstruction error ||P D P^T - A||_F = " << recon_err << "\n";
    std::cout << "orthogonality error ||P^T P - I||_F = " << ortho_err << "\n\n";
}

int main() {
    std::mt19937 rng(12345);

    int sizes[] = {2, 3, 5, 10};

    for (int n : sizes) {
        // 1) Identity
        std::vector<double> flatI(n * n, 0.0);
        for (int i = 0; i < n; ++i) flatI[i * n + i] = 1.0;
        Matrix I(flatI, n, n);
        test_matrix(I, "Identity");

        // 2) Simple diagonal with distinct entries
        std::vector<double> flatD(n * n, 0.0);
        for (int i = 0; i < n; ++i) flatD[i * n + i] = i + 1; // 1,2,3,...
        Matrix D(flatD, n, n);
        test_matrix(D, "Diagonal");

        // 3) Random symmetric
        Matrix R = make_random_symmetric(n, rng);
        test_matrix(R, "Random symmetric");
    }

    return 0;
}

