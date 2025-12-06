#pragma once
#include <vector>

using naive_matrix = std::vector<std::vector<double>>;

naive_matrix make_naive_matrix(int rows, int cols, const std::vector<double>& vals) {
    naive_matrix M(rows, std::vector<double>(cols));
    int idx = 0;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            M[i][j] = vals[idx++];
    return M;
}

naive_matrix operator+(const naive_matrix& A, const naive_matrix& B) {
    int rows = A.size();
    int cols = A[0].size();

    naive_matrix C(rows, std::vector<double>(cols));

    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            C[i][j] = A[i][j] + B[i][j];

    return C;
}
