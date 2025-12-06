#include <iostream>
#include "Matrix.hpp"

void print_vec(const std::vector<double>& v, const std::string& name) {
    std::cout << name << ": ";
    for (double x : v) std::cout << x << " ";
    std::cout << "\n";
}

void test_small() {
    Matrix A({1,3,
              3,4}, 2, 2);
    TridiagonalResult tri = A.householder_tridiagonalize(true);

    print_vec(tri.d, "d");
    print_vec(tri.e, "e");

    std::cout << "Q_house:\n";
    tri.Q_house.print();
}

int main() {
    test_small();
    return 0;
}

