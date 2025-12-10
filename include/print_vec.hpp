#pragma once
#include <iostream>
#include <vector>
#include <iomanip>

inline std::ostream& operator<<(std::ostream& out, const std::vector<double>& v) {
  out << std::fixed << std::setprecision(4);
  for (double x : v) {
    out << x << "\n";
  }
  return out;
}