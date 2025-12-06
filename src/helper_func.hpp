#pragma once
#include <cmath>

// helper function
double SIGN(double a, double b) {
  if ( b >= 0.0 ) {
    return std::fabs(a);
  } else {
    return -std::fabs(a);
  }
}

// helper function
double pythag(const double a, const double b) {
  // computes sqrt(a^2 + b^2) without destructive underflow or overflow
  double absa = std::abs(a);
  double absb = std::abs(b);
  if (absa > absb) {
    return absa * sqrt(1.0 + (absb/absa)*(absb/absa));
  } else {
      return absb * sqrt(1.0 + (absa/absb)*(absa/absb));
  }
}