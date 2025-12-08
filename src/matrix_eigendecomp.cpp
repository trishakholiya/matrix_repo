#include "matrix.h"
#include "helper_func.cpp"
#include <cmath>
#include <stdexcept>
#include <limits>

TridiagonalResult Matrix::householder_tridiagonalize(bool yesvecs) const {
    TridiagonalResult result;
    int l, k, j, i;
    int n = num_rows;
    Matrix z = *this;
    vec e(n), d(n);
    double scale, hh, h, g, f;
    for (i = n - 1; i > 0; i--) {
      l = i - 1;
      h = scale = 0.0;
      if (l > 0) {
        for (k = 0; k < i; k++)
          scale += std::abs(z(i, k));
        if (scale == 0.0) {
          e[i] = z(i, l);
        } else {
          for (k = 0; k < i; k++) {
            z(i, k) /= scale;
            h += z(i, k) * z(i, k);
          }

          f = z(i, l);
          g = (f >= 0.0 ? -std::sqrt(h) : std::sqrt(h));
          e[i] = scale * g;
          h -= f * g;
          z(i, l) = f - g;
          f = 0.0;
          for (j = 0; j < i; j ++) {
            if (yesvecs)
              z(j, i) = z(i, j) / h;
            g = 0.0;
            for (k = 0; k < j+1; k++)
              g += z(j, k) * z(i, k);
            for (k = j + 1; k < i; k++)
              g += z(k, j) * z(i, k);
            e[j] = g / h;
            f += e[j] * z(i, j);
        }

        hh = f / (h + h);

        for (j = 0; j < i; j++) {
          f = z(i, j);
          e[j] = g = e[j] - hh * f;

          for (k = 0; k < j + 1; k++)
            z(j, k) -= (f * e[k] + g * z(i, k));
        }
      }
      } else
          e[i] = z(i, l);
        d[i] = h;
    }
    if (yesvecs)
      d[0] = 0.0;
    e[0] = 0.0;
    for (i = 0; i < n; i++) {
      if (yesvecs) {
        if (d[i] != 0.0) {
          for (j = 0; j < i; j++) {
            g = 0.0;
            for (k = 0; k < i; k++) 
              g += z(i, k) * z(k, j);
            for (k = 0; k < i; k++)
              z(k, j) -= g * z(k, i);
          }
        }
        
        d[i] = z(i, i);
        z(i, i) = 1.0;
        for (j = 0; j < i; j++)
          z(j, i) = z(i, j) = 0.0;

      } else {
        d[i] = z(i, i);
      }        
    }

    if (yesvecs) {
      result.Q_house = z;
    }

    result.d = d;
    result.e = e;

    return result;
}

QLEigenResult Matrix::QL(std::vector<double> d, std::vector<double> e) const {
  QLEigenResult result;
  int n = d.size();
  int m, l, iter, i, k;
  double s, r, p, g, f, dd, c, b;
  const double eps=std::numeric_limits<double>::epsilon();

  Matrix z = Matrix::Identity(n); // to store eigenvectors

  for (i = 1; i < n; i++) {
    e[i-1] = e[i];
  }
  e[n-1] = 0.0;

  for ( l=0; l<n; l++ ) {
    iter = 0;
    do {
      for ( m=l; m<n-1; m++ ) {
        dd = std::abs(d[m]) + std::abs(d[m+1]);
        if ( std::abs(e[m]) <= eps*dd ) break;
      }
      if ( m!=l ) {
        if ( iter++ == 30 ) throw std::runtime_error("Too many iterations in tqli");
        g = (d[l+1] - d[l]) / (2.0 * e[l]);
        r = pythag(g, 1.0);
        g = d[m] - d[l] + e[l]/(g + SIGN(r,g));
        s = 1.0;
        c = 1.0;
        p = 0.0;

        for (i = m-1; i >= l; i--) {
          f = s * e[i];
          b = c * e[i];

          e[i+1] = (r=pythag(f, g));
          if ( r == 0 ) {
            d[i+1] -= p;
            e[m] = 0.0;
            break;
          }
          s = f / r;
          c = g / r;
          g = d[i+1]-p;
          r =  (d[i] - g)*s+2.0*c*b;
          d[i+1]  = g + (p=s*r);
          g = c*r - b;
          
          // FORM EIGENVECTORS
          for ( k=0; k<n; k++ ) {
            f = z(k, i+1);
            z(k, i+1) = s*z(k, i)+c*f;
            z(k, i) = c*z(k, i)-s*f;
          }

        }
        if ( r== 0.0 && i >= l) continue;
        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
      }
    } while ( m != l );
  }

  result.eigenvalues = d;
  result.Q_ql = z;
  return result;
}

EigsymResult Matrix::eigsym() const {
  if (num_rows != num_cols) {
    throw InvalidMatrixSize("householder_tridiagonalize requires a square matrix");
  }

  if (!is_symmetric(1e-8)) {
    throw InvalidMatrixSize("Matrix must be symmetric for eigsym()");
  }
  
  EigsymResult result;
  
  TridiagonalResult tri = householder_tridiagonalize(true);

  QLEigenResult ql = QL(tri.d, tri.e);

  Matrix P = tri.Q_house * ql.Q_ql;

  result.eigenvalues = ql.eigenvalues;
  result.eigenvectors = P;

  return result;
}