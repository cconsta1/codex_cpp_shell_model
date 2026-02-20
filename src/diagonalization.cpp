#include "shellmodel/diagonalization.hpp"

#include <cmath>
#include <stdexcept>

namespace shellmodel {

EigenSystem diagonalize_hermitian(const linalg::Matrix& matrix, int max_iterations, double tolerance) {
  if (matrix.rows() != matrix.cols()) {
    throw std::invalid_argument("Matrix must be square");
  }
  const std::size_t n = matrix.rows();
  linalg::Matrix a = matrix;
  linalg::Matrix v = linalg::identity(n);

  for (int iter = 0; iter < max_iterations; ++iter) {
    std::size_t p = 0;
    std::size_t q = 1;
    double max_offdiag = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = i + 1; j < n; ++j) {
        const double val = std::abs(a(i, j));
        if (val > max_offdiag) {
          max_offdiag = val;
          p = i;
          q = j;
        }
      }
    }

    if (max_offdiag < tolerance) {
      break;
    }

    const double app = a(p, p);
    const double aqq = a(q, q);
    const double apq = a(p, q);
    const double tau = (aqq - app) / (2.0 * apq);
    const double t = (tau >= 0.0 ? 1.0 : -1.0) / (std::abs(tau) + std::sqrt(1.0 + tau * tau));
    const double c = 1.0 / std::sqrt(1.0 + t * t);
    const double s = t * c;

    for (std::size_t k = 0; k < n; ++k) {
      if (k == p || k == q) {
        continue;
      }
      const double aik = a(k, p);
      const double akq = a(k, q);
      a(k, p) = c * aik - s * akq;
      a(p, k) = a(k, p);
      a(k, q) = s * aik + c * akq;
      a(q, k) = a(k, q);
    }

    a(p, p) = c * c * app - 2.0 * s * c * apq + s * s * aqq;
    a(q, q) = s * s * app + 2.0 * s * c * apq + c * c * aqq;
    a(p, q) = 0.0;
    a(q, p) = 0.0;

    for (std::size_t k = 0; k < n; ++k) {
      const double vkp = v(k, p);
      const double vkq = v(k, q);
      v(k, p) = c * vkp - s * vkq;
      v(k, q) = s * vkp + c * vkq;
    }
  }

  linalg::Vector values(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    values[i] = a(i, i);
  }
  linalg::sort_eigensystem(values, v);
  return EigenSystem{values, v};
}

}  // namespace shellmodel
