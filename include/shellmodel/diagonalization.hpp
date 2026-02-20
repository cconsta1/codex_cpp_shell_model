#pragma once

#include "shellmodel/linalg.hpp"

namespace shellmodel {

struct EigenSystem {
  linalg::Vector eigenvalues;
  linalg::Matrix eigenvectors;
};

EigenSystem diagonalize_hermitian(const linalg::Matrix& matrix,
                                  int max_iterations = 100,
                                  double tolerance = 1e-12);

}  // namespace shellmodel
