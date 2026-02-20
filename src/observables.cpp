#include "shellmodel/observables.hpp"

#include <cstdint>

namespace shellmodel {
namespace {

struct OpResult {
  bool valid = false;
  std::uint64_t det = 0;
  int phase = 1;
};

int parity_below(std::uint64_t det, int idx) {
  if (idx <= 0) {
    return 1;
  }
  const std::uint64_t mask = (1ULL << idx) - 1ULL;
  const int count = static_cast<int>(__builtin_popcountll(det & mask));
  return (count % 2 == 0) ? 1 : -1;
}

OpResult annihilate(std::uint64_t det, int idx) {
  if (((det >> idx) & 1ULL) == 0ULL) {
    return {};
  }
  return OpResult{true, det & ~(1ULL << idx), parity_below(det, idx)};
}

OpResult create(std::uint64_t det, int idx) {
  if (((det >> idx) & 1ULL) == 1ULL) {
    return {};
  }
  return OpResult{true, det | (1ULL << idx), parity_below(det, idx)};
}

}  // namespace

linalg::Matrix build_one_body_matrix(const SlaterBasis& basis, const OneBodyOperator& operator_ob) {
  const std::size_t dim = basis.dimension();
  linalg::Matrix matrix = linalg::Matrix::zero(dim, dim);

  for (std::size_t i = 0; i < dim; ++i) {
    const auto bra = basis.determinants()[i];
    for (std::size_t j = 0; j < dim; ++j) {
      const auto ket = basis.determinants()[j];
      double element = 0.0;
      for (int a = 0; a < basis.n_states(); ++a) {
        for (int b = 0; b < basis.n_states(); ++b) {
          const double ob = operator_ob.get(a, b);
          if (ob == 0.0) {
            continue;
          }
          const auto ann = annihilate(ket, b);
          if (!ann.valid) {
            continue;
          }
          const auto crt = create(ann.det, a);
          if (!crt.valid || crt.det != bra) {
            continue;
          }
          element += ob * static_cast<double>(ann.phase * crt.phase);
        }
      }
      matrix(i, j) = element;
    }
  }
  return matrix;
}

double transition_strength(const linalg::Vector& initial_state,
                           const linalg::Vector& final_state,
                           const linalg::Matrix& operator_matrix) {
  const linalg::Vector op_initial = linalg::mat_vec(operator_matrix, initial_state);
  const double amplitude = linalg::dot(final_state, op_initial);
  return amplitude * amplitude;
}

double expectation_value(const linalg::Vector& state, const linalg::Matrix& operator_matrix) {
  const linalg::Vector op_state = linalg::mat_vec(operator_matrix, state);
  return linalg::dot(state, op_state);
}

}  // namespace shellmodel
