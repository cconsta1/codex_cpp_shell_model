#include "shellmodel/hamiltonian.hpp"

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

double one_body_element(std::uint64_t bra, std::uint64_t ket, const ModelSpace& model_space) {
  double value = 0.0;
  const int n = static_cast<int>(model_space.size());
  for (int p = 0; p < n; ++p) {
    const auto ann = annihilate(ket, p);
    if (!ann.valid) {
      continue;
    }
    const auto crt = create(ann.det, p);
    if (!crt.valid || crt.det != bra) {
      continue;
    }
    value += model_space.orbitals()[p].energy * static_cast<double>(ann.phase * crt.phase);
  }
  return value;
}

double two_body_element(std::uint64_t bra,
                        std::uint64_t ket,
                        const TwoBodyOperator& interaction,
                        int n_states) {
  double value = 0.0;
  for (int a = 0; a < n_states; ++a) {
    for (int b = 0; b < n_states; ++b) {
      for (int c = 0; c < n_states; ++c) {
        for (int d = 0; d < n_states; ++d) {
          const double tbme = interaction.get(a, b, c, d);
          if (tbme == 0.0) {
            continue;
          }
          const auto ann_d = annihilate(ket, d);
          if (!ann_d.valid) {
            continue;
          }
          const auto ann_c = annihilate(ann_d.det, c);
          if (!ann_c.valid) {
            continue;
          }
          const auto crt_b = create(ann_c.det, b);
          if (!crt_b.valid) {
            continue;
          }
          const auto crt_a = create(crt_b.det, a);
          if (!crt_a.valid || crt_a.det != bra) {
            continue;
          }
          value += 0.25 * tbme * static_cast<double>(ann_d.phase * ann_c.phase * crt_b.phase * crt_a.phase);
        }
      }
    }
  }
  return value;
}

}  // namespace

linalg::Matrix HamiltonianBuilder::build(const ModelSpace& model_space,
                                         const SlaterBasis& basis,
                                         const TwoBodyOperator& interaction) {
  const std::size_t dim = basis.dimension();
  linalg::Matrix hamiltonian = linalg::Matrix::zero(dim, dim);
  for (std::size_t i = 0; i < dim; ++i) {
    for (std::size_t j = i; j < dim; ++j) {
      const auto bra = basis.determinants()[i];
      const auto ket = basis.determinants()[j];
      const double element = one_body_element(bra, ket, model_space) +
                             two_body_element(bra, ket, interaction, basis.n_states());
      hamiltonian(i, j) = element;
      hamiltonian(j, i) = element;
    }
  }
  return hamiltonian;
}

}  // namespace shellmodel
