#include "shellmodel/basis.hpp"

#include <stdexcept>

namespace shellmodel {

namespace {

void choose_states(int start,
                   int remaining,
                   int n_states,
                   std::uint64_t det,
                   std::vector<std::uint64_t>& out) {
  if (remaining == 0) {
    out.push_back(det);
    return;
  }
  for (int i = start; i <= n_states - remaining; ++i) {
    choose_states(i + 1, remaining - 1, n_states, det | (1ULL << i), out);
  }
}

}  // namespace

SlaterBasis::SlaterBasis(int n_particles, int n_single_particle_states)
    : n_particles_(n_particles), n_states_(n_single_particle_states) {
  if (n_single_particle_states > 63) {
    throw std::invalid_argument("n_single_particle_states must be <= 63 for uint64_t bit encoding");
  }
  if (n_particles < 0 || n_particles > n_single_particle_states) {
    throw std::invalid_argument("Invalid particle count for basis generation");
  }
  generate();
}

void SlaterBasis::generate() {
  determinants_.clear();
  choose_states(0, n_particles_, n_states_, 0ULL, determinants_);
  index_map_.clear();
  for (std::size_t i = 0; i < determinants_.size(); ++i) {
    index_map_[determinants_[i]] = static_cast<int>(i);
  }
}

int SlaterBasis::index_of(std::uint64_t det) const {
  const auto it = index_map_.find(det);
  if (it == index_map_.end()) {
    return -1;
  }
  return it->second;
}

}  // namespace shellmodel
