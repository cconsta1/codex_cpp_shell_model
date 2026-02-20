#pragma once

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace shellmodel {

class SlaterBasis {
 public:
  SlaterBasis() = default;
  SlaterBasis(int n_particles, int n_single_particle_states);

  [[nodiscard]] int n_particles() const { return n_particles_; }
  [[nodiscard]] int n_states() const { return n_states_; }
  [[nodiscard]] const std::vector<std::uint64_t>& determinants() const { return determinants_; }
  [[nodiscard]] std::size_t dimension() const { return determinants_.size(); }

  [[nodiscard]] int index_of(std::uint64_t det) const;

 private:
  int n_particles_ = 0;
  int n_states_ = 0;
  std::vector<std::uint64_t> determinants_;
  std::unordered_map<std::uint64_t, int> index_map_;

  void generate();
};

}  // namespace shellmodel
