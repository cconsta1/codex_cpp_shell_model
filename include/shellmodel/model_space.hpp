#pragma once

#include <string>
#include <vector>

namespace shellmodel {

struct Orbital {
  std::string label;
  int n = 0;
  int l = 0;
  int two_j = 0;
  int two_m = 0;
  int isospin_z = 0;
  double energy = 0.0;
};

class ModelSpace {
 public:
  void add_orbital(Orbital orbital);
  [[nodiscard]] const std::vector<Orbital>& orbitals() const { return orbitals_; }
  [[nodiscard]] std::size_t size() const { return orbitals_.size(); }

 private:
  std::vector<Orbital> orbitals_;
};

}  // namespace shellmodel
