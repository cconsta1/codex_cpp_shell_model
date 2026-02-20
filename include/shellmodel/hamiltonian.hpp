#pragma once

#include "shellmodel/basis.hpp"
#include "shellmodel/linalg.hpp"
#include "shellmodel/model_space.hpp"
#include "shellmodel/operators.hpp"

namespace shellmodel {

class HamiltonianBuilder {
 public:
  static linalg::Matrix build(const ModelSpace& model_space,
                              const SlaterBasis& basis,
                              const TwoBodyOperator& interaction);
};

}  // namespace shellmodel
