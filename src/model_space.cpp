#include "shellmodel/model_space.hpp"

namespace shellmodel {

void ModelSpace::add_orbital(Orbital orbital) { orbitals_.push_back(std::move(orbital)); }

}  // namespace shellmodel
