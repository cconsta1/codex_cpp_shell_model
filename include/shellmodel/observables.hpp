#pragma once

#include "shellmodel/basis.hpp"
#include "shellmodel/linalg.hpp"
#include "shellmodel/operators.hpp"

namespace shellmodel {

linalg::Matrix build_one_body_matrix(const SlaterBasis& basis,
                                     const OneBodyOperator& operator_ob);

double transition_strength(const linalg::Vector& initial_state,
                           const linalg::Vector& final_state,
                           const linalg::Matrix& operator_matrix);

double expectation_value(const linalg::Vector& state,
                         const linalg::Matrix& operator_matrix);

}  // namespace shellmodel
