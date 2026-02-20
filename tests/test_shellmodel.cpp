#include <cmath>
#include <iostream>
#include <stdexcept>

#include "shellmodel/basis.hpp"
#include "shellmodel/diagonalization.hpp"
#include "shellmodel/hamiltonian.hpp"
#include "shellmodel/linalg.hpp"
#include "shellmodel/model_space.hpp"
#include "shellmodel/observables.hpp"
#include "shellmodel/operators.hpp"

namespace {

void expect_true(bool condition, const char* message) {
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void expect_near(double lhs, double rhs, double eps, const char* message) {
  if (std::abs(lhs - rhs) > eps) {
    throw std::runtime_error(message);
  }
}

void test_basis_dimension() {
  shellmodel::SlaterBasis basis(2, 4);
  expect_true(basis.dimension() == 6, "Basis dimension should be C(4,2)=6");
}

void test_non_interacting_ground_energy() {
  using namespace shellmodel;
  ModelSpace space;
  space.add_orbital({"a", 0, 0, 1, -1, +1, 0.0});
  space.add_orbital({"b", 0, 0, 1, +1, +1, 1.0});
  space.add_orbital({"c", 0, 0, 1, -1, +1, 2.0});

  SlaterBasis basis(2, 3);
  TwoBodyOperator interaction;
  const auto h = HamiltonianBuilder::build(space, basis, interaction);
  const auto eig = diagonalize_hermitian(h);
  expect_near(eig.eigenvalues[0], 1.0, 1e-10, "Ground-state energy should be sum of two lowest SP energies");
}

void test_transition_strength() {
  shellmodel::linalg::Matrix op(2, 2, 0.0);
  op(0, 1) = 1.0;
  op(1, 0) = 1.0;
  shellmodel::linalg::Vector psi_i{1.0, 0.0};
  shellmodel::linalg::Vector psi_f{0.0, 1.0};
  const double strength = shellmodel::transition_strength(psi_i, psi_f, op);
  expect_near(strength, 1.0, 1e-12, "Transition strength for basis flip should be 1");
}

}  // namespace

int main() {
  try {
    test_basis_dimension();
    test_non_interacting_ground_energy();
    test_transition_strength();
    std::cout << "All tests passed.\n";
    return 0;
  } catch (const std::exception& ex) {
    std::cerr << "Test failed: " << ex.what() << '\n';
    return 1;
  }
}
