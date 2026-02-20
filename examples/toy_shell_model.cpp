#include <iomanip>
#include <iostream>

#include "shellmodel/basis.hpp"
#include "shellmodel/diagonalization.hpp"
#include "shellmodel/hamiltonian.hpp"
#include "shellmodel/model_space.hpp"
#include "shellmodel/observables.hpp"
#include "shellmodel/operators.hpp"

int main() {
  using namespace shellmodel;

  ModelSpace space;
  space.add_orbital({"p3/2,m=-3/2", 0, 1, 3, -3, +1, 0.0});
  space.add_orbital({"p3/2,m=-1/2", 0, 1, 3, -1, +1, 0.0});
  space.add_orbital({"p3/2,m=+1/2", 0, 1, 3, +1, +1, 1.2});
  space.add_orbital({"p3/2,m=+3/2", 0, 1, 3, +3, +1, 1.2});

  SlaterBasis basis(2, static_cast<int>(space.size()));

  TwoBodyOperator interaction;
  interaction.set(0, 1, 0, 1, -1.0);
  interaction.set(1, 0, 1, 0, -1.0);
  interaction.set(2, 3, 2, 3, -0.7);
  interaction.set(3, 2, 3, 2, -0.7);
  interaction.set(0, 3, 1, 2, -0.3);
  interaction.set(1, 2, 0, 3, -0.3);

  const auto hamiltonian = HamiltonianBuilder::build(space, basis, interaction);
  const EigenSystem eig = diagonalize_hermitian(hamiltonian);

  OneBodyOperator e2_op;
  OneBodyOperator m1_op;
  for (int i = 0; i < static_cast<int>(space.size()); ++i) {
    e2_op.set(i, i, static_cast<double>(space.orbitals()[i].two_m * space.orbitals()[i].two_m) / 4.0);
    m1_op.set(i, i, 0.5 * static_cast<double>(space.orbitals()[i].two_m));
  }
  e2_op.set(1, 2, 0.4);
  e2_op.set(2, 1, 0.4);
  m1_op.set(0, 3, 0.2);
  m1_op.set(3, 0, 0.2);

  const auto e2_matrix = build_one_body_matrix(basis, e2_op);
  const auto m1_matrix = build_one_body_matrix(basis, m1_op);

  const auto ground_state = linalg::column(eig.eigenvectors, 0);
  const auto first_excited = linalg::column(eig.eigenvectors, 1);

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Toy shell-model example (2 particles in 4 states)\n";
  std::cout << "Basis dimension: " << basis.dimension() << "\n\n";
  std::cout << "Lowest energies [arb. units]:\n";
  for (std::size_t i = 0; i < eig.eigenvalues.size(); ++i) {
    std::cout << "  E(" << i << ") = " << eig.eigenvalues[i] << "\n";
  }

  std::cout << "\nToy electromagnetic observables:\n";
  std::cout << "  B(E2; 1->0) = " << transition_strength(first_excited, ground_state, e2_matrix) << "\n";
  std::cout << "  B(M1; 1->0) = " << transition_strength(first_excited, ground_state, m1_matrix) << "\n";
  std::cout << "  Q(ground)   = " << expectation_value(ground_state, e2_matrix) << "\n";
  std::cout << "  mu(ground)  = " << expectation_value(ground_state, m1_matrix) << "\n";
  return 0;
}
