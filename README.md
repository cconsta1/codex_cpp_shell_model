# Minimal C++20 Shell-Model MVP

This repository contains a **minimal but real** shell-model codebase in modern C++20 for personal computers (macOS + Linux).

It is intentionally scoped as an MVP:
- fixed particle number in an **m-scheme Slater determinant basis**,
- one major shell / small valence space,
- toy two-body matrix elements (TBMEs),
- dense diagonalization,
- basic electromagnetic-style observables (`B(E2)`, `B(M1)`, moments) through one-body operators.

## Dependencies
- CMake >= 3.18
- C++20 compiler (Clang or GCC)
- No external linear-algebra dependency: this MVP includes a tiny header-only dense matrix/vector layer (`include/shellmodel/linalg.hpp`) with a Jacobi solver for symmetric matrices.

No heavyweight frameworks are required.

## Project layout

```text
.
├── CMakeLists.txt
├── README.md
├── include/shellmodel/
│   ├── basis.hpp
│   ├── diagonalization.hpp
│   ├── hamiltonian.hpp
│   ├── model_space.hpp
│   ├── observables.hpp
│   └── operators.hpp
├── src/
│   ├── basis.cpp
│   ├── diagonalization.cpp
│   ├── hamiltonian.cpp
│   ├── model_space.cpp
│   ├── observables.cpp
│   └── operators.cpp
├── examples/
│   └── toy_shell_model.cpp
└── tests/
    └── test_shellmodel.cpp
```

## Build and run

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/toy_shell_model
```

## Example output (toy case)

Representative output for the included toy interaction:

```text
Toy shell-model example (2 particles in 4 states)
Basis dimension: 6

Lowest energies [arb. units]:
  E(0) = 0.500000
  E(1) = 1.125000
  E(2) = 1.200000
  E(3) = 1.200000
  E(4) = 1.275000
  E(5) = 2.750000

Toy electromagnetic observables:
  B(E2; 1->0) = 0.000000
  B(M1; 1->0) = 0.000000
  Q(ground)   = 2.500000
  mu(ground)  = -2.000000
```

Values are in arbitrary units for the toy interaction/operator set.

## Theory assumptions and limitations

1. **m-scheme basis only**
   - States are Slater determinants encoded as 64-bit occupancy bit patterns.
   - Current implementation supports up to 63 single-particle states.

2. **Hamiltonian form**
   - One-body part uses orbital single-particle energies.
   - Two-body part uses provided TBMEs in second-quantized form:
     \[
       H_2 = \frac{1}{4}\sum_{abcd} V_{ab,cd} a^\dagger_a a^\dagger_b a_d a_c
     \]

3. **Observables**
   - One-body operators are provided in m-scheme matrix elements.
   - `B(E2)`/`B(M1)` are currently computed as simple squared amplitudes
     \(|\langle f|\hat O|i\rangle|^2\) without full reduced-matrix formalism,
     angular-momentum coupling, or effective charges/g-factors.

4. **Numerics**
   - Dense matrices with an internal Jacobi solver; appropriate only for small spaces.

## Extensibility roadmap

- Add TBME reader for realistic interaction files (e.g., `J,T`-coupled tables converted to m-scheme).
- Add proton-neutron spaces and explicit isospin handling.
- Add angular-momentum projection / coupled-J basis.
- Replace dense diagonalization with Lanczos for larger dimensions.
- Add proper reduced transition rates and electromagnetic operators with effective charges/g-factors.
