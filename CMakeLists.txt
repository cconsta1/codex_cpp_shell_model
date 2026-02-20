cmake_minimum_required(VERSION 3.18)
project(codex_cpp_shell_model VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_library(shellmodel
  src/basis.cpp
  src/diagonalization.cpp
  src/hamiltonian.cpp
  src/model_space.cpp
  src/observables.cpp
  src/operators.cpp
)

target_include_directories(shellmodel PUBLIC include)
target_compile_options(shellmodel PRIVATE -Wall -Wextra -Wpedantic)

add_executable(toy_shell_model examples/toy_shell_model.cpp)
target_link_libraries(toy_shell_model PRIVATE shellmodel)

add_executable(shellmodel_tests tests/test_shellmodel.cpp)
target_link_libraries(shellmodel_tests PRIVATE shellmodel)

enable_testing()
add_test(NAME shellmodel_tests COMMAND shellmodel_tests)
