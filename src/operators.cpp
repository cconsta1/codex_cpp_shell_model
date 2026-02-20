#include "shellmodel/operators.hpp"

namespace shellmodel {

std::size_t TwoBodyKeyHash::operator()(const TwoBodyKey& key) const {
  std::size_t seed = 0;
  seed ^= static_cast<std::size_t>(key.a) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  seed ^= static_cast<std::size_t>(key.b) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  seed ^= static_cast<std::size_t>(key.c) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  seed ^= static_cast<std::size_t>(key.d) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  return seed;
}

namespace {

std::uint64_t pair_key(int a, int b) {
  return (static_cast<std::uint64_t>(a) << 32U) | static_cast<std::uint32_t>(b);
}

}  // namespace

void OneBodyOperator::set(int a, int b, double value) { values_[pair_key(a, b)] = value; }

double OneBodyOperator::get(int a, int b) const {
  const auto it = values_.find(pair_key(a, b));
  return it == values_.end() ? 0.0 : it->second;
}

void TwoBodyOperator::set(int a, int b, int c, int d, double value) {
  values_[TwoBodyKey{a, b, c, d}] = value;
}

double TwoBodyOperator::get(int a, int b, int c, int d) const {
  const auto it = values_.find(TwoBodyKey{a, b, c, d});
  return it == values_.end() ? 0.0 : it->second;
}

}  // namespace shellmodel
