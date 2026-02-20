#pragma once

#include <cstdint>
#include <unordered_map>

namespace shellmodel {

struct TwoBodyKey {
  int a;
  int b;
  int c;
  int d;

  bool operator==(const TwoBodyKey& other) const {
    return a == other.a && b == other.b && c == other.c && d == other.d;
  }
};

struct TwoBodyKeyHash {
  std::size_t operator()(const TwoBodyKey& key) const;
};

class OneBodyOperator {
 public:
  void set(int a, int b, double value);
  [[nodiscard]] double get(int a, int b) const;

 private:
  std::unordered_map<std::uint64_t, double> values_;
};

class TwoBodyOperator {
 public:
  void set(int a, int b, int c, int d, double value);
  [[nodiscard]] double get(int a, int b, int c, int d) const;

 private:
  std::unordered_map<TwoBodyKey, double, TwoBodyKeyHash> values_;
};

}  // namespace shellmodel
