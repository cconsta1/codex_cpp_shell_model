#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace shellmodel::linalg {

class Matrix {
 public:
  Matrix() = default;
  Matrix(std::size_t rows, std::size_t cols, double value = 0.0)
      : rows_(rows), cols_(cols), data_(rows * cols, value) {}

  static Matrix zero(std::size_t rows, std::size_t cols) { return Matrix(rows, cols, 0.0); }

  [[nodiscard]] std::size_t rows() const { return rows_; }
  [[nodiscard]] std::size_t cols() const { return cols_; }

  double& operator()(std::size_t r, std::size_t c) { return data_[r * cols_ + c]; }
  double operator()(std::size_t r, std::size_t c) const { return data_[r * cols_ + c]; }

 private:
  std::size_t rows_ = 0;
  std::size_t cols_ = 0;
  std::vector<double> data_;
};

using Vector = std::vector<double>;

inline Matrix identity(std::size_t n) {
  Matrix id(n, n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    id(i, i) = 1.0;
  }
  return id;
}

inline double dot(const Vector& a, const Vector& b) {
  if (a.size() != b.size()) {
    throw std::invalid_argument("dot size mismatch");
  }
  double sum = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    sum += a[i] * b[i];
  }
  return sum;
}

inline Vector mat_vec(const Matrix& m, const Vector& v) {
  if (m.cols() != v.size()) {
    throw std::invalid_argument("mat_vec size mismatch");
  }
  Vector out(m.rows(), 0.0);
  for (std::size_t r = 0; r < m.rows(); ++r) {
    for (std::size_t c = 0; c < m.cols(); ++c) {
      out[r] += m(r, c) * v[c];
    }
  }
  return out;
}

inline Vector column(const Matrix& m, std::size_t c) {
  Vector out(m.rows(), 0.0);
  for (std::size_t r = 0; r < m.rows(); ++r) {
    out[r] = m(r, c);
  }
  return out;
}

inline void sort_eigensystem(Vector& values, Matrix& vectors) {
  std::vector<std::pair<double, std::size_t>> order;
  order.reserve(values.size());
  for (std::size_t i = 0; i < values.size(); ++i) {
    order.emplace_back(values[i], i);
  }
  std::sort(order.begin(), order.end(), [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; });

  Vector sorted_values(values.size(), 0.0);
  Matrix sorted_vectors(vectors.rows(), vectors.cols(), 0.0);
  for (std::size_t i = 0; i < order.size(); ++i) {
    sorted_values[i] = order[i].first;
    for (std::size_t r = 0; r < vectors.rows(); ++r) {
      sorted_vectors(r, i) = vectors(r, order[i].second);
    }
  }
  values = std::move(sorted_values);
  vectors = std::move(sorted_vectors);
}

}  // namespace shellmodel::linalg
