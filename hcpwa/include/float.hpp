#pragma once

#include <istream>
#include <ostream>

namespace hcpwa {

template <typename T, T epsilon = static_cast<T>(1e-5)>
class CustomFloat {
  static_assert(std::is_floating_point_v<T>,
                "CustomFloat only works with floating-point types");

 public:
  // NOLINTNEXTLINE
  constexpr CustomFloat(T value = static_cast<T>(0))
      : value_(value) {
  }

  // Conversion to float type
  constexpr explicit operator float() const {
    return value_;
  }
  // Conversion to double type
  constexpr explicit operator double() const {
    return value_;
  }

  // Compound assignment operators (must be members)
  constexpr CustomFloat& operator+=(const CustomFloat& rhs) {
    value_ += rhs.value_;
    return *this;
  }
  constexpr CustomFloat& operator-=(const CustomFloat& rhs) {
    value_ -= rhs.value_;
    return *this;
  }
  constexpr CustomFloat& operator*=(const CustomFloat& rhs) {
    value_ *= rhs.value_;
    return *this;
  }
  constexpr CustomFloat& operator/=(const CustomFloat& rhs) {
    value_ /= rhs.value_;
    return *this;
  }

  // Increment/decrement (must be members)
  constexpr CustomFloat& operator++() {
    ++value_;
    return *this;
  }
  constexpr CustomFloat operator++(int) {
    CustomFloat tmp(*this);
    ++value_;
    return tmp;
  }
  constexpr CustomFloat& operator--() {
    --value_;
    return *this;
  }
  constexpr CustomFloat operator--(int) {
    CustomFloat tmp(*this);
    --value_;
    return tmp;
  }
  constexpr CustomFloat operator-() {
    return -value_;
  }

  // Access to underlying value
  constexpr T Get() const {
    return value_;
  }

  // Static epsilon access
  static constexpr T GetEpsilon() {
    return epsilon;
  }

  // Non-member arithmetic operators

  constexpr friend CustomFloat<T, epsilon> operator+(
      CustomFloat<T, epsilon> lhs, const CustomFloat<T, epsilon>& rhs) {
    return lhs += rhs;
  }

  constexpr friend CustomFloat<T, epsilon> operator-(
      CustomFloat<T, epsilon> lhs, const CustomFloat<T, epsilon>& rhs) {
    return lhs -= rhs;
  }

  constexpr friend CustomFloat<T, epsilon> operator*(
      CustomFloat<T, epsilon> lhs, const CustomFloat<T, epsilon>& rhs) {
    return lhs *= rhs;
  }

  constexpr friend CustomFloat<T, epsilon> operator/(
      CustomFloat<T, epsilon> lhs, const CustomFloat<T, epsilon>& rhs) {
    return lhs /= rhs;
  }

  // Non-member comparison operators

  constexpr friend bool operator==(const CustomFloat<T, epsilon>& lhs,
                                   const CustomFloat<T, epsilon>& rhs) {
    return std::abs(lhs.Get() - rhs.Get()) <= epsilon;
  }

  constexpr friend bool operator!=(const CustomFloat<T, epsilon>& lhs,
                                   const CustomFloat<T, epsilon>& rhs) {
    return !(lhs == rhs);
  }

  constexpr friend bool operator<(const CustomFloat<T, epsilon>& lhs,
                                  const CustomFloat<T, epsilon>& rhs) {
    return lhs.Get() < rhs.Get();
  }

  constexpr friend bool operator<=(const CustomFloat<T, epsilon>& lhs,
                                   const CustomFloat<T, epsilon>& rhs) {
    return lhs.Get() <= rhs.Get();
  }

  constexpr friend bool operator>(const CustomFloat<T, epsilon>& lhs,
                                  const CustomFloat<T, epsilon>& rhs) {
    return lhs.Get() > rhs.Get();
  }

  constexpr friend bool operator>=(const CustomFloat<T, epsilon>& lhs,
                                   const CustomFloat<T, epsilon>& rhs) {
    return lhs.Get() >= rhs.Get();
  }

  // Non-member stream operators

  constexpr friend std::ostream& operator<<(std::ostream& os,
                                            const CustomFloat<T, epsilon>& cf) {
    return os << cf.Get();
  }

  constexpr friend std::istream& operator>>(std::istream& is,
                                            CustomFloat<T, epsilon>& cf) {
    T temp;
    is >> temp;
    cf = CustomFloat<T, epsilon>(temp);
    return is;
  }

  constexpr friend CustomFloat abs(const CustomFloat& value) {  // NOLINT
    return std::abs(value.value_);
  }

  T value_;
};

using Float = float;
}  // namespace hcpwa