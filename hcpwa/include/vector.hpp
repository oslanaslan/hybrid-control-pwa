#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <iostream>

namespace linalg {

template<typename T, int N>
class Vector {
    static_assert(N > 0, "Vector dimension must be greater than 0");

private:
    std::array<T, N> data_;

public:
    // Type aliases
    using ValueType = T;
    using SizeType = size_t;
    using Reference = T&;
    using ConstReference = const T&;
    using Iterator = typename std::array<T, N>::iterator;
    using ConstIterator = typename std::array<T, N>::const_iterator;

    // Constructors
    constexpr Vector() : data_{} {
        data_.fill(T{});
    }

    explicit constexpr Vector(const T& value) : data_{} {
        data_.fill(value);
    }

    template<std::convertible_to<T>... Args>
    explicit constexpr Vector(Args... args) : data_{static_cast<T>(args)...} {
        static_assert(sizeof...(Args) <= N, "Number of arguments must be less than or equal to vector dimension");
    }

    explicit constexpr Vector(const std::array<T, N>& arr) : data_(arr) {}

    Vector(std::initializer_list<T> init) {
        std::copy(init.begin(), init.end(), data_.begin());
    }

    template<int M> requires(M <= N)
    Vector(const Vector<T, M>& other) { // NOLINT
        std::copy(other.begin(), other.end(), data_.begin());
    }

    template<int M> requires(M < N)
    explicit Vector(const Vector<T, M>& other, T fill) {
        std::copy(other.begin(), other.end(), data_.begin());
        data_[M] = fill;
    }

    // Copy constructor
    constexpr Vector(const Vector& other) = default;

    // Move constructor
    constexpr Vector(Vector&& other) noexcept = default;

    // Assignment operators
    constexpr Vector& operator=(const Vector& other) = default;
    constexpr Vector& operator=(Vector&& other) noexcept = default;

    // Element access
    constexpr Reference operator[](SizeType index) {
        return data_[index];
    }

    constexpr ConstReference operator[](SizeType index) const {
        return data_[index];
    }

    constexpr Reference At(SizeType index) {
        if (index >= N) {
            throw std::out_of_range("Vector index out of range");
        }
        return data_[index];
    }

    constexpr ConstReference At(SizeType index) const {
        if (index >= N) {
            throw std::out_of_range("Vector index out of range");
        }
        return data_[index];
    }

    constexpr Reference Front() { return data_.front(); }
    constexpr ConstReference Front() const { return data_.front(); }
    constexpr Reference Back() { return data_.back(); }
    constexpr ConstReference Back() const { return data_.back(); }

    // Iterators
    constexpr Iterator begin() { return data_.begin(); } // NOLINT
    constexpr ConstIterator begin() const { return data_.begin(); } // NOLINT
    constexpr ConstIterator cbegin() const { return data_.cbegin(); } // NOLINT
    constexpr Iterator end() { return data_.end(); } // NOLINT
    constexpr ConstIterator end() const { return data_.end(); } // NOLINT
    constexpr ConstIterator cend() const { return data_.cend(); } // NOLINT

    // Capacity
    constexpr SizeType Size() const noexcept { return N; }
    constexpr SizeType MaxSize() const noexcept { return N; }
    constexpr bool Empty() const noexcept { return false; }

    // Data access
    constexpr T* Data() { return data_.data(); }
    constexpr const T* Data() const { return data_.data(); }

    // Mathematical operations

    // Unary operators
    constexpr Vector operator+() const {
        return *this;
    }

    constexpr Vector operator-() const {
        Vector result;
        for (SizeType i = 0; i < N; ++i) {
            result[i] = -data_[i];
        }
        return result;
    }

    // Arithmetic assignment operators
    constexpr Vector& operator+=(const Vector& other) {
        for (SizeType i = 0; i < N; ++i) {
            data_[i] += other.data_[i];
        }
        return *this;
    }

    constexpr Vector& operator-=(const Vector& other) {
        for (SizeType i = 0; i < N; ++i) {
            data_[i] -= other.data_[i];
        }
        return *this;
    }

    constexpr Vector& operator*=(const T& scalar) {
        for (SizeType i = 0; i < N; ++i) {
            data_[i] *= scalar;
        }
        return *this;
    }

    Vector& operator/=(const T& scalar) {
        if (scalar == T{}) {
            throw std::invalid_argument("Division by zero");
        }
        for (SizeType i = 0; i < N; ++i) {
            data_[i] /= scalar;
        }
        return *this;
    }

    // Binary arithmetic operators
    friend constexpr Vector operator+(const Vector& lhs, const Vector& rhs) {
        Vector result = lhs;
        result += rhs;
        return result;
    }

    friend constexpr Vector operator-(const Vector& lhs, const Vector& rhs) {
        Vector result = lhs;
        result -= rhs;
        return result;
    }

    friend constexpr Vector operator*(const Vector& vec, const T& scalar) {
        Vector result = vec;
        result *= scalar;
        return result;
    }

    friend constexpr Vector operator*(const T& scalar, const Vector& vec) {
        return vec * scalar;
    }

    friend constexpr Vector operator/(const Vector& vec, const T& scalar) {
        Vector result = vec;
        result /= scalar;
        return result;
    }

    // Comparison operators
    constexpr bool operator==(const Vector& other) const {
        for (SizeType i = 0; i < N; ++i) {
            if (data_[i] != other.data_[i]) {
                return false;
            }
        }
        return true;
    }

    constexpr bool operator!=(const Vector& other) const {
        return !(*this == other);
    }

    // Mathematical functions

    // Dot product
    constexpr T Dot(const Vector& other) const {
        T result = T{};
        for (SizeType i = 0; i < N; ++i) {
            result += data_[i] * other.data_[i];
        }
        return result;
    }

    // Cross product (only for 3D vectors)
    template<int M = N>
    constexpr std::enable_if_t<M == 3, Vector> Cross(const Vector& other) const {
        return Vector{
            data_[1] * other.data_[2] - data_[2] * other.data_[1],
            data_[2] * other.data_[0] - data_[0] * other.data_[2],
            data_[0] * other.data_[1] - data_[1] * other.data_[0]
        };
    }

    // Magnitude (Euclidean norm)
    T Magnitude() const {
        T sum = T{};
        for (SizeType i = 0; i < N; ++i) {
            sum += data_[i] * data_[i];
        }
        return std::sqrt(sum);
    }

    // Squared magnitude (for efficiency when comparing magnitudes)
    constexpr T MagnitudeSquared() const {
        T sum = T{};
        for (SizeType i = 0; i < N; ++i) {
            sum += data_[i] * data_[i];
        }
        return sum;
    }

    // Normalize vector (unit vector)
    Vector Normalized() const {
        T mag = Magnitude();
        if (mag == T{}) {
            throw std::runtime_error("Cannot normalize zero vector");
        }
        return *this / mag;
    }

    // Normalize in place
    Vector& Normalize() {
        T mag = Magnitude();
        if (mag == T{}) {
            throw std::runtime_error("Cannot normalize zero vector");
        }
        *this /= mag;
        return *this;
    }

    // Distance to another vector
    T Distance(const Vector& other) const {
        return (*this - other).Magnitude();
    }

    // Squared distance to another vector
    constexpr T DistanceSquared(const Vector& other) const {
        return (*this - other).MagnitudeSquared();
    }

    // Angle between vectors (in radians)
    T Angle(const Vector& other) const {
        T mag1 = Magnitude();
        T mag2 = other.Magnitude();
        if (mag1 == T{} || mag2 == T{}) {
            throw std::runtime_error("Cannot compute angle with zero vector");
        }
        T cos_angle = Dot(other) / (mag1 * mag2);
        // Clamp to avoid numerical errors
        cos_angle = std::max(T{-1}, std::min(T{1}, cos_angle));
        return std::acos(cos_angle);
    }

    // Projection onto another vector
    Vector Projection(const Vector& other) const {
        T other_mag_sq = other.MagnitudeSquared();
        if (other_mag_sq == T{}) {
            throw std::runtime_error("Cannot project onto zero vector");
        }
        return other * (Dot(other) / other_mag_sq);
    }

    // Reflection across another vector
    Vector Reflection(const Vector& other) const {
        return *this - T{2} * Projection(other);
    }

    // Linear interpolation
    constexpr Vector Lerp(const Vector& other, const T& t) const {
        return *this + t * (other - *this);
    }

    // Spherical linear interpolation (for unit vectors)
    Vector Slerp(const Vector& other, const T& t) const {
        T dot_product = Dot(other);
        // Clamp to avoid numerical errors
        dot_product = std::max(T{-1}, std::min(T{1}, dot_product));
        
        T theta = std::acos(dot_product) * t;
        Vector relative_vec = other - *this * dot_product;
        relative_vec.Normalize();
        
        return *this * std::cos(theta) + relative_vec * std::sin(theta);
    }

    // Utility functions

    // Check if vector is zero
    constexpr bool IsZero() const {
        for (SizeType i = 0; i < N; ++i) {
            if (data_[i] != T{}) {
                return false;
            }
        }
        return true;
    }

    // Check if vector is normalized (within tolerance)
    bool IsNormalized(const T& tolerance = T{1e-6}) const {
        return std::abs(Magnitude() - T{1}) < tolerance;
    }

    // Get minimum element
    constexpr T Min() const {
        T min_val = data_[0];
        for (SizeType i = 1; i < N; ++i) {
            if (data_[i] < min_val) {
                min_val = data_[i];
            }
        }
        return min_val;
    }

    // Get maximum element
    constexpr T Max() const {
        T max_val = data_[0];
        for (SizeType i = 1; i < N; ++i) {
            if (data_[i] > max_val) {
                max_val = data_[i];
            }
        }
        return max_val;
    }

    // Get sum of all elements
    constexpr T Sum() const {
        T result = T{};
        for (SizeType i = 0; i < N; ++i) {
            result += data_[i];
        }
        return result;
    }

    // Get mean of all elements
    constexpr T Mean() const {
        return Sum() / static_cast<T>(N);
    }

    // Stream output
    friend std::ostream& operator<<(std::ostream& os, const Vector& vec) {
        os << "(";
        for (SizeType i = 0; i < N; ++i) {
            os << vec[i];
            if (i < N - 1) {
                os << ", ";
            }
        }
        os << ")";
        return os;
    }

    template<std::convertible_to<size_t>... Args> requires(sizeof...(Args) > 1)
    constexpr auto operator[](Args... args) const {
        return Vector<T, sizeof...(Args)>(data_[args]...);
    }

    std::vector<T> ToVector() const {
        return std::vector<T>(data_.begin(), data_.end());
    }
};

template<typename T>
Vector<T, 3> cross(const Vector<T, 3>& a, const Vector<T, 3>& b) {  // NOLINT
    return Vector<T, 3>{
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    };
}

template<typename T, int N>
T dot(const Vector<T, N>& a, const Vector<T, N>& b) {  // NOLINT
    T result = T{};
    for (size_t i = 0; i < N; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

} // namespace linalg
