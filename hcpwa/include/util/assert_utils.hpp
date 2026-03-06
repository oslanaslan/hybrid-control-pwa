#ifndef HCPWA_UTIL_ASSERT_UTILS_HPP
#define HCPWA_UTIL_ASSERT_UTILS_HPP

#include <Eigen/Core>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace hcpwa {
namespace util {

inline void assertShape(const Eigen::MatrixXd& matr, int expected_rows, int expected_cols) {
    if (matr.rows() != expected_rows || matr.cols() != expected_cols) {
        std::ostringstream oss;
        oss << "Matrix shape must be (" << expected_rows << ", " << expected_cols
            << "). Got: (" << matr.rows() << ", " << matr.cols() << ")";
        throw std::runtime_error(oss.str());
    }
}

inline void assertShape(const Eigen::VectorXd& vec, int expected_size) {
    if (vec.size() != expected_size) {
        std::ostringstream oss;
        oss << "Vector size must be " << expected_size << ". Got: " << vec.size();
        throw std::runtime_error(oss.str());
    }
}

inline void assertScalar(double value) {
    if (std::isnan(value) || std::isinf(value)) {
        std::ostringstream oss;
        oss << "Value must be a finite scalar. Got: " << value;
        throw std::runtime_error(oss.str());
    }
}

}  // namespace util
}  // namespace hcpwa

#endif  // HCPWA_UTIL_ASSERT_UTILS_HPP
