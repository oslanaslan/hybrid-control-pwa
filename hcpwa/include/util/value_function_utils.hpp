#ifndef HCPWA_UTIL_VALUE_FUNCTION_UTILS_HPP
#define HCPWA_UTIL_VALUE_FUNCTION_UTILS_HPP

#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace hcpwa {
namespace util {

// ValueFunction[phase][r][theta_idx][theta_end_idx] = [V, v] or affine params
using ValueFunction = std::unordered_map<
    int,
    std::unordered_map<
        int,
        std::unordered_map<int, std::unordered_map<int, std::vector<double>>>>>;

inline void setValueFunction(ValueFunction& value_function, int phase, int r,
                             int theta_idx, int theta_end_idx,
                             const std::vector<double>& x_prev,
                             std::size_t expected_size) {
    if (x_prev.size() != expected_size) {
        throw std::invalid_argument("x_prev size must be "
                                    + std::to_string(expected_size));
    }
    auto& target = value_function[phase][r][theta_idx][theta_end_idx];
    if (!target.empty()) {
        throw std::runtime_error(
            "setValueFunction: location already used for phase="
            + std::to_string(phase) + " r=" + std::to_string(r)
            + " theta_idx=" + std::to_string(theta_idx)
            + " theta_end_idx=" + std::to_string(theta_end_idx));
    }
    target = x_prev;
}

inline std::vector<double> getValueFunction(ValueFunction& value_function,
                                            int phase, int r, int theta_idx,
                                            int theta_end_idx) {
    try {
        return value_function[phase][r][theta_idx][theta_end_idx];
    } catch (const std::out_of_range&) {
        throw std::invalid_argument(
            "getValueFunction: location not found for phase="
            + std::to_string(phase) + " r=" + std::to_string(r)
            + " theta_idx=" + std::to_string(theta_idx)
            + " theta_end_idx=" + std::to_string(theta_end_idx));
    }
}

}  // namespace util
}  // namespace hcpwa

#endif  // HCPWA_UTIL_VALUE_FUNCTION_UTILS_HPP
