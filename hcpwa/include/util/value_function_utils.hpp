#ifndef HCPWA_UTIL_VALUE_FUNCTION_UTILS_HPP
#define HCPWA_UTIL_VALUE_FUNCTION_UTILS_HPP

#include <cstddef>
#include <fstream>
#include <sstream>
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

inline std::vector<double> getValueFunction(const ValueFunction& value_function,
                                            int phase, int r, int theta_idx,
                                            int theta_end_idx) {
    try {
        return value_function.at(phase).at(r).at(theta_idx).at(theta_end_idx);
    } catch (const std::out_of_range&) {
        throw std::invalid_argument(
            "getValueFunction: location not found for phase="
            + std::to_string(phase) + " r=" + std::to_string(r)
            + " theta_idx=" + std::to_string(theta_idx)
            + " theta_end_idx=" + std::to_string(theta_end_idx));
    }
}

inline void dumpVectorToJson(const std::vector<double>& vec,
                            const std::string& filepath) {
    std::ostringstream out;
    out << '[';
    for (std::size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) {
            out << ',';
        }
        out << vec[i];
    }
    out << ']';

    std::ofstream f(filepath);
    if (!f) {
        throw std::runtime_error("dumpVectorToJson: cannot open file " +
                                filepath + " for writing");
    }
    f << out.str();
}

namespace detail {

inline void dumpValueFunctionArray(std::ostringstream& out,
                                   const std::vector<double>& vec) {
    out << '[';
    for (std::size_t i = 0; i < vec.size(); ++i) {
        if (i > 0) {
            out << ',';
        }
        out << vec[i];
    }
    out << ']';
}

inline void dumpValueFunctionLevel4(
    std::ostringstream& out,
    const std::unordered_map<int, std::vector<double>>& level, bool first) {
    for (const auto& [theta_end_idx, vec] : level) {
        if (!first) {
            out << ',';
        }
        first = false;
        out << '"' << theta_end_idx << "\":";
        dumpValueFunctionArray(out, vec);
    }
}

inline void dumpValueFunctionLevel3(
    std::ostringstream& out,
    const std::unordered_map<int, std::unordered_map<int, std::vector<double>>>&
        level,
    bool first) {
    for (const auto& [theta_idx, inner] : level) {
        if (!first) {
            out << ',';
        }
        first = false;
        out << '"' << theta_idx << "\":{";
        dumpValueFunctionLevel4(out, inner, true);
        out << '}';
    }
}

inline void dumpValueFunctionLevel2(
    std::ostringstream& out,
    const std::unordered_map<
        int,
        std::unordered_map<int, std::unordered_map<int, std::vector<double>>>>&
        level,
    bool first) {
    for (const auto& [r, inner] : level) {
        if (!first) {
            out << ',';
        }
        first = false;
        out << '"' << r << "\":{";
        dumpValueFunctionLevel3(out, inner, true);
        out << '}';
    }
}

}  // namespace detail

inline void dumpValueFunctionToJson(const ValueFunction& value_function,
                                    const std::string& filepath) {
    std::ostringstream out;
    out << '{';
    bool first_phase = true;
    for (const auto& [phase, level] : value_function) {
        if (!first_phase) {
            out << ',';
        }
        first_phase = false;
        out << '"' << phase << "\":{";
        detail::dumpValueFunctionLevel2(out, level, true);
        out << '}';
    }
    out << '}';

    std::ofstream f(filepath);
    if (!f) {
        throw std::runtime_error("dumpValueFunctionToJson: cannot open file " +
                                filepath + " for writing");
    }
    f << out.str();
}

}  // namespace util
}  // namespace hcpwa

#endif  // HCPWA_UTIL_VALUE_FUNCTION_UTILS_HPP
