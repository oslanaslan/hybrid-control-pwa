#ifndef HCPWA_UTIL_SPARSE_MATRIX_UTILS_HPP
#define HCPWA_UTIL_SPARSE_MATRIX_UTILS_HPP

#include <Eigen/Core>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace hcpwa {
namespace util {

/// Appends one CSR row from a dense row vector.
/// Skips entries with |val| <= eps when eps > 0.
///
/// The row pointer is a vector<int> because that is what HiGHS's addRows takes
/// (HighsInt is 32-bit in this build), so the nonzero count is bounded by
/// INT_MAX no matter how the accumulation is written. It is accumulated in
/// size_t and checked rather than silently wrapped: past the limit `row_ptr`
/// stops being monotone and addRows indexes outside col_idx and values, which
/// is a segfault inside HiGHS or, worse, a quietly corrupted LP.
///
/// This is closer than it looks. One global affine feasibility row carries up
/// to 17 nonzeros, so the ceiling is about 126 million rows, and fixing the
/// polygon vertex truncation multiplied that path's row count by 8.8.
inline void csrAppendRow(std::vector<int>& row_ptr,
                         std::vector<int>& col_idx,
                         std::vector<double>& values,
                         const Eigen::RowVectorXd& row,
                         double eps = 0.0) {
    std::size_t nnz = values.size();
    for (int i = 0; i < row.size(); ++i) {
        double v = row(i);
        if (eps > 0.0 && std::abs(v) <= eps) {
            continue;
        }
        col_idx.push_back(i);
        values.push_back(v);
        ++nnz;
    }
    if (nnz > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error(
            "csrAppendRow: the CSR matrix has more than INT_MAX nonzeros, "
            "which the 32-bit row pointer HiGHS takes cannot address");
    }
    row_ptr.push_back(static_cast<int>(nnz));
}

}  // namespace util
}  // namespace hcpwa

#endif  // HCPWA_UTIL_SPARSE_MATRIX_UTILS_HPP
