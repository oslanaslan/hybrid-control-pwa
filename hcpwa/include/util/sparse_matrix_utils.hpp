#ifndef HCPWA_UTIL_SPARSE_MATRIX_UTILS_HPP
#define HCPWA_UTIL_SPARSE_MATRIX_UTILS_HPP

#include <Eigen/Core>
#include <cmath>
#include <vector>

namespace hcpwa {
namespace util {

/// Appends one CSR row from a dense row vector.
/// Skips entries with |val| <= eps when eps > 0.
inline void csrAppendRow(std::vector<int>& row_ptr,
                         std::vector<int>& col_idx,
                         std::vector<double>& values,
                         const Eigen::RowVectorXd& row,
                         double eps = 0.0) {
    int nnz = static_cast<int>(values.size());
    for (int i = 0; i < row.size(); ++i) {
        double v = row(i);
        if (eps > 0.0 && std::abs(v) <= eps) {
            continue;
        }
        col_idx.push_back(i);
        values.push_back(v);
        ++nnz;
    }
    row_ptr.push_back(nnz);
}

}  // namespace util
}  // namespace hcpwa

#endif  // HCPWA_UTIL_SPARSE_MATRIX_UTILS_HPP
