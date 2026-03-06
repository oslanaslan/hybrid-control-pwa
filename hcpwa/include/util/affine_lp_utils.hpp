#ifndef HCPWA_UTIL_AFFINE_LP_UTILS_HPP
#define HCPWA_UTIL_AFFINE_LP_UTILS_HPP

#include <Eigen/Core>

namespace hcpwa {
namespace util {

/// c^{(j)}(n) = Q_c^{(j)} n + q_c^{(j)}
inline Eigen::VectorXd centerC(const Eigen::MatrixXd& Qc,
                               const Eigen::VectorXd& qc,
                               const Eigen::VectorXd& n) {
    return Qc * n + qc;
}

/// r^{(j)}(n) = Q_r^{(j)} n + q_r^{(j)}
inline Eigen::VectorXd radiusR(const Eigen::MatrixXd& Qr,
                               const Eigen::VectorXd& qr,
                               const Eigen::VectorXd& n) {
    return Qr * n + qr;
}

/// p^{(j)}(n) = (A^{(j)} * dt - I) n + f^{(j)} * dt + dt * c^{(j)}(n)
inline Eigen::VectorXd coeffP(const Eigen::MatrixXd& Aj,
                              const Eigen::VectorXd& fj,
                              const Eigen::MatrixXd& Qc,
                              const Eigen::VectorXd& qc,
                              const Eigen::VectorXd& n,
                              double dt) {
    const int m = static_cast<int>(n.size());
    const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(m, m);
    const Eigen::VectorXd c = centerC(Qc, qc, n);
    return (Aj * dt - I) * n + fj * dt + dt * c;
}

/// g_i^{(j)}(n) = g_n^{(j)T} n + g_0^{(j)}
inline double gAffine(const Eigen::VectorXd& g_n, double g0, const Eigen::VectorXd& n) {
    return g_n.dot(n) + g0;
}

/// kappa_fix(n) = dt * g_i^{(j)}(n)
inline double kappaFixed(const Eigen::VectorXd& g_n,
                         double g0,
                         const Eigen::VectorXd& n,
                         double dt) {
    return dt * gAffine(g_n, g0, n);
}

}  // namespace util
}  // namespace hcpwa

#endif  // HCPWA_UTIL_AFFINE_LP_UTILS_HPP
