#include "global_affine_approximator.h"

#include <algo.hpp>
#include <algorithm>
#include <cstddef>
#include <format>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>
#include "cddwrap/lineareq.hpp"
#include "morph.hpp"
#include "types.hpp"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Highs.h>
#include <spdlog/spdlog.h>
#include "spdlog/sinks/stdout_color_sinks.h"

namespace global_affine_approximator {

const std::vector<int> kInIds = {2 - 1, 3 - 1, 5 - 1, 8 - 1};
const std::vector<int> kOutIds = {1 - 1, 4 - 1, 6 - 1, 7 - 1};
const int kPhases = 2;

void setValueFunction(ValueFunction& value_function, int phase, int r, int theta_idx, int theta_end_idx, const std::vector<double>& x_prev) {
    // x_prev is a flattened vector if size 2 * (SPACE_DIM + 1), consisting of [V^+, V^-, v^+, v^-]
    if (x_prev.size() != 2 * (kSpaceDim + 1)) {
        throw std::invalid_argument("x_prev size must be " + std::to_string(kSpaceDim));
    }
    // Check if this value_function location is already used (not empty)
    auto& target = value_function[phase][r][theta_idx][theta_end_idx];
    if (target.size() > 0) {
        throw std::runtime_error("setValueFunction: location already used for phase=" + std::to_string(phase) +
                                 " r=" + std::to_string(r) +
                                 " theta_idx=" + std::to_string(theta_idx) +
                                 " theta_end_idx=" + std::to_string(theta_end_idx));
    }
    target = x_prev;
}

std::vector<double> getValueFunction(ValueFunction& value_function, int phase, int r, int theta_idx, int theta_end_idx) {
    try {
        return value_function[phase][r][theta_idx][theta_end_idx];
    } catch (const std::out_of_range& e) {
        throw std::invalid_argument("getValueFunction: location not found for phase=" + std::to_string(phase) +
                                 " r=" + std::to_string(r) +
                                 " theta_idx=" + std::to_string(theta_idx) +
                                 " theta_end_idx=" + std::to_string(theta_end_idx));
    }
}

GlobalAffineApproximator::GlobalAffineApproximator(double t_max,
                                       int t_split_count,
                                       int max_switches,
                                       double tau_min,
                                       double tau_max,
                                       const SystemParams& system_params)
    : t_max_(t_max),
      t_split_count_(t_split_count),
      max_switches_(max_switches),
      tau_min_(tau_min),
      tau_max_(tau_max),
      system_params_(system_params) {
    if (t_split_count < 1) {
        throw std::invalid_argument(
            "t_split_count must be at least 1 in GlobalAffineApproximator constructor.");
    }
    if (tau_min >= tau_max || tau_min < 0.0 || tau_max < 0.0) {
        throw std::invalid_argument(
            "tau_min must be less than tau_max and both must be >= 0 in GlobalAffineApproximator constructor.");
    }
    if (max_switches < 1) {
        throw std::invalid_argument(
            "max_switches must be at least 1 in GlobalAffineApproximator constructor.");
    }
    t_range_.resize(t_split_count);
    t_delta_ = t_split_count > 1 ? std::abs(t_max / (t_split_count - 1)) : 0.0;
    for (int i = 0; i < t_split_count; ++i) {
        t_range_[i] = (t_split_count == 1) ? t_max : i * t_delta_;
    }
    t_index_.resize(t_split_count);
    for (int i = 0; i < t_split_count; ++i) {
        t_index_[i] = i;
    }
    int n_vertices = 1 << kSpaceDim;
    cube_angle_vertices_.clear();
    for (int vert = 0; vert < n_vertices; ++vert) {
        Eigen::VectorXd v(kSpaceDim);
        for (int d = 0; d < kSpaceDim; ++d) {
            v(d) = ((vert & (1 << d)) != 0) ? this->system_params_.N : 0.0;
        }
        cube_angle_vertices_.push_back(v);
    }
    logger_ = spdlog::stdout_color_mt("affine_approximator");
    logger_->set_level(spdlog::level::info);
}

std::tuple<double, double, double, double> GlobalAffineApproximator::getMainSystemParams() const {
        return std::make_tuple(
            system_params_.N,
            system_params_.F,
            system_params_.v,
            system_params_.w
        );
    }

double GlobalAffineApproximator::getBetaParamForAxis(int i, int j) const {
        if (i == 5 - 1 && j == 1 - 1) {
            return system_params_.b51;
        } else if (i == 5 - 1 && j == 7 - 1) {
            return system_params_.b57;
        } else if (i == 8 - 1 && j == 4 - 1) {
            return system_params_.b84;
        } else if (i == 8 - 1 && j == 6 - 1) {
            return system_params_.b86;
        } else if (i == 3 - 1 && j == 1 - 1) {
            return system_params_.b31;
        } else if (i == 3 - 1 && j == 6 - 1) {
            return system_params_.b36;
        } else if (i == 2 - 1 && j == 4 - 1) {
            return system_params_.b24;
        } else if (i == 2 - 1 && j == 7 - 1) {
            return system_params_.b27;
        } else {
            throw std::invalid_argument(
                std::string("No such beta params for axis (") +
                std::to_string(i) + ", " + std::to_string(j) + ")"
            );
        }
    }

std::pair<double, double> GlobalAffineApproximator::getFMinMaxForAxis(int i) const {
        if (i == 2 - 1) {
            return std::make_pair(system_params_.f2min, system_params_.f2max);
        } else if (i == 3 - 1) {
            return std::make_pair(system_params_.f3min, system_params_.f3max);
        } else if (i == 5 - 1) {
            return std::make_pair(system_params_.f5min, system_params_.f5max);
        } else if (i == 8 - 1) {
            return std::make_pair(system_params_.f8min, system_params_.f8max);
        } else {
            throw std::invalid_argument(
                std::string("No such f min max for axis ") + std::to_string(i)
            );
        }
    }

void GlobalAffineApproximator::getIntersectionPoints() {
    hcpwa::AreasVerticesResult areas_vertices = hcpwa::compute_areas_vertices(
        system_params_.N,
        system_params_.F,
        system_params_.v,
        system_params_.w,
        system_params_.b51,
        system_params_.b57,
        system_params_.b84,
        system_params_.b86,
        system_params_.b31,
        system_params_.b36,
        system_params_.b24,
        system_params_.b27,
        system_params_.f2min,
        system_params_.f3min,
        system_params_.f5min,
        system_params_.f8min,
        system_params_.f2max,
        system_params_.f3max,
        system_params_.f5max,
        system_params_.f8max);

    intersection_points_.resize(2);

    auto convert_phase = [](const std::vector<std::vector<hcpwa::Vec<8>>>& src,
                            std::vector<std::vector<Eigen::VectorXd>>& dst) {
        dst.resize(src.size());
        for (size_t a = 0; a < src.size(); ++a) {
            const auto& area = src[a];
            dst[a].resize(area.size());
            for (size_t v = 0; v < area.size(); ++v) {
                Eigen::VectorXd vec(kSpaceDim);
                for (int i = 0; i < kSpaceDim; ++i) {
                    vec(i) = static_cast<double>(area[v][i]);
                }
                dst[a][v] = std::move(vec);
            }
        }
    };

    convert_phase(areas_vertices.intersection_points_phase0, intersection_points_[0]);
    convert_phase(areas_vertices.intersection_points_phase1, intersection_points_[1]);
}

std::pair<Eigen::RowVectorXd, Eigen::RowVectorXd> GlobalAffineApproximator::getFIJMinResolution(int i, int j, const Eigen::VectorXd& n) const {
        // f_i_j_min_a = min{beta_i_j * F, beta_i_j * v * n_i, w(N − n_j)}
        // Where f_matr_row is row vector [length SPACE_DIM], f_vec_row is 1x1 (just value)
        const double N = system_params_.N;
        const double F = system_params_.F;
        const double v = system_params_.v;
        const double w = system_params_.w;
        double beta_i_j = getBetaParamForAxis(i, j);
        double n_i = n(i);
        double n_j = n(j);
        Eigen::RowVectorXd f_matr_row = Eigen::RowVectorXd::Zero(kSpaceDim); // 1 x SPACE_DIM
        Eigen::RowVectorXd f_vec_row = Eigen::RowVectorXd::Zero(1); // 1 x 1

        double a = beta_i_j * v * n_i;
        double b = beta_i_j * F;
        double c = w * (N - n_j);

        if (b < a + kEps && b < c + kEps) {
            // Case: b is smallest
            f_vec_row(0) = beta_i_j * F;
        } else if (c < a + kEps && c < b + kEps) {
            // Case: c is smallest
            f_matr_row(j) = -w;
            f_vec_row(0) = w * N;
        } else if (a < b + kEps && a < c + kEps) {
            // Case: a is smallest
            f_matr_row(i) = beta_i_j * v;
            // f_vec_row remains zero
        } else {
            std::ostringstream oss;
            oss << "No such case: a=" << a << ", b=" << b << ", c=" << c << " (getFIJMinResolution)";
            throw std::invalid_argument(oss.str());
        }
        return std::make_pair(f_matr_row, f_vec_row);
    }

Eigen::VectorXd GlobalAffineApproximator::areaCentroidCoords(int j, int phase) const {
        // Check if intersection_points are computed for given phase
        if (intersection_points_.empty() || phase < 0 || phase >= static_cast<int>(intersection_points_.size()) || intersection_points_[phase].empty()) {
            throw std::runtime_error("Intersection points are not computed for the specified phase. Please compute intersection points before calling areaCentroidCoords.");
        }
        const std::vector<std::vector<Eigen::VectorXd>>& intersection_points_per_phase = intersection_points_[phase];
        if (j < 0 || j >= static_cast<int>(intersection_points_per_phase.size())) {
            throw std::invalid_argument("Invalid area index j in areaCentroidCoords");
        }
        
        const auto& area = intersection_points_per_phase[j];
        if (area.empty()) {
            throw std::invalid_argument("Area is empty in areaCentroidCoords");
        }
        
        Eigen::VectorXd centroid = Eigen::VectorXd::Zero(kSpaceDim);
        for (const auto& point : area) {
            centroid += point;
        }
        centroid /= static_cast<double>(area.size());
        
        assertShape(centroid, kSpaceDim);
        return centroid;
    }

void GlobalAffineApproximator::assertShape(const Eigen::MatrixXd& matr, int expected_rows, int expected_cols) {
        if (matr.rows() != expected_rows || matr.cols() != expected_cols) {
            std::ostringstream oss;
            oss << "Matrix shape must be (" << expected_rows << ", " << expected_cols 
                << "). Got: (" << matr.rows() << ", " << matr.cols() << ")";
            throw std::runtime_error(oss.str());
        }
    }

void GlobalAffineApproximator::assertShape(const Eigen::VectorXd& vec, int expected_size) {
        if (vec.size() != expected_size) {
            std::ostringstream oss;
            oss << "Vector size must be " << expected_size << ". Got: " << vec.size();
            throw std::runtime_error(oss.str());
        }
    }

void GlobalAffineApproximator::assertScalar(double value) {
        // In C++, double is always scalar, but we can check for NaN/Inf
        if (std::isnan(value) || std::isinf(value)) {
            std::ostringstream oss;
            oss << "Value must be a finite scalar. Got: " << value;
            throw std::runtime_error(oss.str());
        }
    }

std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, double>
GlobalAffineApproximator::getAMatrFVecGVecAndGScalJ(int j, int phase) const {
        Eigen::VectorXd n = areaCentroidCoords(j, phase);
        
        Eigen::MatrixXd a_matr;
        Eigen::VectorXd b_vec;
        Eigen::VectorXd g_vec;
        double g_scal;
        
        if (phase == 0) {
            auto [f31_matr_row, f31_vec_row] = getFIJMinResolution(3 - 1, 1 - 1, n);
            auto [f36_matr_row, f36_vec_row] = getFIJMinResolution(3 - 1, 6 - 1, n);
            auto [f24_matr_row, f24_vec_row] = getFIJMinResolution(2 - 1, 4 - 1, n);
            auto [f27_matr_row, f27_vec_row] = getFIJMinResolution(2 - 1, 7 - 1, n);
            
            // Build A_matr: SPACE_DIM x SPACE_DIM
            a_matr.resize(kSpaceDim, kSpaceDim);
            a_matr.row(0) = f31_matr_row;
            a_matr.row(1) = -f24_matr_row - f27_matr_row;
            a_matr.row(2) = -f31_matr_row - f36_matr_row;
            a_matr.row(3) = f24_matr_row;
            a_matr.row(4) = Eigen::RowVectorXd::Zero(kSpaceDim);
            a_matr.row(5) = f36_matr_row;
            a_matr.row(6) = f27_matr_row;
            a_matr.row(7) = Eigen::RowVectorXd::Zero(kSpaceDim);
            
            // Build b_vec: SPACE_DIM x 1
            b_vec.resize(kSpaceDim);
            b_vec(0) = f31_vec_row(0);
            b_vec(1) = -f24_vec_row(0) - f27_vec_row(0);
            b_vec(2) = -f31_vec_row(0) - f36_vec_row(0);
            b_vec(3) = f24_vec_row(0);
            b_vec(4) = 0.0;
            b_vec(5) = f36_vec_row(0);
            b_vec(6) = f27_vec_row(0);
            b_vec(7) = 0.0;
            
            // Build g_vec: SPACE_DIM x 1
            g_vec = (f31_matr_row + f36_matr_row + f24_matr_row + f27_matr_row).transpose();
            
            // Build g_scal: scalar
            g_scal = f31_vec_row(0) + f36_vec_row(0) + f24_vec_row(0) + f27_vec_row(0);
        } else if (phase == 1) {
            auto [f51_matr_row, f51_vec_row] = getFIJMinResolution(5 - 1, 1 - 1, n);
            auto [f57_matr_row, f57_vec_row] = getFIJMinResolution(5 - 1, 7 - 1, n);
            auto [f84_matr_row, f84_vec_row] = getFIJMinResolution(8 - 1, 4 - 1, n);
            auto [f86_matr_row, f86_vec_row] = getFIJMinResolution(8 - 1, 6 - 1, n);
            
            // Build A_matr: SPACE_DIM x SPACE_DIM
            a_matr.resize(kSpaceDim, kSpaceDim);
            a_matr.row(0) = f51_matr_row;
            a_matr.row(1) = Eigen::RowVectorXd::Zero(kSpaceDim);
            a_matr.row(2) = Eigen::RowVectorXd::Zero(kSpaceDim);
            a_matr.row(3) = f84_matr_row;
            a_matr.row(4) = -f51_matr_row - f57_matr_row;
            a_matr.row(5) = f86_matr_row;
            a_matr.row(6) = f57_matr_row;
            a_matr.row(7) = -f84_matr_row - f86_matr_row;
            
            // Build b_vec: SPACE_DIM x 1
            b_vec.resize(kSpaceDim);
            b_vec(0) = f51_vec_row(0);
            b_vec(1) = 0.0;
            b_vec(2) = 0.0;
            b_vec(3) = f84_vec_row(0);
            b_vec(4) = -f51_vec_row(0) - f57_vec_row(0);
            b_vec(5) = f86_vec_row(0);
            b_vec(6) = f57_vec_row(0);
            b_vec(7) = -f84_vec_row(0) - f86_vec_row(0);
            
            // Build g_vec: SPACE_DIM x 1
            g_vec = (f51_matr_row + f57_matr_row + f84_matr_row + f86_matr_row).transpose();
            
            // Build g_scal: scalar
            g_scal = f51_vec_row(0) + f57_vec_row(0) + f84_vec_row(0) + f86_vec_row(0);
        } else {
            std::ostringstream oss;
            oss << "No such phase: " << phase;
            throw std::invalid_argument(oss.str());
        }
        
        // Assert shapes
        assertShape(b_vec, kSpaceDim);
        assertShape(a_matr, kSpaceDim, kSpaceDim);
        assertShape(g_vec, kSpaceDim);
        assertScalar(g_scal);
        
        return std::make_tuple(a_matr, b_vec, g_vec, g_scal);
    }

std::pair<Eigen::RowVectorXd, double> GlobalAffineApproximator::getMaxEstimForOutIds(
    int i, double n_i, bool v_sign) const {
        const double F = system_params_.F;
        const double v = system_params_.v;

        if (v_sign) {
            return std::make_pair(Eigen::RowVectorXd::Zero(kSpaceDim), 0.0);
        }

        Eigen::RowVectorXd q_mat_row = Eigen::RowVectorXd::Zero(kSpaceDim);
        double q_vec_row;

        if (v * n_i <= F) {
            q_mat_row(i) = -v;
            q_vec_row = 0.0;
        } else {
            q_mat_row(i) = 0.0;
            q_vec_row = -F;
        }

        return std::make_pair(q_mat_row, q_vec_row);
    }

std::pair<Eigen::RowVectorXd, double> GlobalAffineApproximator::getMaxEstimForInIds(
    int i, double n_i, bool v_sign) const {
        Eigen::RowVectorXd q_mat_row = Eigen::RowVectorXd::Zero(kSpaceDim);
        double q_vec_row = 0.0;
        
        const double N = system_params_.N;
        const double w = system_params_.w;
        auto [f_min, f_max] = getFMinMaxForAxis(i);

        if (v_sign) {
            if (w * (N - n_i) <= f_max) {
                q_mat_row(i) = -w;
                q_vec_row = w * N;
            } else {
                q_mat_row(i) = 0.0;
                q_vec_row = f_max;
            }
        } else {
            if (w * (N - n_i) <= f_min) {
                q_mat_row(i) = -w;
                q_vec_row = w * N;
            } else {
                q_mat_row(i) = 0.0;
                q_vec_row = f_min;
            }
        }

        return std::make_pair(q_mat_row, q_vec_row);
    }

std::pair<Eigen::MatrixXd, Eigen::VectorXd> GlobalAffineApproximator::getQQForArea(
    int j, int phase) const {
        // Get representative point in Ω(j)
        Eigen::VectorXd n0 = areaCentroidCoords(j, phase);

        Eigen::MatrixXd Q_matr = Eigen::MatrixXd::Zero(2 * kSpaceDim, kSpaceDim);
        Eigen::VectorXd q_vec = Eigen::VectorXd::Zero(2 * kSpaceDim);

        // Process IN_IDS
        for (int i : kInIds) {
            // max V_i * f_{i, in} = V_i^+ * u_i(n) - V_i^- * l_i(n)
            double n_i = n0(i);
            
            // V^+ coefs (Right border of the segment)
            auto [q_mat_row, q_vec_row] = getMaxEstimForInIds(i, n_i, true);
            Q_matr.row(i) = q_mat_row;
            q_vec(i) = q_vec_row;
            
            // V^- coefs (Left border of the segment)
            auto [q_mat_row_neg, q_vec_row_neg] = getMaxEstimForInIds(i, n_i, false);
            Q_matr.row(i + kSpaceDim) = -q_mat_row_neg;
            q_vec(i + kSpaceDim) = -q_vec_row_neg;
        }

        // Process OUT_IDS
        for (int i : kOutIds) {
            // max V_i * f_{i, out} = V_i^- * u_i(n)
            double n_i = n0(i);
            
            // V^- coefs
            auto [q_mat_row, q_vec_row] = getMaxEstimForOutIds(i, n_i, false);  // TODO Why no True option here? Is it u bug?
            Q_matr.row(i) = q_mat_row;
            q_vec(i) = q_vec_row;
        }

        return std::make_pair(Q_matr, q_vec);
    }

std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, double>
GlobalAffineApproximator::getBMatAndBVecJ(int j, int phase) const {
        // Get A, f, g, g_scal for area j and phase
        auto [A_j_matr, f_j_vec, g_j_vec, g_j_scal] = getAMatrFVecGVecAndGScalJ(j, phase);
        
        // V = [V^+, V^-]^T => V^T A = (V^+ - V^-)^T*A = [V^+, V^-]*[A, -A]^T
        // Concatenate A_j_matr with -A_j_matr vertically
        Eigen::MatrixXd A_j_matr_concat(2 * kSpaceDim, kSpaceDim);
        A_j_matr_concat.topRows(kSpaceDim) = A_j_matr;
        A_j_matr_concat.bottomRows(kSpaceDim) = -A_j_matr;
        assertShape(A_j_matr_concat, 2 * kSpaceDim, kSpaceDim);
        
        // Concatenate f_j_vec with -f_j_vec vertically
        Eigen::VectorXd f_j_vec_concat(2 * kSpaceDim);
        f_j_vec_concat.head(kSpaceDim) = f_j_vec;
        f_j_vec_concat.tail(kSpaceDim) = -f_j_vec;
        assertShape(f_j_vec_concat, 2 * kSpaceDim);
        
        // Get Q and q for area j and phase
        auto [Q_j_matr, q_j_vec] = getQQForArea(j, phase);
        assertShape(Q_j_matr, 2 * kSpaceDim, kSpaceDim);
        assertShape(q_j_vec, 2 * kSpaceDim);
        
        // B_part_matr = (A_j_matr + Q_j_matr) * t_delta
        Eigen::MatrixXd B_part_matr = (A_j_matr_concat + Q_j_matr) * t_delta_;
        
        // Subtract [I; -I] from B_part_matr
        // ... -V = ... -(V^+ - V^-) = ... [V^+, V^-] * [I, -I]^T
        Eigen::MatrixXd identity_neg_identity(2 * kSpaceDim, kSpaceDim);
        identity_neg_identity.topRows(kSpaceDim) = Eigen::MatrixXd::Identity(kSpaceDim, kSpaceDim);
        identity_neg_identity.bottomRows(kSpaceDim) = -Eigen::MatrixXd::Identity(kSpaceDim, kSpaceDim);
        B_part_matr -= identity_neg_identity;
        
        // Concatenate B_part_matr with zeros(2, SPACE_DIM) to get B_j_matr
        // Add 2 zeros for the last component v = v^+ - v^- that are not used in this equation
        Eigen::MatrixXd B_j_matr(2 * kSpaceDim + 2, kSpaceDim);
        B_j_matr.topRows(2 * kSpaceDim) = B_part_matr;
        B_j_matr.bottomRows(2) = Eigen::MatrixXd::Zero(2, kSpaceDim);
        assertShape(B_j_matr, 2 * kSpaceDim + 2, kSpaceDim);
        
        // V^T * (f + q) * \delta_t − v = (V^+ - V^-)^T * (f + q) * \delta_t − (v^+ - v^-)
        Eigen::VectorXd b_j_part_vec = (f_j_vec_concat + q_j_vec) * t_delta_;
        
        // Concatenate b_j_part_vec with [-1, 1]^T to get b_j_vec
        // Add last two dims for the "− (v^+ - v^-)" component
        Eigen::VectorXd b_j_vec(2 * kSpaceDim + 2);
        b_j_vec.head(2 * kSpaceDim) = b_j_part_vec;
        b_j_vec(2 * kSpaceDim) = -1.0;
        b_j_vec(2 * kSpaceDim + 1) = 1.0;
        assertShape(b_j_vec, 2 * kSpaceDim + 2);
        
        return std::make_tuple(B_j_matr, b_j_vec, g_j_vec, g_j_scal);
    }

std::tuple<std::vector<Eigen::MatrixXd>, std::vector<Eigen::VectorXd>, std::vector<Eigen::VectorXd>, std::vector<double>>
GlobalAffineApproximator::precomputeBMatAndBVec(int phase) {
        auto intersection_points_phase = intersection_points_[phase];
        int n_areas = intersection_points_phase.size();
        std::vector<Eigen::MatrixXd> B_j_matrs_lst;
        std::vector<Eigen::VectorXd> b_j_vecs_lst;
        std::vector<Eigen::VectorXd> g_j_vecs_lst;
        std::vector<double> g_j_scals_lst;
        B_j_matrs_lst.reserve(n_areas);
        b_j_vecs_lst.reserve(n_areas);
        g_j_vecs_lst.reserve(n_areas);
        g_j_scals_lst.reserve(n_areas);

        for (int j = 0; j < n_areas; ++j) {
            auto [B_j_matr, b_j_vec, g_j_vec, g_j_scal] = getBMatAndBVecJ(j, phase);
            B_j_matrs_lst.emplace_back(B_j_matr);
            b_j_vecs_lst.emplace_back(b_j_vec);
            g_j_vecs_lst.emplace_back(g_j_vec);
            g_j_scals_lst.emplace_back(g_j_scal);
        }

        return std::make_tuple(B_j_matrs_lst, b_j_vecs_lst, g_j_vecs_lst, g_j_scals_lst);
    }

std::pair<std::vector<Eigen::VectorXd>, std::vector<double>>
GlobalAffineApproximator::getCVecDScalLists(int phase,
                                      const Eigen::VectorXd& x_prev_vec,
                                      const std::vector<Eigen::VectorXd>& g_j_vecs,
                                      const std::vector<double>& g_j_scals) const {
        if (x_prev_vec.size() != 2 * kSpaceDim + 2) {
            throw std::invalid_argument(
                "getCVecDScalLists: x_prev_vec must have size 2*kSpaceDim+2.");
        }
        if (g_j_vecs.size() != g_j_scals.size()) {
            throw std::invalid_argument(
                "getCVecDScalLists: g_j_vecs and g_j_scals must have the same size.");
        }

        int n_areas = static_cast<int>(g_j_vecs.size());
        std::vector<Eigen::VectorXd> c_j_vecs_lst;
        c_j_vecs_lst.reserve(n_areas);
        std::vector<double> d_j_scals_lst;
        d_j_scals_lst.reserve(n_areas);

        for (int j = 0; j < n_areas; ++j) {
            const Eigen::VectorXd& g_j_vec = g_j_vecs[j];
            double g_j_scal = g_j_scals[j];
            // c_j_vec = V_{t+1}^T + GT * \delta_t
            Eigen::VectorXd c_j_vec =
                x_prev_vec.head(kSpaceDim)
                - x_prev_vec.segment(kSpaceDim, kSpaceDim)
                + g_j_vec * t_delta_;
            assertShape(c_j_vec, kSpaceDim);
            double d_j_scal =
                x_prev_vec(2 * kSpaceDim) - x_prev_vec(2 * kSpaceDim + 1)
                + g_j_scal * t_delta_;
            c_j_vecs_lst.emplace_back(c_j_vec);
            d_j_scals_lst.emplace_back(d_j_scal);
        }

        return std::make_pair(c_j_vecs_lst, d_j_scals_lst);
    }

/**

Prepare matrices for LP solver.

Solves the following LP:

 *   Objective: maximize z over z and x
 *
 *   Subject to:
 *     max_{n in Omega^(j)} [x^T * B^(j) * n + x^T * b^(j) + c^(j)^T * n + d^(j)] <= 0,  for all j = 1...M
 *     z <= min_{n in Omega^(j)} [x^T * B^(j) * n + x^T * b^(j) + c^(j)^T * n + d^(j)],  for all j = 1...M

@param self: pointer to GlobalAffineApproximator object
@param phase: phase index
@return: tuple of vectors:
    - starts: list of indices where each constraint starts
    - col_index: list of column indices for each non-zero element
    - value: list of values for each non-zero element
    - row_upper: list of upper bounds for each constraint
    - c_vec: list of coefficients for the objective function
*/
std::tuple<std::vector<int>, std::vector<int>, std::vector<double>, std::vector<double>, Eigen::RowVectorXd>
GlobalAffineApproximator::prepareLpMatrices(int phase) {
        const std::vector<Eigen::MatrixXd>& b_j_matrs = this->b_j_matrs_[phase];
        const std::vector<Eigen::VectorXd>& b_j_vecs = this->b_j_vecs_[phase];
        const std::vector<Eigen::VectorXd>& g_j_vecs = this->g_j_vecs_[phase];
        const std::vector<double>& g_j_scals = this->g_j_scals_[phase];
        const auto& intersection_points_phase = intersection_points_[phase];
        int n_areas = intersection_points_phase.size();
        std::vector<int> starts;
        std::vector<int> col_index;
        std::vector<double> value;
        std::vector<double> row_upper;
        Eigen::RowVectorXd c_vec = Eigen::RowVectorXd::Zero(kLPCols);
        c_vec(kLPCols - 2) = 1.0; // z^+
        c_vec(kLPCols - 1) = -1.0; // z^-
        size_t nnz = 0;

        for (size_t j = 0; j < n_areas; ++j) {
            const Eigen::MatrixXd& b_j_matr = b_j_matrs[j];
            const Eigen::VectorXd& b_j_vec = b_j_vecs[j];
            const Eigen::VectorXd& g_j_vec = g_j_vecs[j];
            const double& g_j_scal = g_j_scals[j];

            Eigen::VectorXd c_j_vec = g_j_vec * t_delta_;

            for (const auto& vertex : intersection_points_phase[j]) {
                Eigen::VectorXd vertex_vec = Eigen::Map<const Eigen::VectorXd>(vertex.data(), kSpaceDim);
                Eigen::RowVectorXd alpha_j = b_j_matr * vertex_vec + b_j_vec;
                // True beta_j = c_j_vec.dot(vertex_vec) + d_j_scal = c_j_vec^T * vertex_vec + v_prev^+ - v_prev^- + g_j_scal * \delta_t
                // But for faster RHS updates we use tilde_beta_j = c_j_vec.dot(vertex_vec) * t_delta_ + g_j_scal * t_delta_ = c_j_vec^T * vertex_vec * t_delta_ + g_j_scal * \delta_t
                // The v_prev^+ - v_prev must be added manually in the RHS updates.
                double tilde_beta_j = c_j_vec.dot(vertex_vec) + g_j_scal * t_delta_;
                assertShape(alpha_j, kLPCols - 2);

                // Add constraint: [alpha_j, 0, 0] * x <= -beta_j, where alpha_j is a row of A_matr stored in CSR format
                // This represents the constraint: max_{n in Omega^(j)} [x^T * B^(j) * n + x^T * b^(j) + c^(j)^T * n + d^(j)] <= 0,  for all j = 1...M
                starts.push_back(nnz);
                for (size_t i = 0; i < kLPCols - 2; ++i) {
                    if (std::abs(alpha_j(i)) > kEps) {
                        col_index.push_back(i);
                        value.push_back(alpha_j(i));
                        nnz++;
                    }
                }
                row_upper.push_back(-tilde_beta_j);

                // Add constraint: [-alpha_j, 1, -1] * x <= beta_j, where beta_j is a row of A_matr stored in CSR format
                // This represents the constraint: z <= min_{n in Omega^(j)} [x^T * B^(j) * n + x^T * b^(j) + c^(j)^T * n + d^(j)],  for all j = 1...M
                starts.push_back(nnz);
                for (size_t i = 0; i < kLPCols - 2; ++i) {
                    if (std::abs(alpha_j(i)) > kEps) {
                        col_index.push_back(i);
                        value.push_back(-alpha_j(i));
                        nnz++;
                    }
                }
                col_index.push_back(kLPCols - 2);
                value.push_back(1.0);
                col_index.push_back(kLPCols - 1);
                value.push_back(-1.0);
                nnz += 2;
                row_upper.push_back(tilde_beta_j);
            }
        }
        starts.push_back(nnz);

        return std::make_tuple(std::move(starts), std::move(col_index), std::move(value), std::move(row_upper), std::move(c_vec));
    }

double GlobalAffineApproximator::getBorderFuncValuesAtN(int r, int theta_idx, int theta_end_idx, int phase, const Eigen::VectorXd& n) {
        try {
            std::vector<double> x = getValueFunction(this->value_function_, phase, r, theta_idx, theta_end_idx);
            // x is [V^+, V^-, v^+, v^-] of size 2 * (kSpaceDim + 1)
            // We need to compute: [n, 1]^T @ [V^+ - V^-, v^+ - v^-]
            // where V^+ - V^- is the first kSpaceDim elements of (x[:kSpaceDim] - x[kSpaceDim:2*kSpaceDim])
            // and v^+ - v^- is (x[2*kSpaceDim] - x[2*kSpaceDim+1])
            
            double result = 0.0;
            for (int i = 0; i < kSpaceDim; ++i) {
                result += (x[i] - x[kSpaceDim + i]) * n(i);
            }
            result += x[2 * kSpaceDim] - x[2 * kSpaceDim + 1];
            return result;
        } catch (const std::invalid_argument&) {
            return -std::numeric_limits<double>::infinity();
        }
    }

double GlobalAffineApproximator::getMaxBorderFuncValuesAtN(int theta_idx,
                                                    const std::vector<int>& theta_end_ids,
                                                    int max_switches,
                                                    int phase,
                                                    const Eigen::VectorXd& n) {
        double max_val = -std::numeric_limits<double>::infinity();
        for (int theta_end_idx : theta_end_ids) {
            for (int r = 0; r < max_switches; ++r) {
                double val = getBorderFuncValuesAtN(r, theta_idx, theta_end_idx, phase, n);
                max_val = std::max(max_val, val);
            }
        }
        return max_val;
    }

std::vector<double> GlobalAffineApproximator::getBorderConditions(int switch_phase, int theta_idx, double theta, int switch_cnt) {
        // If switch_cnt == 0, return zeros
        if (switch_cnt == 0) {
            return std::vector<double>(2 * (kSpaceDim + 1), 0.0);
        }

        // Calculate theta_min and theta_max
        double theta_min = std::min(theta + this->tau_min_, this->t_max_);
        double theta_max = std::min(theta + this->tau_max_, this->t_max_);

        // Find theta_range_ids using loop (replacing NumPy boolean indexing)
        // avoid numerical errors by adding a half-size step to the theta_min and theta_max
        double half_step = this->t_delta_ / 2.0;
        theta_min -= half_step;
        theta_max += half_step;
        std::vector<int> theta_range_ids;
        for (size_t i = 0; i < this->t_range_.size(); ++i) {
            if (this->t_range_[i] >= theta_min && this->t_range_[i] <= theta_max) {
                theta_range_ids.push_back(this->t_index_[i]);
            }
        }

        if (theta_range_ids.empty()) {
            throw std::runtime_error("getBorderConditions: No theta range found");
        }

        // Initialize c_vec as zeros
        std::vector<double> c_vec(2 * (kSpaceDim + 1), 0.0);
        
        // Build a_matr and b_vec
        std::vector<std::vector<double>> a_matr_lst;
        std::vector<double> b_vec_lst;
        
        for (const auto& vertex : this->cube_angle_vertices_) {
            double f_scal = getMaxBorderFuncValuesAtN(
                theta_idx, theta_range_ids, switch_cnt, switch_phase, vertex
            );
            
            // Build row: [vertex, -vertex, 1, -1]
            std::vector<double> row(2 * (kSpaceDim + 1));
            for (int i = 0; i < kSpaceDim; ++i) {
                row[i] = vertex(i);                    // V^+
                row[kSpaceDim + i] = -vertex(i);       // V^-
            }
            row[2 * kSpaceDim] = 1.0;                  // v^+
            row[2 * kSpaceDim + 1] = -1.0;             // v^-
            
            // Accumulate into c_vec
            for (size_t i = 0; i < c_vec.size(); ++i) {
                c_vec[i] += row[i];
            }
            
            b_vec_lst.push_back(f_scal);
            a_matr_lst.push_back(row);
        }

        // Convert a_matr_lst to a single matrix (for Highs sparse format)
        const int m = static_cast<int>(a_matr_lst.size());  // 2^kSpaceDim
        const int n = static_cast<int>(c_vec.size());       // 2 * (kSpaceDim + 1)

        // Build sparse matrix representation for Highs
        // For constraints: -a_matr @ x <= -b_vec  (from Python: A_ub=-a_matr, b_ub=-b_vec)
        // In Highs format: we need A x <= b, so we use -a_matr and -b_vec
        std::vector<int> starts(m + 1, 0);
        std::vector<int> col_index;
        std::vector<double> value;
        
        int nnz = 0;
        for (int i = 0; i < m; ++i) {
            starts[i] = nnz;
            for (int j = 0; j < n; ++j) {
                const double a_val = -a_matr_lst[i][j];  // Negate for constraint format
                if (std::abs(a_val) > 1e-10) {  // Skip near-zero values
                    col_index.push_back(j);
                    value.push_back(a_val);
                    ++nnz;
                }
            }
        }
        starts[m] = nnz;

        // Prepare row bounds: -a_matr @ x <= -b_vec means row_upper = -b_vec, row_lower = -inf
        std::vector<double> row_lower(m, -std::numeric_limits<double>::infinity());
        std::vector<double> row_upper(m);
        for (int i = 0; i < m; ++i) {
            row_upper[i] = -b_vec_lst[i];  // Negate b_vec for constraint format
        }

        // Set up and solve LP with Highs
        Highs highs;
        highs.setOptionValue("solver", "simplex");
        highs.setOptionValue("presolve", "on");
        
        const double inf = highs.getInfinity();

        // Add columns (variables)
        std::vector<double> col_lower(n, 0.0);
        std::vector<double> col_upper(n, inf);
        
        HighsStatus st = highs.addCols(
            n,
            c_vec.data(),
            col_lower.data(),
            col_upper.data(),
            0,
            nullptr,
            nullptr,
            nullptr
        );
        if (st != HighsStatus::kOk) {
            throw std::runtime_error("getBorderConditions: highs.addCols failed");
        }

        // Add rows (constraints)
        st = highs.addRows(
            m,
            row_lower.data(),
            row_upper.data(),
            nnz,
            starts.data(),
            col_index.data(),
            value.data()
        );
        if (st != HighsStatus::kOk) {
            throw std::runtime_error("getBorderConditions: highs.addRows failed");
        }

        // Solve
        st = highs.run();
        if (st != HighsStatus::kOk) {
            throw std::runtime_error("getBorderConditions: highs.run() failed");
        }

        if (highs.getModelStatus() != HighsModelStatus::kOptimal) {
            throw std::runtime_error("getBorderConditions: LP solution not found");
        }

        // Extract solution
        const auto& solution = highs.getSolution();
        std::vector<double> result(solution.col_value.begin(), solution.col_value.end());
        
        return result;
    }

std::tuple<std::unique_ptr<Highs>, std::vector<double>, std::vector<double>>
GlobalAffineApproximator::initializeHighs(int phase) {
        auto [starts, col_index, value, row_upper, c_vec] = prepareLpMatrices(phase);
        const int m = row_upper.size();
        const int n = c_vec.size();

        std::unique_ptr<Highs> highs = std::make_unique<Highs>();

        // ---- Solver options (set once) ----
        highs->setOptionValue("solver", "simplex");
        // Depending on HiGHS version, some options may not exist; handle failures if needed.
        highs->setOptionValue("presolve", "on");
        highs->changeObjectiveSense(ObjSense::kMaximize); // Must be maximization because need to maximize the lower bound of the non-positive value (z) for Chebyshev approximation.


        const double inf = highs->getInfinity();

        // ---- Add columns (variables) ----
        std::vector<double> col_lower(n, 0.0);
        std::vector<double> col_upper(n, inf);
        std::vector<double> row_lower(m, -inf);

        // Add all columns at once, with no matrix coefficients yet (we'll add rows next).
        // Signature: addCols(num_new_col, cost, lower, upper, num_nz, start, index, value)
        HighsStatus st = highs->addCols(
            n,
            c_vec.data(),
            col_lower.data(),
            col_upper.data(),
            /*num_nz=*/0,
            /*start=*/nullptr,
            /*index=*/nullptr,
            /*value=*/nullptr
        );
        if (st != HighsStatus::kOk) {
            throw std::runtime_error("InitializeHighs: highs.addCols failed.");
        }

        // ---- Add rows with their coefficients ----
        // Signature: addRows(num_new_row, lower, upper, num_nz, start, index, value)
        st = highs->addRows(
            m,
            row_lower.data(),
            row_upper.data(),
            value.size(),
            starts.data(),
            col_index.data(),
            value.data()
        );
        if (st != HighsStatus::kOk) {
            throw std::runtime_error("InitializeHighs: highs.addRows failed.");
        }

        return std::make_tuple<std::unique_ptr<Highs>, std::vector<double>, std::vector<double>>(std::move(highs), std::move(row_lower), std::move(row_upper));
    }

void GlobalAffineApproximator::updateHighsRhsUpperBounds(int phase, const std::vector<double>& v_prev_vec) {
        if (v_prev_vec.size() != kVDeltaDim) {
            throw std::invalid_argument("v_prev_vec size must be " + std::to_string(kVDeltaDim));
        }
        const auto& row_lower = row_lowers_[phase];
        const auto& row_upper = row_uppers_[phase];
        auto& highs_solver = highs_solvers_[phase];
        std::vector<double> new_row_upper(row_lower.size());
        const auto& intersection_points_phase = intersection_points_[phase];
        std::vector<double> v_prev_diff(kSpaceDim);
        std::vector<int> row_ids(row_upper.size());
        std::iota(row_ids.begin(), row_ids.end(), 0);

        for (size_t i = 0; i < kSpaceDim; ++i) {
            v_prev_diff[i] = v_prev_vec[i] - v_prev_vec[kSpaceDim + i]; // V_prev = V_prev^+ - V_prev^-
        }

        int j = 0;
        for (auto& area : intersection_points_phase) {
            for (auto& vertex : area) {
                double b_upd = 0.0;
                for (size_t i = 0; i < kSpaceDim; ++i) {
                    b_upd += v_prev_diff[i] * vertex(i);
                }
                b_upd += v_prev_vec[2 * kSpaceDim] - v_prev_vec[2 * kSpaceDim + 1]; // v_prev = v_prev^+ - v_prev^-
                new_row_upper[j] = row_upper[j] - b_upd; // [alpha_j, 0, 0] * x <= -tilde_beta_j - upd
                j++;
                new_row_upper[j] = row_upper[j] + b_upd; // [-alpha_j, 1, -1] * x <= tilde_beta_j + upd
                j++;
            }
        }

        auto st = highs_solver->changeRowsBounds(row_upper.size(), row_ids.data(), row_lower.data(), new_row_upper.data());
        if (st != HighsStatus::kOk) {
            throw std::runtime_error("updateHighsRhsUpperBounds: highs.changeRowsUpperBounds failed.");
        }
    }

std::vector<double> GlobalAffineApproximator::solveLp(int phase) {
        auto& highs_solver = highs_solvers_[phase];
        
        // Run the solver
        HighsStatus run_status = highs_solver->run();
        if (run_status != HighsStatus::kOk) {
            throw std::runtime_error("solveLp: highs_solver.run() failed with status " + std::to_string(static_cast<int>(run_status)));
        }
        
        // Check if solution is optimal
        if (highs_solver->getModelStatus() != HighsModelStatus::kOptimal) {
            throw std::runtime_error("solveLp: LP solution not found for phase " + std::to_string(phase) + ", model status: " + std::to_string(static_cast<int>(highs_solver->getModelStatus())));
        }
        
        // Extract solution
        const auto& solution = highs_solver->getSolution();
        std::vector<double> result(solution.col_value.begin(), solution.col_value.end());

        // x_next = [V^+, V^-, v^+, v^-, z^+, z^-], extract only v_next = [V^+, V^-, v^+, v^-]
        std::vector<double> v_next(kVDeltaDim);
        std::copy(result.begin(), result.begin() + kVDeltaDim, v_next.begin());
        double z_next = result[kVDeltaDim] - result[kVDeltaDim + 1];

        // z is a lower bound of the non-positive value, so if z apper to be non-negative, it means that something went wrong.
        if (z_next > -kEps) {
            throw std::runtime_error("solveLp: z_next is non-negative: " + std::to_string(z_next) + ", should be less than -kEps: " + std::to_string(-kEps));
        }
        
        return v_next;
    }

void GlobalAffineApproximator::precomputeMatrices() {
        logger_->info("Starting precomputeMatrices");
        if (intersection_points_.empty()) {
            throw std::runtime_error(
                "precomputeMatrices: intersection points not computed. Call getIntersectionPoints() before precomputeMatrices().");
        }
        if (!b_j_matrs_.empty() || !b_j_vecs_.empty() || !g_j_vecs_.empty() || !g_j_scals_.empty()) {
            throw std::runtime_error("precomputeMatrices: b_j_matrs, b_j_vecs, g_j_vecs, or g_j_scals are not empty.");
        }

        for (int phase = 0; phase < kPhases; ++phase) {
            auto [B_j_matrs, b_j_vecs, g_j_vecs, g_j_scals] = precomputeBMatAndBVec(phase);
            this->b_j_matrs_.push_back(std::move(B_j_matrs));
            this->b_j_vecs_.push_back(std::move(b_j_vecs));
            this->g_j_vecs_.push_back(std::move(g_j_vecs));
            this->g_j_scals_.push_back(std::move(g_j_scals));
        }
        logger_->info("Finished precomputeMatrices");
    }

void GlobalAffineApproximator::run() {
        logger_->info("Starting affine approximator");
        getIntersectionPoints();
        precomputeMatrices();

        // Initialize HiGHS solvers for each phase
        logger_->info("Start initializing Highs solvers");
        for (int phase = 0; phase < kPhases; ++phase) {
            auto [highs_solver, row_lower, row_upper] = initializeHighs(phase);
            highs_solvers_.push_back(std::move(highs_solver));
            row_lowers_.push_back(std::move(row_lower));
            row_uppers_.push_back(std::move(row_upper));
        }
        logger_->info("Done initializing Highs solvers");

        const double half_step = this->t_delta_ / 2;

        // Update HiGHS solvers for each switch count
        for (int switch_cnt = 0; switch_cnt <= max_switches_; ++switch_cnt) {
            logger_->info("Computing value function for switch count {}", switch_cnt);
            for (int phase = 0; phase < kPhases; ++phase) {
                logger_->info("Computing value function for phase {}", phase);

                double t_min, t_max;
                if (switch_cnt == 0) {
                    t_min = max(0.0, this->t_max_ - this->tau_max_);
                    t_max = this->t_max_;
                } else {
                    t_min = max(0.0, this->t_max_ - switch_cnt * this->tau_max_);
                    t_max = this->t_max_ - switch_cnt * this->tau_min_;
                }

                double theta_min = t_min;
                double theta_max = std::min(t_max + this->tau_max_, this->t_max_);
                
                // Replace NumPy-style boolean indexing with C++ loop
                std::vector<int> theta_range_ids;
                for (size_t i = 0; i < this->t_range_.size(); ++i) {
                    // t_theta_mina and theta might appear to be aligned with "t" grid nodes and that can cause instable behaviour due to numerical errors.
                    // Make a t_delta half-size step to ensure not a small gap between the grid nodes and possible same values representing segment borders
                    if (this->t_range_[i] >= (theta_min - half_step) && this->t_range_[i] <= (theta_max + half_step)) {
                        theta_range_ids.push_back(this->t_index_[i]);
                    }
                }
                
                std::vector<double> theta_range;
                theta_range.reserve(theta_range_ids.size());
                for (int idx : theta_range_ids) {
                    theta_range.push_back(this->t_range_[idx]);
                }

                const auto n = static_cast<std::ptrdiff_t>(theta_range.size());
                for (std::ptrdiff_t itheta_idx = n - 1; itheta_idx >= 0; --itheta_idx) {
                    int theta_idx = theta_range_ids[itheta_idx];
                    double theta = theta_range[itheta_idx];
                    logger_->info("Computing value function for theta {}", theta);
                    double t_theta_min = max(theta - this->tau_max_, t_min);
                    
                    // Replace NumPy-style boolean indexing with C++ loop
                    std::vector<int> t_theta_range_ids;
                    for (size_t i = 0; i < this->t_range_.size(); ++i) {
                        // t_theta_mina and theta might appear to be aligned with "t" grid nodes and that can cause instable behaviour due to numerical errors.
                        // Make a t_delta half-size step to ensure not a small gap between the grid nodes and possible same values representing segment borders
                        if (this->t_range_[i] >= (t_theta_min - half_step) && this->t_range_[i] < (theta - half_step)) {
                            t_theta_range_ids.push_back(this->t_index_[i]);
                        }
                    }
                    
                    int switch_phase = phase == 0 ? 1 : 0;
                    std::vector<double> v_prev = getBorderConditions(switch_phase, theta_idx, theta, switch_cnt);
                    // v_prev = [V^+, V^-, v^+, v^-], and x_prev = [V^+, V^-, v^+, v^-, z^+, z^-]
                    if (v_prev.size() != kVDeltaDim) {
                        throw std::invalid_argument("v_prev size must be " + std::to_string(kVDeltaDim));
                    }
                    setValueFunction(value_function_, phase, switch_cnt, theta_idx, theta_idx, v_prev);
                    updateHighsRhsUpperBounds(phase, v_prev);

                    if (t_theta_range_ids.empty()) {
                        throw std::runtime_error("t_theta_range_ids is empty: No previous timesteps to compute value function for.");
                    }

                    const auto n_t = static_cast<std::ptrdiff_t>(t_theta_range_ids.size());
                    for (std::ptrdiff_t i_t_idx = n_t - 1; i_t_idx >= 0; --i_t_idx) {
                        int t_idx = t_theta_range_ids[i_t_idx];
                        logger_->info("Computing value function for t {}", t_idx);
                        std::vector<double> v_next = solveLp(phase);
                        // v_next = [V^+, V^-, v^+, v^-]
                        if (v_next.size() != kVDeltaDim) {
                            throw std::invalid_argument("v_next size must be " + std::to_string(kVDeltaDim));
                        }
                        setValueFunction(value_function_, phase, switch_cnt, t_idx, theta_idx, v_next);
                        updateHighsRhsUpperBounds(phase, v_next);
                        v_prev = v_next;
                    }
                }
            }
        }
    }
}  // namespace global_affine_approximator