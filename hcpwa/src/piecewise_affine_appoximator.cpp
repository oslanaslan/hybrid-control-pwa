#include "piecewise_affine_approximator.h"

#include "util/assert_utils.hpp"

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
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Highs.h>
#include <spdlog/spdlog.h>
#include "spdlog/sinks/stdout_color_sinks.h"

namespace piecewise_affine_approximator {

const std::vector<int> kInIds = {2 - 1, 3 - 1, 5 - 1, 8 - 1};
const std::vector<int> kOutIds = {1 - 1, 4 - 1, 6 - 1, 7 - 1};
const int kPhases = 2;

VarLayout::VarLayout(const std::array<int, kSubsystemCount>& Ms)
    : M_s(Ms) {
    for (int s = 0; s < kSubsystemCount; ++s) {
        if (M_s[s] < 0) {
            throw std::invalid_argument("VarLayout: M_s must be non-negative");
        }
        numV += M_s[s] * kSpaceDim;
    }
    numX = numV + 1;
}

int VarLayout::idxV(int s, int k, int i) const {
    if (s < 0 || s >= kSubsystemCount) {
        throw std::invalid_argument("VarLayout::idxV: invalid subsystem id");
    }
    if (k < 0 || k >= M_s[s]) {
        throw std::invalid_argument(
            "VarLayout::idxV: invalid simplex id for subsystem");
    }
    if (i < 0 || i >= kSpaceDim) {
        throw std::invalid_argument(
            "VarLayout::idxV: invalid vector component id");
    }
    int offset = 0;
    for (int ss = 0; ss < s; ++ss) {
        offset += M_s[ss] * kSpaceDim;
    }
    return offset + k * kSpaceDim + i;
}

PiecewiseAffineApproximator::PiecewiseAffineApproximator(
    double t_max, int t_split_count, int max_switches, double tau_min,
    double tau_max, const SystemParams& system_params)
    : t_max_(t_max),
      t_split_count_(t_split_count),
      max_switches_(max_switches),
      tau_min_(tau_min),
      tau_max_(tau_max),
      system_params_(system_params) {
    if (t_split_count < 1) {
        throw std::invalid_argument(
            "t_split_count must be at least 1 in PiecewiseAffineApproximator "
            "constructor.");
    }
    if (tau_min >= tau_max || tau_min < 0.0 || tau_max < 0.0) {
        throw std::invalid_argument(
            "tau_min must be less than tau_max and both must be >= 0 in "
            "PiecewiseAffineApproximator constructor.");
    }
    if (max_switches < 1) {
        throw std::invalid_argument(
            "max_switches must be at least 1 in PiecewiseAffineApproximator "
            "constructor.");
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
    const int n_vertices = 1 << kSpaceDim;
    cube_angle_vertices_.clear();
    for (int vert = 0; vert < n_vertices; ++vert) {
        Eigen::VectorXd v(kSpaceDim);
        for (int d = 0; d < kSpaceDim; ++d) {
            v(d) = ((vert & (1 << d)) != 0) ? this->system_params_.N : 0.0;
        }
        cube_angle_vertices_.push_back(v);
    }
    logger_ = spdlog::stdout_color_mt("piecewise_affine_approximator");
    logger_->set_level(spdlog::level::info);
}

std::tuple<double, double, double, double>
PiecewiseAffineApproximator::getMainSystemParams() const {
    return std::make_tuple(system_params_.N, system_params_.F, system_params_.v,
                           system_params_.w);
}

double PiecewiseAffineApproximator::getBetaParamForAxis(int i, int j) const {
    if (i == 5 - 1 && j == 1 - 1) {
        return system_params_.b51;
    }
    if (i == 5 - 1 && j == 7 - 1) {
        return system_params_.b57;
    }
    if (i == 8 - 1 && j == 4 - 1) {
        return system_params_.b84;
    }
    if (i == 8 - 1 && j == 6 - 1) {
        return system_params_.b86;
    }
    if (i == 3 - 1 && j == 1 - 1) {
        return system_params_.b31;
    }
    if (i == 3 - 1 && j == 6 - 1) {
        return system_params_.b36;
    }
    if (i == 2 - 1 && j == 4 - 1) {
        return system_params_.b24;
    }
    if (i == 2 - 1 && j == 7 - 1) {
        return system_params_.b27;
    }
    throw std::invalid_argument(
        std::format("No such beta params for axis ({}, {})", i, j));
}

std::pair<double, double> PiecewiseAffineApproximator::getFMinMaxForAxis(
    int i) const {
    if (i == 2 - 1) {
        return std::make_pair(system_params_.f2min, system_params_.f2max);
    }
    if (i == 3 - 1) {
        return std::make_pair(system_params_.f3min, system_params_.f3max);
    }
    if (i == 5 - 1) {
        return std::make_pair(system_params_.f5min, system_params_.f5max);
    }
    if (i == 8 - 1) {
        return std::make_pair(system_params_.f8min, system_params_.f8max);
    }
    throw std::invalid_argument("No such f min max for axis "
                                + std::to_string(i));
}

void PiecewiseAffineApproximator::getIntersectionPoints() {
    hcpwa::TriangleAreasVerticesResult areas_vertices = hcpwa::compute_triangle_areas_vertices(
        system_params_.N, system_params_.F, system_params_.v, system_params_.w,
        system_params_.b51, system_params_.b57, system_params_.b84,
        system_params_.b86, system_params_.b31, system_params_.b36,
        system_params_.b24, system_params_.b27, system_params_.f2min,
        system_params_.f3min, system_params_.f5min, system_params_.f8min,
        system_params_.f2max, system_params_.f3max, system_params_.f5max,
        system_params_.f8max);

    intersection_points_.resize(kPhases);
    prism_indices_.resize(kPhases);
    subsystem_axes_.resize(kPhases);
    layouts_.clear();
    layouts_.reserve(kPhases);

    auto convert_phase = [](const std::vector<std::vector<hcpwa::Vec<8>>>& src,
                            std::vector<std::vector<Eigen::VectorXd>>& dst) {
        dst.resize(src.size());
        for (std::size_t a = 0; a < src.size(); ++a) {
            const auto& area = src[a];
            dst[a].resize(area.size());
            for (std::size_t v = 0; v < area.size(); ++v) {
                Eigen::VectorXd vec(kSpaceDim);
                for (int i = 0; i < kSpaceDim; ++i) {
                    vec(i) = static_cast<double>(area[v][i]);
                }
                dst[a][v] = std::move(vec);
            }
        }
    };

    convert_phase(areas_vertices.intersection_points_phase0,
                  intersection_points_[0]);
    convert_phase(areas_vertices.intersection_points_phase1,
                  intersection_points_[1]);

    prism_indices_[0] = areas_vertices.intersection_prism_indices_phase0;
    prism_indices_[1] = areas_vertices.intersection_prism_indices_phase1;
    for (int phase = 0; phase < kPhases; ++phase) {
        for (const auto& region_idx : prism_indices_[phase]) {
            if (region_idx.size() != kSubsystemCount) {
                throw std::runtime_error(
                    "Expected 5 simplex indices per region");
            }
        }
    }

    subsystem_axes_[0] = {{{2, 0}, {2, 5}, {1, 3}, {1, 6}, {4, 7}}};
    subsystem_axes_[1] = {{{4, 0}, {4, 6}, {7, 3}, {7, 5}, {1, 2}}};

    layouts_.emplace_back(std::array<int, kSubsystemCount>{
        static_cast<int>(areas_vertices.triangles31.size()),
        static_cast<int>(areas_vertices.triangles36.size()),
        static_cast<int>(areas_vertices.triangles24.size()),
        static_cast<int>(areas_vertices.triangles27.size()),
        static_cast<int>(areas_vertices.triangles58.size()),
    });
    layouts_.emplace_back(std::array<int, kSubsystemCount>{
        static_cast<int>(areas_vertices.triangles51.size()),
        static_cast<int>(areas_vertices.triangles57.size()),
        static_cast<int>(areas_vertices.triangles84.size()),
        static_cast<int>(areas_vertices.triangles86.size()),
        static_cast<int>(areas_vertices.triangles23.size()),
    });
}

std::pair<Eigen::RowVectorXd, Eigen::RowVectorXd>
PiecewiseAffineApproximator::getFIJMinResolution(
    int i, int j, const Eigen::VectorXd& n) const {
    const double N = system_params_.N;
    const double F = system_params_.F;
    const double v = system_params_.v;
    const double w = system_params_.w;
    const double beta_i_j = getBetaParamForAxis(i, j);
    const double n_i = n(i);
    const double n_j = n(j);
    Eigen::RowVectorXd f_matr_row = Eigen::RowVectorXd::Zero(kSpaceDim);
    Eigen::RowVectorXd f_vec_row = Eigen::RowVectorXd::Zero(1);

    const double a = beta_i_j * v * n_i;
    const double b = beta_i_j * F;
    const double c = w * (N - n_j);

    if (b < a + kEps && b < c + kEps) {
        f_vec_row(0) = beta_i_j * F;
    } else if (c < a + kEps && c < b + kEps) {
        f_matr_row(j) = -w;
        f_vec_row(0) = w * N;
    } else if (a < b + kEps && a < c + kEps) {
        f_matr_row(i) = beta_i_j * v;
    } else {
        throw std::invalid_argument(std::format(
            "No such case: a={}, b={}, c={} (getFIJMinResolution)", a, b, c));
    }
    return std::make_pair(f_matr_row, f_vec_row);
}

Eigen::VectorXd PiecewiseAffineApproximator::areaCentroidCoords(
    int j, int phase) const {
    if (intersection_points_.empty() || phase < 0
        || phase >= static_cast<int>(intersection_points_.size())
        || intersection_points_[phase].empty()) {
        throw std::runtime_error(
            "Intersection points are not computed for the specified phase");
    }
    const auto& intersection_points_per_phase = intersection_points_[phase];
    if (j < 0 || j >= static_cast<int>(intersection_points_per_phase.size())) {
        throw std::invalid_argument(
            "Invalid area index j in areaCentroidCoords");
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
    hcpwa::util::assertShape(centroid, kSpaceDim);
    return centroid;
}

std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, double>
PiecewiseAffineApproximator::getAMatrFVecGVecAndGScalJ(int j, int phase) const {
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

        a_matr.resize(kSpaceDim, kSpaceDim);
        a_matr.row(0) = f31_matr_row;
        a_matr.row(1) = -f24_matr_row - f27_matr_row;
        a_matr.row(2) = -f31_matr_row - f36_matr_row;
        a_matr.row(3) = f24_matr_row;
        a_matr.row(4) = Eigen::RowVectorXd::Zero(kSpaceDim);
        a_matr.row(5) = f36_matr_row;
        a_matr.row(6) = f27_matr_row;
        a_matr.row(7) = Eigen::RowVectorXd::Zero(kSpaceDim);

        b_vec.resize(kSpaceDim);
        b_vec(0) = f31_vec_row(0);
        b_vec(1) = -f24_vec_row(0) - f27_vec_row(0);
        b_vec(2) = -f31_vec_row(0) - f36_vec_row(0);
        b_vec(3) = f24_vec_row(0);
        b_vec(4) = 0.0;
        b_vec(5) = f36_vec_row(0);
        b_vec(6) = f27_vec_row(0);
        b_vec(7) = 0.0;

        g_vec = (f31_matr_row + f36_matr_row + f24_matr_row + f27_matr_row)
                    .transpose();
        g_scal
            = f31_vec_row(0) + f36_vec_row(0) + f24_vec_row(0) + f27_vec_row(0);
    } else if (phase == 1) {
        auto [f51_matr_row, f51_vec_row] = getFIJMinResolution(5 - 1, 1 - 1, n);
        auto [f57_matr_row, f57_vec_row] = getFIJMinResolution(5 - 1, 7 - 1, n);
        auto [f84_matr_row, f84_vec_row] = getFIJMinResolution(8 - 1, 4 - 1, n);
        auto [f86_matr_row, f86_vec_row] = getFIJMinResolution(8 - 1, 6 - 1, n);

        a_matr.resize(kSpaceDim, kSpaceDim);
        a_matr.row(0) = f51_matr_row;
        a_matr.row(1) = Eigen::RowVectorXd::Zero(kSpaceDim);
        a_matr.row(2) = Eigen::RowVectorXd::Zero(kSpaceDim);
        a_matr.row(3) = f84_matr_row;
        a_matr.row(4) = -f51_matr_row - f57_matr_row;
        a_matr.row(5) = f86_matr_row;
        a_matr.row(6) = f57_matr_row;
        a_matr.row(7) = -f84_matr_row - f86_matr_row;

        b_vec.resize(kSpaceDim);
        b_vec(0) = f51_vec_row(0);
        b_vec(1) = 0.0;
        b_vec(2) = 0.0;
        b_vec(3) = f84_vec_row(0);
        b_vec(4) = -f51_vec_row(0) - f57_vec_row(0);
        b_vec(5) = f86_vec_row(0);
        b_vec(6) = f57_vec_row(0);
        b_vec(7) = -f84_vec_row(0) - f86_vec_row(0);

        g_vec = (f51_matr_row + f57_matr_row + f84_matr_row + f86_matr_row)
                    .transpose();
        g_scal
            = f51_vec_row(0) + f57_vec_row(0) + f84_vec_row(0) + f86_vec_row(0);
    } else {
        throw std::invalid_argument(std::format("No such phase: {}", phase));
    }

    hcpwa::util::assertShape(b_vec, kSpaceDim);
    hcpwa::util::assertShape(a_matr, kSpaceDim, kSpaceDim);
    hcpwa::util::assertShape(g_vec, kSpaceDim);
    hcpwa::util::assertScalar(g_scal);

    return std::make_tuple(a_matr, b_vec, g_vec, g_scal);
}

std::tuple<std::vector<int>, std::vector<int>, std::vector<double>,
           std::vector<double>, Eigen::RowVectorXd>
PiecewiseAffineApproximator::prepareLpMatrices(int phase) {
    const auto& layout = layouts_.at(phase);
    const auto& intersection_points_phase = intersection_points_.at(phase);
    const auto& prism_indices_phase = prism_indices_.at(phase);
    const auto& subsystem_axes = subsystem_axes_.at(phase);
    const auto& a_j_matrs = a_j_matrs_.at(phase);
    const auto& f_j_vecs = f_j_vecs_.at(phase);
    const auto& g_j_vecs = g_j_vecs_.at(phase);
    const auto& g_j_scals = g_j_scals_.at(phase);

    if (intersection_points_phase.size() != prism_indices_phase.size()) {
        throw std::runtime_error(
            "prepareLpMatrices: inconsistent geometry and prism index counts");
    }

    std::vector<int> starts;
    std::vector<int> col_index;
    std::vector<double> value;
    std::vector<double> row_upper;
    std::vector<double> b0_upper;
    Eigen::RowVectorXd c_vec = Eigen::RowVectorXd::Zero(layout.numX);
    c_vec(layout.idxS()) = 1.0;

    beta_terms_[phase].clear();
    std::size_t nnz = 0;

    for (std::size_t j = 0; j < intersection_points_phase.size(); ++j) {
        const auto& Aj = a_j_matrs[j];
        const auto& fj = f_j_vecs[j];
        const auto& gj = g_j_vecs[j];
        const double g_scal = g_j_scals[j];
        const auto& j_s_size_t = prism_indices_phase[j];
        if (j_s_size_t.size() != kSubsystemCount) {
            throw std::runtime_error(
                "prepareLpMatrices: each region must have 5 simplex indices");
        }

        std::array<int, kSubsystemCount> j_s{};
        for (int s = 0; s < kSubsystemCount; ++s) {
            j_s[s] = static_cast<int>(j_s_size_t[s]);
            if (j_s[s] < 0 || j_s[s] >= layout.M_s[s]) {
                throw std::runtime_error(
                    "prepareLpMatrices: simplex index out of range for "
                    "subsystem");
            }
        }

        for (const auto& nu : intersection_points_phase[j]) {
            Eigen::VectorXd z = Aj * nu + fj;
            const double g_nu = gj.dot(nu) + g_scal;

            std::array<Eigen::Matrix<double, kSpaceDim, 1>, kSubsystemCount>
                q_arr{};
            std::array<Eigen::Matrix<double, kSpaceDim, 1>, kSubsystemCount>
                r_arr{};
            for (int s = 0; s < kSubsystemCount; ++s) {
                q_arr[s].setZero();
                r_arr[s].setZero();
                const int axis_a = subsystem_axes[s][0];
                const int axis_b = subsystem_axes[s][1];
                q_arr[s](axis_a) = nu(axis_a);
                q_arr[s](axis_b) = nu(axis_b);
                r_arr[s](axis_a) = z(axis_a) - q_arr[s](axis_a) / t_delta_;
                r_arr[s](axis_b) = z(axis_b) - q_arr[s](axis_b) / t_delta_;
            }

            // Row (I): -sum_s (v_m^{s, j_s})^T r_s <= g + beta
            starts.push_back(static_cast<int>(nnz));
            for (int s = 0; s < kSubsystemCount; ++s) {
                const int k = j_s[s];
                for (int i = 0; i < kSpaceDim; ++i) {
                    const double coeff = -r_arr[s](i);
                    if (std::abs(coeff) <= kEps) {
                        continue;
                    }
                    col_index.push_back(layout.idxV(s, k, i));
                    value.push_back(coeff);
                    ++nnz;
                }
            }
            row_upper.push_back(g_nu);
            b0_upper.push_back(g_nu);
            beta_terms_[phase].push_back(RHSBetaTerm{
                static_cast<int>(row_upper.size() - 1), j_s, q_arr, +1.0});

            // Row (II): sum_s (v_m^{s, j_s})^T r_s - s <= -g - beta
            starts.push_back(static_cast<int>(nnz));
            for (int s = 0; s < kSubsystemCount; ++s) {
                const int k = j_s[s];
                for (int i = 0; i < kSpaceDim; ++i) {
                    const double coeff = r_arr[s](i);
                    if (std::abs(coeff) <= kEps) {
                        continue;
                    }
                    col_index.push_back(layout.idxV(s, k, i));
                    value.push_back(coeff);
                    ++nnz;
                }
            }
            col_index.push_back(layout.idxS());
            value.push_back(-1.0);
            ++nnz;
            row_upper.push_back(-g_nu);
            b0_upper.push_back(-g_nu);
            beta_terms_[phase].push_back(RHSBetaTerm{
                static_cast<int>(row_upper.size() - 1), j_s, q_arr, -1.0});
        }
    }
    starts.push_back(static_cast<int>(nnz));
    b0_uppers_[phase] = b0_upper;
    return std::make_tuple(std::move(starts), std::move(col_index),
                           std::move(value), std::move(row_upper),
                           std::move(c_vec));
}

double PiecewiseAffineApproximator::getBorderFuncValuesAtN(
    int r, int theta_idx, int theta_end_idx, int phase,
    const Eigen::VectorXd& n) {
    try {
        const auto x = value_function_.get(phase, r, theta_idx, theta_end_idx);
        const auto& layout = layouts_.at(phase);
        if (x.size() != static_cast<std::size_t>(layout.numV)) {
            throw std::invalid_argument(
                "ValueFunction vector has invalid size");
        }
        const auto& prism_indices_phase = prism_indices_.at(phase);
        const auto& subsystem_axes = subsystem_axes_.at(phase);
        double best = -std::numeric_limits<double>::infinity();
        for (const auto& region : prism_indices_phase) {
            if (region.size() != kSubsystemCount) {
                continue;
            }
            double value_at_region = 0.0;
            for (int s = 0; s < kSubsystemCount; ++s) {
                const int k = static_cast<int>(region[s]);
                const int axis_a = subsystem_axes[s][0];
                const int axis_b = subsystem_axes[s][1];
                value_at_region += x[layout.idxV(s, k, axis_a)] * n(axis_a);
                value_at_region += x[layout.idxV(s, k, axis_b)] * n(axis_b);
            }
            best = std::max(best, value_at_region);
        }
        return best;
    } catch (const std::invalid_argument&) {
        return -std::numeric_limits<double>::infinity();
    }
}

double PiecewiseAffineApproximator::getMaxBorderFuncValuesAtN(
    int theta_idx, const std::vector<int>& theta_end_ids, int max_switches,
    int phase, const Eigen::VectorXd& n) {
    double max_val = -std::numeric_limits<double>::infinity();
    for (int theta_end_idx : theta_end_ids) {
        for (int r = 0; r < max_switches; ++r) {
            const double val
                = getBorderFuncValuesAtN(r, theta_idx, theta_end_idx, phase, n);
            max_val = std::max(max_val, val);
        }
    }
    return max_val;
}

std::vector<double> PiecewiseAffineApproximator::getBorderConditions(
    int switch_phase, int theta_idx, double theta, int switch_cnt) {
    const auto& layout = layouts_.at(switch_phase);
    if (switch_cnt == 0) {
        return std::vector<double>(layout.numV, 0.0);
    }

    double theta_min = std::min(theta + this->tau_min_, this->t_max_);
    double theta_max = std::min(theta + this->tau_max_, this->t_max_);
    const double half_step = this->t_delta_ / 2.0;
    theta_min -= half_step;
    theta_max += half_step;

    std::vector<int> theta_range_ids;
    for (std::size_t i = 0; i < this->t_range_.size(); ++i) {
        if (this->t_range_[i] >= theta_min && this->t_range_[i] <= theta_max) {
            theta_range_ids.push_back(this->t_index_[i]);
        }
    }
    if (theta_range_ids.empty()) {
        throw std::runtime_error("getBorderConditions: No theta range found");
    }

    std::vector<int> starts;
    std::vector<int> col_index;
    std::vector<double> value;
    std::vector<double> row_upper;
    Eigen::RowVectorXd c_vec = Eigen::RowVectorXd::Zero(layout.numV);
    std::size_t nnz = 0;

    const auto& prism_indices_phase = prism_indices_.at(switch_phase);
    const auto& subsystem_axes = subsystem_axes_.at(switch_phase);
    for (const auto& vertex : this->cube_angle_vertices_) {
        const double f_scal = getMaxBorderFuncValuesAtN(
            theta_idx, theta_range_ids, switch_cnt, switch_phase, vertex);
        for (const auto& region : prism_indices_phase) {
            starts.push_back(static_cast<int>(nnz));
            for (int s = 0; s < kSubsystemCount; ++s) {
                const int k = static_cast<int>(region[s]);
                const int axis_a = subsystem_axes[s][0];
                const int axis_b = subsystem_axes[s][1];
                const int col_a = layout.idxV(s, k, axis_a);
                const int col_b = layout.idxV(s, k, axis_b);
                const double coeff_a = -vertex(axis_a);
                const double coeff_b = -vertex(axis_b);
                if (std::abs(coeff_a) > kEps) {
                    col_index.push_back(col_a);
                    value.push_back(coeff_a);
                    ++nnz;
                    c_vec(col_a) += vertex(axis_a);
                }
                if (std::abs(coeff_b) > kEps) {
                    col_index.push_back(col_b);
                    value.push_back(coeff_b);
                    ++nnz;
                    c_vec(col_b) += vertex(axis_b);
                }
            }
            row_upper.push_back(-f_scal);
        }
    }
    starts.push_back(static_cast<int>(nnz));

    const int m = static_cast<int>(row_upper.size());
    const int n = layout.numV;
    Highs highs;
    highs.setOptionValue("solver", "simplex");
    highs.setOptionValue("presolve", "on");
    highs.changeObjectiveSense(ObjSense::kMaximize);
    const double inf = highs.getInfinity();
    const double border = std::max(1.0, system_params_.N);
    std::vector<double> col_lower(n, -border);
    std::vector<double> col_upper(n, border);
    std::vector<double> row_lower(m, -inf);

    HighsStatus st
        = highs.addCols(n, c_vec.data(), col_lower.data(), col_upper.data(), 0,
                        nullptr, nullptr, nullptr);
    if (st != HighsStatus::kOk) {
        throw std::runtime_error("getBorderConditions: highs.addCols failed");
    }

    st = highs.addRows(m, row_lower.data(), row_upper.data(),
                       static_cast<int>(value.size()), starts.data(),
                       col_index.data(), value.data());
    if (st != HighsStatus::kOk) {
        throw std::runtime_error("getBorderConditions: highs.addRows failed");
    }

    st = highs.run();
    if (st != HighsStatus::kOk
        || highs.getModelStatus() != HighsModelStatus::kOptimal) {
        throw std::runtime_error("getBorderConditions: LP solution not found");
    }
    const auto& solution = highs.getSolution();
    return std::vector<double>(solution.col_value.begin(),
                               solution.col_value.end());
}

std::tuple<std::unique_ptr<Highs>, std::vector<double>, std::vector<double>>
PiecewiseAffineApproximator::initializeHighs(int phase) {
    auto [starts, col_index, value, row_upper, c_vec]
        = prepareLpMatrices(phase);
    const int m = static_cast<int>(row_upper.size());
    const int n = static_cast<int>(c_vec.size());
    const auto& layout = layouts_.at(phase);

    std::unique_ptr<Highs> highs = std::make_unique<Highs>();
    highs->setOptionValue("solver", "simplex");
    highs->setOptionValue("presolve", "on");
    highs->changeObjectiveSense(ObjSense::kMinimize);
    const double inf = highs->getInfinity();

    std::vector<double> col_lower(n, -inf);
    std::vector<double> col_upper(n, inf);
    col_lower[layout.idxS()] = 0.0;
    std::vector<double> row_lower(m, -inf);

    HighsStatus st
        = highs->addCols(n, c_vec.data(), col_lower.data(), col_upper.data(), 0,
                         nullptr, nullptr, nullptr);
    if (st != HighsStatus::kOk) {
        throw std::runtime_error("initializeHighs: highs.addCols failed.");
    }

    st = highs->addRows(m, row_lower.data(), row_upper.data(),
                        static_cast<int>(value.size()), starts.data(),
                        col_index.data(), value.data());
    if (st != HighsStatus::kOk) {
        throw std::runtime_error("initializeHighs: highs.addRows failed.");
    }

    return std::make_tuple(std::move(highs), std::move(row_lower),
                           std::move(row_upper));
}

void PiecewiseAffineApproximator::updateHighsRhsUpperBounds(
    int phase, const std::vector<double>& v_prev_vec) {
    const auto& layout = layouts_.at(phase);
    if (v_prev_vec.size() != static_cast<std::size_t>(layout.numV)) {
        throw std::invalid_argument("v_prev_vec size must be "
                                    + std::to_string(layout.numV));
    }
    const auto& row_lower = row_lowers_.at(phase);
    const auto& b0_upper = b0_uppers_.at(phase);
    const auto& terms = beta_terms_.at(phase);
    auto& highs_solver = highs_solvers_.at(phase);
    if (row_lower.size() != b0_upper.size()
        || row_lower.size() != terms.size()) {
        throw std::runtime_error(
            "updateHighsRhsUpperBounds: invalid precomputed RHS buffers");
    }

    std::vector<double> new_row_upper = b0_upper;
    for (const auto& term : terms) {
        double beta = 0.0;
        for (int s = 0; s < kSubsystemCount; ++s) {
            const int k = term.j_s[s];
            for (int i = 0; i < kSpaceDim; ++i) {
                beta += v_prev_vec[layout.idxV(s, k, i)] * term.q[s](i)
                        / t_delta_;
            }
        }
        new_row_upper[term.row] += term.sign * beta;
    }

    std::vector<int> row_ids(new_row_upper.size());
    std::iota(row_ids.begin(), row_ids.end(), 0);
    auto st = highs_solver->changeRowsBounds(static_cast<int>(row_ids.size()),
                                             row_ids.data(), row_lower.data(),
                                             new_row_upper.data());
    if (st != HighsStatus::kOk) {
        throw std::runtime_error(
            "updateHighsRhsUpperBounds: highs.changeRowsBounds failed.");
    }
}

std::vector<double> PiecewiseAffineApproximator::solveLp(int phase) {
    auto& highs_solver = highs_solvers_.at(phase);
    const auto& layout = layouts_.at(phase);
    HighsStatus run_status = highs_solver->run();
    if (run_status != HighsStatus::kOk) {
        throw std::runtime_error(
            "solveLp: highs_solver.run() failed with status "
            + std::to_string(static_cast<int>(run_status)));
    }
    if (highs_solver->getModelStatus() != HighsModelStatus::kOptimal) {
        throw std::runtime_error(
            "solveLp: LP solution not found for phase " + std::to_string(phase)
            + ", model status: "
            + std::to_string(static_cast<int>(highs_solver->getModelStatus())));
    }

    const auto& solution = highs_solver->getSolution();
    std::vector<double> result(solution.col_value.begin(),
                               solution.col_value.end());
    std::vector<double> v_next(layout.numV);
    std::copy(result.begin(), result.begin() + layout.numV, v_next.begin());
    return v_next;
}

void PiecewiseAffineApproximator::precomputeMatrices() {
    logger_->info("Starting precomputeMatrices");
    if (intersection_points_.empty() || prism_indices_.empty()
        || layouts_.size() != kPhases) {
        throw std::runtime_error(
            "precomputeMatrices: intersection points are not initialized");
    }
    a_j_matrs_.assign(kPhases, {});
    f_j_vecs_.assign(kPhases, {});
    g_j_vecs_.assign(kPhases, {});
    g_j_scals_.assign(kPhases, {});
    beta_terms_.assign(kPhases, {});
    b0_uppers_.assign(kPhases, {});

    for (int phase = 0; phase < kPhases; ++phase) {
        const int n_areas
            = static_cast<int>(intersection_points_[phase].size());
        a_j_matrs_[phase].reserve(n_areas);
        f_j_vecs_[phase].reserve(n_areas);
        g_j_vecs_[phase].reserve(n_areas);
        g_j_scals_[phase].reserve(n_areas);
        for (int j = 0; j < n_areas; ++j) {
            auto [A_j, f_j, g_j_vec, g_j_scal]
                = getAMatrFVecGVecAndGScalJ(j, phase);
            a_j_matrs_[phase].emplace_back(std::move(A_j));
            f_j_vecs_[phase].emplace_back(std::move(f_j));
            g_j_vecs_[phase].emplace_back(std::move(g_j_vec));
            g_j_scals_[phase].emplace_back(g_j_scal);
        }
    }
    logger_->info("Finished precomputeMatrices");
}

void PiecewiseAffineApproximator::run() {
    logger_->info("Starting piecewise affine approximator");
    getIntersectionPoints();
    precomputeMatrices();

    highs_solvers_.clear();
    row_lowers_.clear();
    row_uppers_.clear();
    highs_solvers_.reserve(kPhases);
    row_lowers_.reserve(kPhases);
    row_uppers_.reserve(kPhases);

    logger_->info("Start initializing Highs solvers");
    for (int phase = 0; phase < kPhases; ++phase) {
        auto [highs_solver, row_lower, row_upper] = initializeHighs(phase);
        highs_solvers_.push_back(std::move(highs_solver));
        row_lowers_.push_back(std::move(row_lower));
        row_uppers_.push_back(std::move(row_upper));
    }
    logger_->info("Done initializing Highs solvers");

    const double half_step = this->t_delta_ / 2;
    for (int switch_cnt = 0; switch_cnt <= max_switches_; ++switch_cnt) {
        logger_->info("Computing value function for switch count {}",
                      switch_cnt);
        for (int phase = 0; phase < kPhases; ++phase) {
            logger_->info("Computing value function for phase {}", phase);
            const auto& layout = layouts_[phase];

            double t_min;
            double t_max;
            if (switch_cnt == 0) {
                t_min = std::max(0.0, this->t_max_ - this->tau_max_);
                t_max = this->t_max_;
            } else {
                t_min
                    = std::max(0.0, this->t_max_ - switch_cnt * this->tau_max_);
                t_max = this->t_max_ - switch_cnt * this->tau_min_;
            }

            const double theta_min = t_min;
            const double theta_max
                = std::min(t_max + this->tau_max_, this->t_max_);
            std::vector<int> theta_range_ids;
            for (std::size_t i = 0; i < this->t_range_.size(); ++i) {
                if (this->t_range_[i] >= (theta_min - half_step)
                    && this->t_range_[i] <= (theta_max + half_step)) {
                    theta_range_ids.push_back(this->t_index_[i]);
                }
            }
            std::vector<double> theta_range;
            theta_range.reserve(theta_range_ids.size());
            for (int idx : theta_range_ids) {
                theta_range.push_back(this->t_range_[idx]);
            }

            const auto n = static_cast<std::ptrdiff_t>(theta_range.size());
            for (std::ptrdiff_t itheta_idx = n - 1; itheta_idx >= 0;
                 --itheta_idx) {
                const int theta_idx = theta_range_ids[itheta_idx];
                const double theta = theta_range[itheta_idx];
                logger_->info("Computing value function for theta {}", theta);
                const double t_theta_min
                    = std::max(theta - this->tau_max_, t_min);

                std::vector<int> t_theta_range_ids;
                for (std::size_t i = 0; i < this->t_range_.size(); ++i) {
                    if (this->t_range_[i] >= (t_theta_min - half_step)
                        && this->t_range_[i] < (theta - half_step)) {
                        t_theta_range_ids.push_back(this->t_index_[i]);
                    }
                }

                const int switch_phase = phase == 0 ? 1 : 0;
                std::vector<double> v_prev = getBorderConditions(
                    switch_phase, theta_idx, theta, switch_cnt);
                if (v_prev.size() != static_cast<std::size_t>(layout.numV)) {
                    throw std::invalid_argument("v_prev size must be "
                                                + std::to_string(layout.numV));
                }
                value_function_.set(phase, switch_cnt, theta_idx, theta_idx,
                                    v_prev, layout.numV);
                updateHighsRhsUpperBounds(phase, v_prev);

                if (t_theta_range_ids.empty()) {
                    throw std::runtime_error(
                        "t_theta_range_ids is empty: No previous timesteps to "
                        "compute "
                        "value function for.");
                }

                const auto n_t
                    = static_cast<std::ptrdiff_t>(t_theta_range_ids.size());
                for (std::ptrdiff_t i_t_idx = n_t - 1; i_t_idx >= 0;
                     --i_t_idx) {
                    const int t_idx = t_theta_range_ids[i_t_idx];
                    logger_->info("Computing value function for t {}", t_idx);
                    std::vector<double> v_next = solveLp(phase);
                    if (v_next.size()
                        != static_cast<std::size_t>(layout.numV)) {
                        throw std::invalid_argument(
                            "v_next size must be "
                            + std::to_string(layout.numV));
                    }
                    value_function_.set(phase, switch_cnt, t_idx, theta_idx,
                                        v_next, layout.numV);
                    updateHighsRhsUpperBounds(phase, v_next);
                    v_prev = v_next;
                }
            }
        }
    }
}

}  // namespace piecewise_affine_approximator