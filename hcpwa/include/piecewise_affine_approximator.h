#ifndef HCPWA_PIECEWISE_AFFINE_APPROXIMATOR_H
#define HCPWA_PIECEWISE_AFFINE_APPROXIMATOR_H

#include <array>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Highs.h>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <spdlog/spdlog.h>

#include "util/assert_utils.hpp"
#include "util/value_function_utils.hpp"

namespace piecewise_affine_approximator {

// Constants
extern const std::vector<int> kInIds;
extern const std::vector<int> kOutIds;
extern const int kPhases;
constexpr int kSpaceDim = 8;
constexpr int kSubsystemCount = 5;
constexpr double kEps = 1e-6;

using ValueFunction = hcpwa::util::ValueFunction;

struct VarLayout {
    std::array<int, kSubsystemCount> M_s{};
    int numV = 0;
    int numX = 1;

    VarLayout() = default;
    explicit VarLayout(const std::array<int, kSubsystemCount>& Ms);

    int idxV(int s, int k, int i) const;
    int idxS() const { return numX - 1; }
};

struct RHSBetaTerm {
    int row = 0;
    std::array<int, kSubsystemCount> j_s{};
    std::array<Eigen::Matrix<double, kSpaceDim, 1>, kSubsystemCount> q{};
    double sign = 0.0;
};

struct SystemParams {
    double N;
    double F;
    double v;
    double w;
    // axis ids in paper start from 1, but in code start from 0
    double b51;
    double b57;
    double b84;
    double b86;
    double b31;
    double b36;
    double b24;
    double b27;
    double f2min;
    double f3min;
    double f5min;
    double f8min;
    double f2max;
    double f3max;
    double f5max;
    double f8max;
};

class PiecewiseAffineApproximator {
private:
    double t_max_;
    int t_split_count_;
    int max_switches_;
    double tau_min_;
    double tau_max_;
    double t_delta_;
    SystemParams system_params_;
    std::vector<double> t_range_;
    std::vector<int> t_index_;

    ValueFunction value_function_;
    std::vector<Eigen::VectorXd> cube_angle_vertices_;
    std::vector<std::vector<std::vector<Eigen::VectorXd>>> intersection_points_;

    std::vector<VarLayout> layouts_;
    std::vector<std::array<std::array<int, 2>, kSubsystemCount>> subsystem_axes_;
    std::vector<std::vector<std::vector<std::size_t>>> prism_indices_;

    std::vector<std::vector<Eigen::MatrixXd>> a_j_matrs_;
    std::vector<std::vector<Eigen::VectorXd>> f_j_vecs_;
    std::vector<std::vector<Eigen::VectorXd>> g_j_vecs_;
    std::vector<std::vector<double>> g_j_scals_;

    std::vector<std::vector<RHSBetaTerm>> beta_terms_;

    std::shared_ptr<spdlog::logger> logger_;
    std::vector<std::unique_ptr<Highs>> highs_solvers_;
    std::vector<std::vector<double>> row_lowers_;
    std::vector<std::vector<double>> row_uppers_;
    std::vector<std::vector<double>> b0_uppers_;

public:
    PiecewiseAffineApproximator(double t_max,
                                int t_split_count,
                                int max_switches,
                                double tau_min,
                                double tau_max,
                                const SystemParams& system_params);

    std::tuple<double, double, double, double> getMainSystemParams() const;

    double getBetaParamForAxis(int i, int j) const;

    std::pair<double, double> getFMinMaxForAxis(int i) const;

    void getIntersectionPoints();

    std::pair<Eigen::RowVectorXd, Eigen::RowVectorXd> getFIJMinResolution(
        int i, int j, const Eigen::VectorXd& n) const;

    Eigen::VectorXd areaCentroidCoords(int j, int phase) const;

    std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, double>
    getAMatrFVecGVecAndGScalJ(int j, int phase) const;

    std::tuple<std::vector<int>,
               std::vector<int>,
               std::vector<double>,
               std::vector<double>,
               Eigen::RowVectorXd>
    prepareLpMatrices(int phase);

    double getBorderFuncValuesAtN(int r,
                                  int theta_idx,
                                  int theta_end_idx,
                                  int phase,
                                  const Eigen::VectorXd& n);

    double getMaxBorderFuncValuesAtN(int theta_idx,
                                     const std::vector<int>& theta_end_ids,
                                     int max_switches,
                                     int phase,
                                     const Eigen::VectorXd& n);

    std::vector<double> getBorderConditions(int switch_phase,
                                           int theta_idx,
                                           double theta,
                                           int switch_cnt);

    std::tuple<std::unique_ptr<Highs>,
               std::vector<double>,
               std::vector<double>>
    initializeHighs(int phase);

    void updateHighsRhsUpperBounds(int phase,
                                   const std::vector<double>& x_prev_vec);

    std::vector<double> solveLp(int phase);

    void precomputeMatrices();

    void run();
};

}  // namespace piecewise_affine_approximator

#endif  // HCPWA_PIECEWISE_AFFINE_APPROXIMATOR_H
