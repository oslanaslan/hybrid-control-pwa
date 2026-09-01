#ifndef HCPWA_GLOBAL_AFFINE_APPROXIMATOR_H
#define HCPWA_GLOBAL_AFFINE_APPROXIMATOR_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Highs.h>
#include <cstddef>
#include <fstream>
#include <memory>
#include <mutex>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <spdlog/spdlog.h>

#include "util/assert_utils.hpp"
#include "util/value_function_utils.hpp"
#include "util/global_block_reduction_lp.hpp"
#include "interval_building.hpp"

namespace global_affine_approximator {

// Constants
extern const std::vector<int> kInIds;
extern const std::vector<int> kOutIds;
extern const int kPhases;
constexpr int kSpaceDim = 8;
constexpr int kVDeltaDim = kSpaceDim + 1;            // [V, v]
// The LP's column count now lives with the assembler that defines the layout:
// hcpwa::util::global_block_lp::kNumCols, which is [V(8) | v | s(8) | mu(3)].
// The old kLPCols = 17 is deliberately not kept here; it would contradict it.
constexpr double kEps = 1e-5;
// Tiny positive cost on s to pin s_i = |V_i| in the L1-residual objective
// (s no longer appears in the residual objective; keep << |c_vec| entries).
constexpr double kSPinWeight = 1e-6;

// Number of coordinate blocks per phase: A, B, C.
constexpr int kBlockCount = hcpwa::util::global_block_lp::kBlockCount;

// Tolerance for that check. It sits above the HiGHS feasibility tolerance so it
// does not fire on solutions the solver legitimately calls optimal.
constexpr double kResidualValidationTol = 1e-4;

// Geometry of one coordinate block of one phase.
//
// An 8D area is the Cartesian product of three low-dimensional block polytopes
// (notes/step7_block_reduction.md, Lemma 1), so carrying the three block vertex
// sets replaces materialising their product. The product has 108-192 vertices
// per area across roughly 1.36M areas; the blocks have about 509 regions in
// total. The 8D vertex lists this replaces were also truncated to 12 vertices
// per area, which is why the old LP was not a bound at all: see
// docs/barycentric_block_reduction_context.md part I.
//
// There are no barycentric charts on this path, so unlike the barycentric
// BlockGeometry there is no layer or local_axis table here.
struct GlobalBlockGeometry {
    std::array<int, 3> coords{};  // state coordinates, ascending
    int coord_count = 0;          // 3 for A and B, 2 for C
    int num_regions = 0;
    // Per block region, all vertices, each of length coord_count.
    std::vector<std::vector<Eigen::VectorXd>> vertices;

    // Per block region, restricted to the block's own coordinates.
    std::vector<Eigen::MatrixXd> a_matr;
    std::vector<Eigen::VectorXd> f_vec;
    std::vector<Eigen::MatrixXd> q_c_matr;
    std::vector<Eigen::VectorXd> q_c_vec;
    std::vector<Eigen::MatrixXd> q_r_matr;
    std::vector<Eigen::VectorXd> q_r_vec;
    std::vector<Eigen::VectorXd> g_vec;
    std::vector<double> g_scal;
};

using ValueFunction = hcpwa::util::ValueFunction;

enum class ApproximationMode { Upper, Lower };

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

class GlobalAffineApproximator {
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

    bool highs_verbose_ = false;
    ApproximationMode approximation_mode_ = ApproximationMode::Upper;
    // Off by default: the bound check runs on every solved step. It is cheap
    // per triple but the sample is large, so it is opt-in.
    bool validate_ = false;

    ValueFunction value_function_;
    std::vector<Eigen::VectorXd> cube_angle_vertices_;
    // The 8D area vertex lists and the per-area CTM matrices that used to live
    // here are gone. Both were indexed by full 8D area, of which there are
    // M_A*M_B*M_C; everything is now per block region, of which there are
    // M_A+M_B+M_C.
    std::array<std::array<GlobalBlockGeometry, kBlockCount>, 2> blocks_;
    std::array<int, 2> num_regions_{};

    // The assembled reduced LP per phase. Owns the constraint matrix, the
    // objective and the per-row data needed to refresh right-hand sides.
    std::array<hcpwa::util::global_block_lp::GlobalReducedLp, 2> reduced_lps_;

    std::shared_ptr<spdlog::logger> logger_;
    std::vector<std::unique_ptr<Highs>> highs_solvers_;
    std::vector<std::vector<double>> row_lowers_;
      std::vector<std::unique_ptr<std::mutex>> solver_mutexes_;

    interval_building::ThetaTIndexLists theta_t_index_lists_;

   public:
    GlobalAffineApproximator(double t_max, int t_split_count,
                             double tau_min, double tau_max,
                             const SystemParams& system_params,
                             bool highs_verbose = false,
                             ApproximationMode mode = ApproximationMode::Upper);

    // Enables the per-step bound check of validateStepResiduals().
    void setValidate(bool validate) { validate_ = validate; }

    double getBetaParamForAxis(int i, int j) const;

    std::pair<double, double> getFMinMaxForAxis(int i) const;

    void getIntersectionPoints();

    std::pair<Eigen::RowVectorXd, Eigen::RowVectorXd> getFIJMinResolution(
        int i, int j, const Eigen::VectorXd& n) const;

    // Centroid of one block region, from its own COMPLETE vertex set. The old
    // whole-area centroid averaged a truncated vertex set whose hull was
    // segment x segment x polygon, a 4-dimensional slice of an 8-dimensional
    // area; midpoints of two vertices sharing a facet land on the area
    // boundary, where the min is tied, so branch resolution was being decided
    // inside the kEps tie window or throwing outright.
    Eigen::VectorXd blockCentroidCoords(int phase, int block_id,
                                        int j_block) const;

    // Expands a block point into an 8-vector, padding coordinates outside the
    // block with NaN. See the definition for why NaN and not zero.
    Eigen::VectorXd blockPointToFullState(int phase, int block_id,
                                          const Eigen::VectorXd& nu) const;

    std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::VectorXd, double>
    getBlockAMatrFVecGVecAndGScalJ(int phase, int block_id, int j_block) const;

    std::tuple<Eigen::MatrixXd, Eigen::VectorXd, Eigen::MatrixXd,
               Eigen::VectorXd>
    getBlockQQ(int phase, int block_id, int j_block) const;

    void precomputeSystemMatrices(int phase);

    hcpwa::util::global_block_lp::GlobalReducedLp prepareLpMatrices(int phase);

    // Checks the property the whole construction exists to provide: that
    // s * F(a,b,c) <= 0 on the ORIGINAL row set, sampled over triples from
    // R_A x R_B x R_C. O(1) per triple, no LP and no 8D vertex. This is the
    // check that would have caught the vertex truncation immediately.
    void validateStepResiduals(int phase,
                               const std::vector<double>& v_prev,
                               const std::vector<double>& v_cur) const;

    double getBorderFuncValuesAtN(int r, int theta_idx, int theta_end_idx,
                                  int phase, const Eigen::VectorXd& n);

    double getMaxBorderFuncValuesAtN(int theta_idx,
                                     const std::vector<int>& theta_end_ids,
                                     int max_switches, int phase,
                                     const Eigen::VectorXd& n);

    double getMinBorderFuncValuesAtN(int theta_idx,
                                     const std::vector<int>& theta_end_ids,
                                     int max_switches, int phase,
                                     const Eigen::VectorXd& n);

    std::vector<double> getBorderConditions(int switch_phase, int theta_idx,
                                            double theta, int switch_cnt);

    std::tuple<std::unique_ptr<Highs>, std::vector<double>, std::vector<double>>
    initializeHighs(int phase);

    void updateHighsRhsUpperBounds(int phase, int solver_index,
                                   const std::vector<double>& x_prev_vec);

    std::vector<double> solveLp(int solver_index);

    void precomputeMatrices();

    void run(const std::string& output_folder_path, int n_threads = 2);

    void dumpInitParamsToJson(const std::string& filepath) const;

    // buildLPSegmentForJ was deleted with the per-area LP assembly; the
    // reduced LP is built by hcpwa::util::global_block_lp::assembleGlobalReducedLp.

};

}  // namespace global_affine_approximator

#endif  // HCPWA_GLOBAL_AFFINE_APPROXIMATOR_H
