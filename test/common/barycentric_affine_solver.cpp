#include <gtest/gtest.h>
#include <hcpwa.hpp>
#include <algo.hpp>
#include "cddwrap/cdd.hpp"
#include <uniqie_pool.hpp>
#include <symbolic.hpp>
#include <Eigen/Dense>
#include <barycentric_affine_approximator.hpp>
#include "utility.hpp"

#include "morph.hpp"
#include "types.hpp"

static barycentric_affine_approximator::BarycentricAffineApproximator
create_barycentric_approximator() {
  cddwrap::global_init();
  defer _ = &cddwrap::global_free;
  // Same parameters as user_algo_tests (compute_areas_vertices)
  constexpr double N = 160;
  constexpr double F = 0.5;
  constexpr double v = 0.017;
  constexpr double w = 0.0055;
  constexpr double b51 = 0.6;
  constexpr double b57 = 0.4;
  constexpr double b84 = 0.8;
  constexpr double b86 = 0.2;
  constexpr double b31 = 0.6;
  constexpr double b36 = 0.4;
  constexpr double b24 = 0.7;
  constexpr double b27 = 0.3;
  // constexpr double f2min = 0.44;
  // constexpr double f3min = 0.23;
  // constexpr double f5min = 0.09;
  // constexpr double f8min = 0.23;
  // constexpr double f2max = 0.44;
  // constexpr double f3max = 0.23;
  // constexpr double f5max = 0.09;
  // constexpr double f8max = 0.23;
  constexpr double f2min = 0.44;
  constexpr double f3min = 0.23;
  constexpr double f5min = 0.09;
  constexpr double f8min = 0.23;
  constexpr double f2max = 0.46;
  constexpr double f3max = 0.27;
  constexpr double f5max = 0.11;
  constexpr double f8max = 0.27;

  barycentric_affine_approximator::SystemParams system_params{
      N,   F,   v,     w,     b51,   b57,   b84,   b86,   b31,   b36,
      b24, b27, f2min, f3min, f5min, f8min, f2max, f3max, f5max, f8max};

  // Same timing / mode params as affine_solver.cpp
  constexpr double t_max = 1200.0;
  // constexpr double tau_min = 0.06;
  // constexpr double tau_max = 0.12;
  constexpr double tau_min = 10;
  constexpr double tau_max = 50;
  // constexpr double t_max = 1.0;
  constexpr int t_split_count = 240;
  // constexpr int max_switches = 5;

  barycentric_affine_approximator::ApproximationMode approximation_mode
      = barycentric_affine_approximator::ApproximationMode::Lower;
  bool highs_verbose = true;
  return barycentric_affine_approximator::BarycentricAffineApproximator(
      t_max, t_split_count, tau_min, tau_max, system_params, highs_verbose,
      approximation_mode);
}

TEST(common, barycentric_affine_solver) {
  auto barycentric_approximator = create_barycentric_approximator();

  barycentric_approximator.run("/root/gitlab/hybrid-control-pwa/results/barycentric/lower/first_run/", 16);
}
