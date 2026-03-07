#include <gtest/gtest.h>
#include <hcpwa.hpp>
#include <algo.hpp>
#include "cddwrap/cdd.hpp"
#include <uniqie_pool.hpp>
#include <symbolic.hpp>
#include <Eigen/Dense>
#include <global_affine_approximator.h>
#include "utility.hpp"

#include "morph.hpp"
#include "types.hpp"

static global_affine_approximator::GlobalAffineApproximator
create_linear_approximator() {
    cddwrap::global_init();
    defer _ = &cddwrap::global_free;
    // Same parameters as user_algo_tests (compute_areas_vertices)
    constexpr double N = 100;
    constexpr double F = 15;
    constexpr double v = 0.2;
    constexpr double w = 0.5;
    constexpr double b51 = 0.5;
    constexpr double b57 = 0.5;
    constexpr double b84 = 0.5;
    constexpr double b86 = 0.5;
    constexpr double b31 = 0.5;
    constexpr double b36 = 0.5;
    constexpr double b24 = 0.5;
    constexpr double b27 = 0.5;
    constexpr double f2min = 5;
    constexpr double f3min = 5;
    constexpr double f5min = 5;
    constexpr double f8min = 5;
    constexpr double f2max = 10;
    constexpr double f3max = 10;
    constexpr double f5max = 10;
    constexpr double f8max = 10;

    global_affine_approximator::SystemParams system_params{
        N,   F,   v,     w,     b51,   b57,   b84,   b86,   b31,   b36,
        b24, b27, f2min, f3min, f5min, f8min, f2max, f3max, f5max, f8max};

    // LinearApproximator constructor parameters
    constexpr double t_max = 10.0;
    constexpr int t_split_count = 100;
    constexpr int max_switches = 2;
    constexpr double tau_min = 0.0;
    constexpr double tau_max = 1.0;

    return global_affine_approximator::GlobalAffineApproximator(
        t_max, t_split_count, max_switches, tau_min, tau_max, system_params);
}

TEST(common, affine_solver) {
    auto linear_approximator = create_linear_approximator();

    linear_approximator.run();
}