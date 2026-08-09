// Compile-and-link smoke test for the barycentric approximator.
//
// It deliberately touches neither the geometry nor HiGHS: the heavy end-to-end
// run lives in barycentric_affine_solver.cpp and is executed separately on a
// machine with enough resources. What this test protects is the build itself —
// a header that stopped compiling, a changed constructor signature, a missing
// ApproximationMode definition, an unresolved symbol at link time.

#include <gtest/gtest.h>

#include <stdexcept>

#include <barycentric_affine_approximator.hpp>

namespace {

barycentric_affine_approximator::SystemParams makeParams() {
  // Same numbers as barycentric_affine_solver.cpp, so the two tests stay
  // comparable if the heavy one is ever run next to this one.
  return barycentric_affine_approximator::SystemParams{
      /*N=*/160.0,   /*F=*/0.5,     /*v=*/0.017,   /*w=*/0.0055,
      /*b51=*/0.6,   /*b57=*/0.4,   /*b84=*/0.8,   /*b86=*/0.2,
      /*b31=*/0.6,   /*b36=*/0.4,   /*b24=*/0.7,   /*b27=*/0.3,
      /*f2min=*/0.44, /*f3min=*/0.23, /*f5min=*/0.09, /*f8min=*/0.23,
      /*f2max=*/0.46, /*f3max=*/0.27, /*f5max=*/0.11, /*f8max=*/0.27};
}

}  // namespace

TEST(barycentric_smoke, constructs_in_both_modes) {
  using barycentric_affine_approximator::ApproximationMode;
  using barycentric_affine_approximator::BarycentricAffineApproximator;

  // Two instances in one process: the logger is fetched with spdlog::get before
  // being created, so the second construction must not throw.
  BarycentricAffineApproximator upper(
      /*t_max=*/300.0, /*t_split_count=*/10, /*tau_min=*/60.0,
      /*tau_max=*/120.0, makeParams(), /*highs_verbose=*/false,
      ApproximationMode::Upper);
  BarycentricAffineApproximator lower(300.0, 10, 60.0, 120.0, makeParams(),
                                      false, ApproximationMode::Lower);

  EXPECT_EQ(upper.approximationMode(), ApproximationMode::Upper);
  EXPECT_EQ(lower.approximationMode(), ApproximationMode::Lower);

  // Cheap parameter accessors: no geometry, no LP, no file system.
  EXPECT_DOUBLE_EQ(upper.getBetaParamForAxis(3 - 1, 1 - 1), 0.6);
  EXPECT_DOUBLE_EQ(lower.getBetaParamForAxis(5 - 1, 7 - 1), 0.4);

  const auto [f_min, f_max] = upper.getFMinMaxForAxis(2 - 1);
  EXPECT_DOUBLE_EQ(f_min, 0.44);
  EXPECT_DOUBLE_EQ(f_max, 0.46);
}

TEST(barycentric_smoke, defaults_to_upper_mode) {
  using barycentric_affine_approximator::ApproximationMode;
  using barycentric_affine_approximator::BarycentricAffineApproximator;

  // The mode argument is last and defaulted, so the pre-existing six-argument
  // call sites keep compiling unchanged.
  BarycentricAffineApproximator approximator(300.0, 10, 60.0, 120.0,
                                             makeParams(), false);
  EXPECT_EQ(approximator.approximationMode(), ApproximationMode::Upper);

  approximator.setValidate(true);
  approximator.setValidate(false);
}

TEST(barycentric_smoke, rejects_invalid_construction) {
  using barycentric_affine_approximator::BarycentricAffineApproximator;

  EXPECT_THROW(BarycentricAffineApproximator(300.0, /*t_split_count=*/0, 60.0,
                                             120.0, makeParams()),
               std::invalid_argument);
  EXPECT_THROW(BarycentricAffineApproximator(300.0, 10, /*tau_min=*/120.0,
                                             /*tau_max=*/60.0, makeParams()),
               std::invalid_argument);
}
