#include <gtest/gtest.h>
#include <hcpwa.hpp>
#include <uniqie_pool.hpp>

TEST(common, pool) {
  using Comp = hcpwa::FixedPrecisionComparator<hcpwa::Float, 0.1f>;
  hcpwa::UniquePool<hcpwa::Float> pool(Comp{});
  hcpwa::Float a = 0.1;
  hcpwa::Float b = 0.2;
  hcpwa::Float c = 0.3;
  hcpwa::Float d = 0.3001;
  ASSERT_EQ(pool.Unique(a), a);
  ASSERT_EQ(pool.Unique(b), b);
  ASSERT_EQ(pool.Unique(c), c);
  ASSERT_EQ(pool.Unique(d), c);
  // check pool state
  ASSERT_EQ(pool.Unique(a), a);
  ASSERT_EQ(pool.Unique(b), b);
  ASSERT_EQ(pool.Unique(c), c);
  ASSERT_EQ(pool.Unique(d), c);
}
