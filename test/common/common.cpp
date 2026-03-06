#include <gtest/gtest.h>
#include <hcpwa.hpp>
#include <uniqie_pool.hpp>
#include <symbolic.hpp>
#include "morph.hpp"
#include "types.hpp"

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

TEST(common, symbolic) {
    // NOLINTNEXTLINE
    using namespace hcpwa::symbols;
    {
        hcpwa::Line<1> a = hcpwa::kZeroVec;
        hcpwa::Line<2> b = hcpwa::kZeroVec;
        hcpwa::Line<3> c = hcpwa::kZeroVec;
    }
    {
        auto v = Linearize(X<0>{});
        hcpwa::Line<1> l1 = {1, 0};
        ASSERT_EQ(v, l1);
    }
    {
        auto v = Linearize(X<0>{} + X<1>{});
        hcpwa::Line<2> ex = {1, 1, 0};
        ASSERT_EQ(v, ex);
    }
    {
        auto v = Linearize(X<0>{} + X<1>{} - C{2});
        hcpwa::Line<2> ex = {1, 1, -2};
        ASSERT_EQ(v, ex);
    }
    {
        auto v = SymMin(X<0>{} + C{2}, C{1} - X<0>{});
        auto lines = MinOptions(v);
        hcpwa::LineSet<1> ex = {{1, 2}, {-1, 1}};
        ASSERT_EQ(lines, ex);
    }
    {
        auto v = SymMin(X<0>{} + C{2}, C{1} - X<0>{});
        auto resolutions = MinResolutions(v);
        constexpr int kDim
            = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
        auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);

        ASSERT_EQ(masks.size(), 2);
        ASSERT_EQ(masks[0], 0b01);
        ASSERT_EQ(masks[1], 0b10);
    }
    {
        auto v = SymMin(C{3} * X<0>{}, C{6} - C{3} * X<0>{});
        auto resolutions = MinResolutions(v);
        constexpr int kDim
            = hcpwa::VectorSize<decltype(resolutions.front().second)>() - 1;
        auto [lines, masks] = hcpwa::ResolutionsToMasks<kDim>(resolutions);

        ASSERT_EQ(masks.size(), 2);
        ASSERT_EQ(masks[0], 0b01);
        ASSERT_EQ(masks[1], 0b10);
    }
}