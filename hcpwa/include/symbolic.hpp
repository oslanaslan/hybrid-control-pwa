#pragma once

#include <morph.hpp>

namespace hcpwa::symbols {

namespace internal {
consteval auto Max(auto a, auto b) {
  return a < b ? b : a;
}
template <typename T, typename... Ts>
  requires(sizeof...(Ts) > 1)
consteval auto Max(T first, Ts... args) {
  return Max(first, Max(args...));
}
}  // namespace internal

template <typename T>
concept SymExpr = requires(T e) {
  { Linearize(e) } -> Vector;
};

template <typename T>
concept SymMinExpr = requires(T e) {
  {
    e.apply([](SymExpr auto ex) {
      return ex;
    })
  };
  { e.Dim() } -> std::same_as<int>;
};

template <typename T>
  requires requires(T v) {
    { v.Linearize() } -> Vector;
  }
constexpr decltype(auto) Linearize(T v) {
  return v.Linearize();
}

constexpr auto MakeSymExpr(auto lambda) {
  struct Expr {
    decltype(lambda) f;
    constexpr auto Linearize() {
      return f();
    }
  };
  return Expr{.f = std::move(lambda)};
}

template <int N>
  requires(N >= 0)
struct X {
  consteval friend int Dim(X) {
    return N + 1;
  }

  constexpr auto Linearize() {
    Line<N + 1> result = kZeroVec;
    result[N] = 1;
    return result;
  }
};

struct C {
  Float value;
  constexpr auto Linearize() {
    return Line<1>{0, value};
  }
};

constexpr auto operator+(SymExpr auto a, SymExpr auto b) {
  return MakeSymExpr([a, b] {
    auto la = Linearize(a);
    auto lb = Linearize(b);
    constexpr int kA = VectorSize<decltype(la)>() - 1;
    constexpr int kB = VectorSize<decltype(lb)>() - 1;
    constexpr int kD = internal::Max(kA, kB);
    return Expand<kD>(la) + Expand<kD>(lb);
  });
}

constexpr auto operator-(SymExpr auto a, SymExpr auto b) {
  return MakeSymExpr([a, b] {
    auto la = Linearize(a);
    auto lb = Linearize(b);
    constexpr int kA = VectorSize<decltype(la)>() - 1;
    constexpr int kB = VectorSize<decltype(lb)>() - 1;
    constexpr int kD = internal::Max(kA, kB);
    return Expand<kD>(la) - Expand<kD>(lb);
  });
}

constexpr auto operator*(SymExpr auto a, C b) {
  return MakeSymExpr([a, b] {
    auto la = Linearize(a);
    la *= b.value;
    return la;
  });
}

constexpr auto operator*(C b, SymExpr auto a) {
  return MakeSymExpr([a, b] {
    auto la = Linearize(a);
    la *= b.value;
    return la;
  });
}

template <SymExpr... Args>
constexpr auto SymMin(Args... args) {
  constexpr int kMaxD
      = internal::Max(VectorSize<decltype(Linearize(args))>()...);

  struct Resize {
    constexpr auto Linearize() {
      return (Vec<kMaxD>)kZeroVec;
    }
  };

  auto l = [args...](auto f) {
    return SymMin(f(args + Resize{})...);
  };
  struct Min {
    decltype(l) apply;
    consteval int Dim() {
      return kMaxD - 1;
    }
  };
  return Min{.apply = std::move(l)};
}

template <SymMinExpr... Args>
constexpr auto SymMin(Args... mins) {
  constexpr int kMaxD = internal::Max(mins.Dim()...);

  auto l = [mins...](auto f) {
    return SymMin(mins.apply(f)...);
  };
  struct Min {
    decltype(l) apply;
    consteval int Dim() {
      return kMaxD;
    }
  };
  return Min{.apply = std::move(l)};
}

constexpr auto operator+(SymMinExpr auto a, SymExpr auto b) {
  return a.apply([b](SymExpr auto e) {
    return e + b;
  });
}
constexpr auto operator-(SymMinExpr auto a, SymExpr auto b) {
  return a.apply([b](SymExpr auto e) {
    return e - b;
  });
}
constexpr auto operator+(SymExpr auto b, SymMinExpr auto a) {
  return a.apply([b](SymExpr auto e) {
    return b + e;
  });
}
constexpr auto operator-(SymExpr auto b, SymMinExpr auto a) {
  return a.apply([b](SymExpr auto e) {
    return b - e;
  });
}

constexpr auto MinOptions(SymMinExpr auto expr) {
  constexpr int kDim = expr.Dim();
  LineSet<kDim> result;
  expr.apply([&result](SymExpr auto ex) {
    result.push_back(Linearize(ex));
    return ex;
  });
  return result;
}

}  // namespace hcpwa::symbols