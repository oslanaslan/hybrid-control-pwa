#include <cdd.hpp>
#include <lineareq.hpp>
#include <string>
#include <stdexcept>
#include <utility.hpp>

matrix<double> GetHullPoints(const matrix<double>& inequalities) {
  dd_MatrixPtr m = dd_CreateMatrix(inequalities.rows(), inequalities.cols());
  m->representation = dd_Inequality;
  m->numbtype = dd_Real;

  // Convert inequalities to CDD format (a*x + b*y + c*z - d <= 0 becomes -a*x -
  // b*y - c*z + d >= 0)
  // Format b|A solving b+Ax>0
  for (size_t row = 0; row < inequalities.rows(); row++) {
    for (size_t col = 0; col < inequalities.cols() - 1; col++) {
      const double value = inequalities[row, col];
      dd_set_d(m->matrix[row][col + 1], -value);
    }
    const double value = inequalities[row].back();
    dd_set_d(m->matrix[row][0], value);
  }
  dd_ErrorType err;
  dd_PolyhedraPtr poly = dd_DDMatrix2Poly(m, &err);
  // defer _ = [&m, &poly] {
  //   dd_FreeMatrix(m);
  //   dd_FreePolyhedra(poly);
  // };
  if (err != dd_NoError) {
    throw std::runtime_error("error " + std::to_string((int)err));
  }
  bool feasible = (poly->child->CompStatus == dd_AllFound);

  if (!feasible) {
    throw std::runtime_error("not feasible");
  }

  dd_MatrixPtr v = dd_CopyGenerators(poly);
  // defer _ = [&v] {
  //   dd_FreeMatrix(v);
  // };

  matrix<double> result(inequalities.cols() - 1, 0);
  for (int i = 0; i < v->rowsize; ++i) {
    if (dd_get_d(v->matrix[i][0]) == 1.0) {
      result.emplace_back();
      for (int j = 0; j < result.cols(); j++) {
        // result.push_back(std::span(v->matrix[i] + 1, result.cols())
        //                  | std::ranges::views::transform([](auto& x) {
        //                      return dd_get_d(x);
        //                    }));
        // const auto value = dd_get_d(v->matrix[i][j + 1]);
        // result[i, j] = value;
        const auto value = dd_get_d(v->matrix[i][j + 1]);
        result.back()[j] = value;
      }
    }
  }
  return result;
}