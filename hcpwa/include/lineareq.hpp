#pragma once

#include <matrix.hpp>

namespace cddwrap {
/**
 * @brief Get the Hull Points object
 * 
 * @param inequalities in format 1* + b*y < c
 * @return matrix<double> 
 */
matrix<double> GetHullPoints(const matrix<double>& inequalities);
}  // namespace cddwrap