#pragma once

#include <matrix.hpp>
// expected row format: a*x + b*y < c
matrix<double> GetHullPoints(const matrix<double>& inequalities);