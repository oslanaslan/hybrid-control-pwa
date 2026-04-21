#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>
#include <iostream>
#include <iomanip>

namespace interval_building {

struct ThetaTIndexLists {
  // original_t_by_k_theta[k][theta_idx] = список допустимых t_idx из раздела 5
  std::unordered_map<int, std::unordered_map<int, std::vector<int>>>
      original_t_by_k_theta;

  // expanded_t_by_k_theta[k][theta_idx] = исходный список + расширение до theta
  // (не включая саму theta, если она была недостижима из раздела 5)
  std::unordered_map<int, std::unordered_map<int, std::vector<int>>>
      expanded_t_by_k_theta;
};

namespace detail {

inline void validateGrid(const std::vector<double>& grid) {
  if (grid.empty()) {
    throw std::invalid_argument("Grid must not be empty.");
  }
  for (size_t i = 1; i < grid.size(); ++i) {
    if (grid[i] < grid[i - 1]) {
      throw std::invalid_argument(
          "Grid must be sorted in nondecreasing order.");
    }
  }
}

inline std::vector<int> collectIndicesInClosedInterval(
    const std::vector<double>& grid, double low, double high, double eps) {
  std::vector<int> result;
  if (high < low - eps) {
    return result;
  }

  for (int i = 0; i < static_cast<int>(grid.size()); ++i) {
    const double x = grid[i];
    if (x < low - eps) {
      continue;
    }
    if (x <= high + eps) {
      result.push_back(i);
    } else {
      break;  // grid is sorted
    }
  }
  return result;
}

inline std::vector<int> collectIndicesInOpenInterval(
    const std::vector<double>& grid, double low, double high, double eps) {
  // (low, high), обе границы открытые
  std::vector<int> result;
  if (high <= low + eps) {
    return result;
  }

  for (int i = 0; i < static_cast<int>(grid.size()); ++i) {
    const double x = grid[i];
    if (x <= low + eps) {
      continue;
    }
    if (x < high - eps) {
      result.push_back(i);
    } else {
      break;  // grid is sorted
    }
  }
  return result;
}

inline void sortAndUnique(std::vector<int>& ids) {
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
}

}  // namespace detail

// Строит структуры только для k = 1..K.
// Это соответствует слоям, где theta действительно играет роль момента
// переключения. Для k = 0 в разделе 5 theta фиктивный, поэтому здесь он
// сознательно не включён.
inline ThetaTIndexLists buildThetaToTIndexLists(
    double T, double tau_min, double tau_max,
    const std::vector<double>& time_grid, int K, double eps = 1e-12) {
  using detail::collectIndicesInClosedInterval;
  using detail::collectIndicesInOpenInterval;
  using detail::sortAndUnique;
  using detail::validateGrid;

  if (T < 0.0) {
    throw std::invalid_argument("T must be nonnegative.");
  }
  if (tau_min <= 0.0 || tau_max <= 0.0) {
    throw std::invalid_argument("tau_min and tau_max must be positive.");
  }
  if (tau_min > tau_max) {
    throw std::invalid_argument("Expected tau_min <= tau_max.");
  }

  validateGrid(time_grid);

  // Желательно, чтобы сетка была внутри [0, T].
  if (time_grid.front() < -eps || time_grid.back() > T + eps) {
    throw std::invalid_argument("time_grid must lie inside [0, T].");
  }

  ThetaTIndexLists out;

  // =========================
  // Шаг 1. Цикл по слоям k
  // =========================
  for (int k = 1; k <= K; ++k) {
    const double t_min_k = std::max(0.0, T - static_cast<double>(k) * tau_max);
    const double t_max_k = T - static_cast<double>(k) * tau_min;

    // Если слой пустой, просто создаём пустые контейнеры и идём дальше.
    if (t_max_k < t_min_k - eps) {
      out.original_t_by_k_theta[k] = {};
      out.expanded_t_by_k_theta[k] = {};
      continue;
    }

    // ============================================================
    // Шаг 2. Найти все индексы t, попавшие в [t_min^(k), t_max^(k)]
    // ============================================================
    const std::vector<int> t_ids_in_layer
        = collectIndicesInClosedInterval(time_grid, t_min_k, t_max_k, eps);

    // Временная структура "для фиксированного t -> список допустимых theta"
    std::unordered_map<int, std::vector<int>> theta_ids_by_t_idx;
    theta_ids_by_t_idx.reserve(t_ids_in_layer.size());

    // =======================================================================
    // Шаг 3. Для каждого допустимого t построить отрезок theta из раздела 5:
    //         theta in [t, min(t + tau_max, T)]
    // и собрать индексы theta на той же сетке
    // =======================================================================
    for (int t_idx : t_ids_in_layer) {
      const double t_value = time_grid[t_idx];
      const double theta_low = t_value;
      const double theta_high = std::min(t_value + tau_max, T);

      theta_ids_by_t_idx[t_idx] = collectIndicesInClosedInterval(
          time_grid, theta_low, theta_high, eps);
    }

    // =========================================================================
    // Шаг 4. Переупорядочить "для фиксированного t -> список theta"
    // в "для фиксированного theta -> список t"
    // =========================================================================
    std::unordered_map<int, std::vector<int>> original_theta_to_t;
    for (const auto& [t_idx, theta_ids] : theta_ids_by_t_idx) {
      for (int theta_idx : theta_ids) {
        original_theta_to_t[theta_idx].push_back(t_idx);
      }
    }

    // Упорядочить и убрать дубликаты, чтобы списки были аккуратными.
    for (auto& [theta_idx, t_ids] : original_theta_to_t) {
      sortAndUnique(t_ids);
    }

    out.original_t_by_k_theta[k] = original_theta_to_t;

    // ====================================================================================
    // Шаг 5. Построить расширенный список:
    // для фиксированного theta добавить индексы t, лежащие между
    // последним допустимым t и theta, не включая саму theta.
    //
    // То есть:
    //   expanded = original U { t_idx : t in (t_right_original, theta) }
    //
    // где
    //   t_right_original = min(t_max^(k), theta)
    // ====================================================================================
    std::unordered_map<int, std::vector<int>> expanded_theta_to_t;
    expanded_theta_to_t.reserve(original_theta_to_t.size());

    for (const auto& [theta_idx, original_t_ids] : original_theta_to_t) {
      const double theta_value = time_grid[theta_idx];

      // Скопировать исходный список
      std::vector<int> expanded_t_ids = original_t_ids;

      // Исходный правый конец непрерывного допустимого сегмента по t:
      const double t_right_original = std::min(t_max_k, theta_value);

      // Добавляем только "хвост" между t_right_original и theta, без theta.
      // const std::vector<int> extra_ids = collectIndicesInOpenInterval(
          // time_grid, t_right_original, theta_value, eps);
      const std::vector<int> extra_ids = collectIndicesInClosedInterval(
          time_grid, t_right_original, theta_value, eps);
      expanded_t_ids.insert(expanded_t_ids.end(), extra_ids.begin(),
                            extra_ids.end());
      sortAndUnique(expanded_t_ids);

      expanded_theta_to_t[theta_idx] = std::move(expanded_t_ids);
    }

    out.expanded_t_by_k_theta[k] = std::move(expanded_theta_to_t);
  }

  return out;
}

inline void prettyPrintThetaTLists(const ThetaTIndexLists& data,
                                   const std::vector<double>& grid,
                                   int precision = 3) {
  std::cout << std::fixed << std::setprecision(precision);

  for (const auto& [k, theta_map] : data.original_t_by_k_theta) {
    std::cout << "=== k = " << k << " ===\n";

    const auto& expanded_map = data.expanded_t_by_k_theta.at(k);

    for (const auto& [theta_idx, t_ids] : theta_map) {
      double theta_val = grid[theta_idx];

      std::cout << "theta[" << theta_idx << "] = " << theta_val << ":\n";

      // --- original ---
      std::cout << "    t     = [";
      for (size_t i = 0; i < t_ids.size(); ++i) {
        std::cout << grid[t_ids[i]];
        if (i + 1 < t_ids.size()) {
          std::cout << ", ";
        }
      }
      std::cout << "]\n";

      // --- expanded ---
      const auto& t_ext_ids = expanded_map.at(theta_idx);

      std::cout << "    t_ext = [";
      for (size_t i = 0; i < t_ext_ids.size(); ++i) {
        std::cout << grid[t_ext_ids[i]];
        if (i + 1 < t_ext_ids.size()) {
          std::cout << ", ";
        }
      }
      std::cout << "]\n\n";
    }

    std::cout << "\n";
  }
}
}  // namespace interval_building