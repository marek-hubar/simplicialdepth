#include <Rcpp.h>
#include <limits>

#include "helpers.h"

using namespace Rcpp;

namespace {

struct RelativeOrderEntry {
    double slope;     // -x/y in the local (ray-aligned) coordinate system
    int is_upper;         // 1 if relative angle in [0, PI), else 0
};

struct GlobalOrderEntry {
    double slope;     // -x/y in the global coordinate system
    int is_upper;       // 1 if atan2(point) in [0, PI), else 0
    int original_index;
};

int upper_half_label(double x, double y) {
    if (y > EPS) return 1;
    if (y < -EPS) return 0;
    return (x >= 0.0) ? 1 : 0;
}

double slope_neg_x_over_y(double x, double y) {
    if (std::fabs(y) <= EPS) {
        return (x >= 0.0) ? -std::numeric_limits<double>::infinity()
                          : std::numeric_limits<double>::infinity();
    }
    return -x / y;
}

std::vector<Point2D> matrix_to_points(NumericMatrix points_mat) {
    const int n = points_mat.nrow();
    std::vector<Point2D> points(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        points[static_cast<std::size_t>(i)] = {points_mat(i, 0), points_mat(i, 1)};
    }
    return points;
}

long long circular_asd_single_query_01(Point2D ray, const std::vector<Point2D>& points) {
    const int m = static_cast<int>(points.size());
    if (m < 2) return 0;
    int lower_count = 0;

    std::vector<RelativeOrderEntry> ordered_points(static_cast<std::size_t>(m));
    for (int i = 0; i < m; ++i) {
        const Point2D& p = points[static_cast<std::size_t>(i)];
        const double x_local = dot(ray, p);
        const double y_local = ray.x * p.y - ray.y * p.x;
        const int is_upper = upper_half_label(x_local, y_local);
        const double slope = slope_neg_x_over_y(x_local, y_local);
        ordered_points[static_cast<std::size_t>(i)] = {slope, is_upper};
        if (is_upper == 0) ++lower_count;
    }

    parallel_sort(ordered_points, [](const RelativeOrderEntry& a, const RelativeOrderEntry& b) {
        if (a.slope < b.slope) return true;
        if (a.slope > b.slope) return false;
        return a.is_upper > b.is_upper;
    });

    // In 01-sequence terms, this counts the number of subsequences "10".
    long long depth = 0;
    for (const auto& entry : ordered_points) {
        if (entry.is_upper == 1) {
            depth += lower_count;
        } else {
            --lower_count;
        }
    }
    return depth;
}

std::vector<long long> circular_asd_all_points_01(const std::vector<Point2D>& points) {
    const int n = static_cast<int>(points.size());
    std::vector<long long> depths(static_cast<std::size_t>(n), 0);
    if (n <= 2) return depths;

    // Sort once by -x/y (equivalent to angle modulo PI), then reuse this shared order for
    // all candidate rays.
    std::vector<GlobalOrderEntry> entries(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        const Point2D& p = points[static_cast<std::size_t>(i)];
        entries[static_cast<std::size_t>(i)] = {
            slope_neg_x_over_y(p.x, p.y),
            upper_half_label(p.x, p.y),
            i
        };
    }

    parallel_sort(entries, [](const GlobalOrderEntry& a, const GlobalOrderEntry& b) {
        if (a.slope < b.slope) return true;
        if (a.slope > b.slope) return false;
        if (a.is_upper != b.is_upper) return a.is_upper > b.is_upper;  // deterministic ties
        return a.original_index < b.original_index;
    });

    std::vector<int> upper_flags(static_cast<std::size_t>(n));
    std::vector<int> sorted_to_original(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        upper_flags[static_cast<std::size_t>(i)] = entries[static_cast<std::size_t>(i)].is_upper;
        sorted_to_original[static_cast<std::size_t>(i)] = entries[static_cast<std::size_t>(i)].original_index;
    }

    std::vector<long long> ones_prefix(static_cast<std::size_t>(n + 1), 0);
    std::vector<long long> zeros_prefix(static_cast<std::size_t>(n + 1), 0);
    std::vector<long long> inv10_prefix(static_cast<std::size_t>(n + 1), 0);
    std::vector<long long> inv01_prefix(static_cast<std::size_t>(n + 1), 0);

    for (int i = 0; i < n; ++i) {
        ones_prefix[static_cast<std::size_t>(i + 1)] = ones_prefix[static_cast<std::size_t>(i)] + upper_flags[static_cast<std::size_t>(i)];
        zeros_prefix[static_cast<std::size_t>(i + 1)] = zeros_prefix[static_cast<std::size_t>(i)] + (1 - upper_flags[static_cast<std::size_t>(i)]);
        inv10_prefix[static_cast<std::size_t>(i + 1)] =
            inv10_prefix[static_cast<std::size_t>(i)] +
            (upper_flags[static_cast<std::size_t>(i)] == 0 ? ones_prefix[static_cast<std::size_t>(i)] : 0);
        inv01_prefix[static_cast<std::size_t>(i + 1)] =
            inv01_prefix[static_cast<std::size_t>(i)] +
            (upper_flags[static_cast<std::size_t>(i)] == 1 ? zeros_prefix[static_cast<std::size_t>(i)] : 0);
    }

    std::vector<long long> ones_suffix(static_cast<std::size_t>(n + 1), 0);
    std::vector<long long> zeros_suffix(static_cast<std::size_t>(n + 1), 0);
    std::vector<long long> inv10_suffix(static_cast<std::size_t>(n + 1), 0);
    std::vector<long long> inv01_suffix(static_cast<std::size_t>(n + 1), 0);

    for (int i = n - 1; i >= 0; --i) {
        ones_suffix[static_cast<std::size_t>(i)] = ones_suffix[static_cast<std::size_t>(i + 1)] + upper_flags[static_cast<std::size_t>(i)];
        zeros_suffix[static_cast<std::size_t>(i)] = zeros_suffix[static_cast<std::size_t>(i + 1)] + (1 - upper_flags[static_cast<std::size_t>(i)]);
        inv10_suffix[static_cast<std::size_t>(i)] =
            inv10_suffix[static_cast<std::size_t>(i + 1)] +
            (upper_flags[static_cast<std::size_t>(i)] == 1 ? zeros_suffix[static_cast<std::size_t>(i + 1)] : 0);
        inv01_suffix[static_cast<std::size_t>(i)] =
            inv01_suffix[static_cast<std::size_t>(i + 1)] +
            (upper_flags[static_cast<std::size_t>(i)] == 0 ? ones_suffix[static_cast<std::size_t>(i + 1)] : 0);
    }

    // For each candidate ray k in sorted order, the remaining points appear as:
    //   suffix (k+1..n-1), then prefix (0..k-1).
    // Using prefix/suffix inversion counts, we can compute each depth in O(1), giving O(n log n)
    // overall (dominated by the initial sort).
    for (int k = 0; k < n; ++k) {
        const int ray_upper = upper_flags[static_cast<std::size_t>(k)];

        const long long ones_pref = ones_prefix[static_cast<std::size_t>(k)];
        const long long zeros_pref = zeros_prefix[static_cast<std::size_t>(k)];
        const long long ones_suf = ones_suffix[static_cast<std::size_t>(k + 1)];
        const long long zeros_suf = zeros_suffix[static_cast<std::size_t>(k + 1)];

        long long depth = 0;
        if (ray_upper == 1) {
            depth = inv10_suffix[static_cast<std::size_t>(k + 1)] + inv01_prefix[static_cast<std::size_t>(k)] + ones_suf * ones_pref;
        } else {
            depth = inv01_suffix[static_cast<std::size_t>(k + 1)] + inv10_prefix[static_cast<std::size_t>(k)] + zeros_suf * zeros_pref;
        }

        const int original_index = sorted_to_original[static_cast<std::size_t>(k)];
        depths[static_cast<std::size_t>(original_index)] = depth;
    }

    return depths;
}

}  // namespace

// [[Rcpp::export]]
long long circular_asd(NumericMatrix X, NumericVector x) {
    const std::vector<Point2D> points = matrix_to_points(X);
    const Point2D query = {x[0], x[1]};
    return circular_asd_single_query_01(query, points);
}

// [[Rcpp::export]]
DataFrame circular_asd_all_points(NumericMatrix points_mat) {
    const std::vector<Point2D> points = matrix_to_points(points_mat);
    const std::vector<long long> depth_counts = circular_asd_all_points_01(points);

    IntegerVector depth(depth_counts.size());
    for (std::size_t i = 0; i < depth_counts.size(); ++i) {
        depth[i] = depth_counts[i];
    }

    return DataFrame::create(Named("depth") = depth);
}
