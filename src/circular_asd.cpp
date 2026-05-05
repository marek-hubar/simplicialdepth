#include <Rcpp.h>
#include "helpers.h"
#include <cassert>
#include <random>
#include <chrono>
using namespace Rcpp;

// Helper struct for sorting 2D points by angle
struct Point2D_indexed {
    Point2D p;
    double angle;
    double slope;
    double angle_mod_PI;
    bool label;
    int original_index;
};

// An O(n log(n)) algorithm for angular simplicial depth
// Very similar to the regular algorithm, but the sweeping line only goes until it meets the -ray
long long circular_asd(Point2D ray, const std::vector<Point2D>& points) {
    int m = points.size();
    if (m < 2) return 0;

    double ray_angle = atan2(ray.y, ray.x);
    ray_angle = ray_angle < 0 ? ray_angle + 2 * M_PI : ray_angle;

    std::vector<Point2D_indexed> sorted_2d_points(m);
    for(int i=0; i<m; ++i) {
        sorted_2d_points[i].p = points[i] - ray;
        double angle = atan2(points[i].y, points[i].x);
        sorted_2d_points[i].angle = angle < 0 ? angle + 2 * M_PI : angle;
        sorted_2d_points[i].angle -= ray_angle;
        sorted_2d_points[i].angle = sorted_2d_points[i].angle < 0 ? sorted_2d_points[i].angle + 2 * M_PI : sorted_2d_points[i].angle;
    }

    double ray_antipodal_angle = atan2(-ray.y, -ray.x);
    ray_antipodal_angle = ray_antipodal_angle < 0 ? ray_antipodal_angle + 2 * M_PI : ray_antipodal_angle;

    // std::sort(sorted_2d_points.begin(), sorted_2d_points.end(), [](const Point2D_indexed& a, const Point2D_indexed& b){
    //     return a.angle < b.angle;
    // });
    parallel_sort(sorted_2d_points, [](const Point2D_indexed& a, const Point2D_indexed& b){
        return a.angle < b.angle;
    });

    long long non_containing_triangles = 0;
    int k = 1, i = 0; // Two-pointer for finding the antipodal point
    for (i = 0; i < m && sorted_2d_points[i].angle < M_PI; ++i) {
        if (k < i + 1) k = i + 1;
        while (k < i + m) {
            double angle_diff = sorted_2d_points[k % m].angle - sorted_2d_points[i].angle;
            if (angle_diff < 0) angle_diff += 2 * M_PI;
            if (angle_diff >= M_PI - EPS) break;
            k++;
        }
        int num_in_half_plane = k - i - 1;
        non_containing_triangles += combinations(num_in_half_plane, 1);
    }
    while (k < i + m + 10) {
        double angle_diff = sorted_2d_points[k % m].angle - M_PI;
        if (angle_diff < 0) angle_diff += 2 * M_PI;
        if (angle_diff >= M_PI - EPS) break;
        k++;
    }
    int num_in_half_plane = k - i;
    non_containing_triangles += combinations(num_in_half_plane, 2);

    return combinations(m, 2) - non_containing_triangles;
}

// An O(n log(n)) algorithm for angular simplicial depth based on the 01 sequence method
long long circular_asd_01(Point2D ray, const std::vector<Point2D>& points) {
    int m = points.size();
    if (m < 2) return 0;

    double ray_angle = atan2(ray.y, ray.x);
    ray_angle = ray_angle < 0 ? ray_angle + 2 * M_PI : ray_angle;

    int N0 = 0, N1 = 0;

    std::vector<Point2D_indexed> sorted_2d_points(m);
    for(int i=0; i<m; ++i) {
        sorted_2d_points[i].p = points[i];
        // Find the polar angle with respect to the ray
        sorted_2d_points[i].angle = atan2pi(points[i]) - ray_angle;
        sorted_2d_points[i].angle = sorted_2d_points[i].angle < 0 ? sorted_2d_points[i].angle + 2 * M_PI : sorted_2d_points[i].angle;
        sorted_2d_points[i].label = (sorted_2d_points[i].angle < M_PI);
        if (sorted_2d_points[i].label == 0) N0++;
        else N1++;
        // Angles modulo PI
        sorted_2d_points[i].angle_mod_PI = sorted_2d_points[i].angle < M_PI ? sorted_2d_points[i].angle : sorted_2d_points[i].angle - M_PI;
    }

    // std::sort(sorted_2d_points.begin(), sorted_2d_points.end(), [](const Point2D_indexed& a, const Point2D_indexed& b){
    //     return a.angle_mod_PI < b.angle_mod_PI;
    // });
    parallel_sort(sorted_2d_points, [](const Point2D_indexed& a, const Point2D_indexed& b){
        return a.angle_mod_PI < b.angle_mod_PI;
    });

    long long subsequences_10 = 0;
    // Count the number of "10" subsequences
    for (const auto& p : sorted_2d_points) {
        if (p.label == 1) subsequences_10 += N0;
        else N0--;
    }
    return subsequences_10;
}

// [[Rcpp::export]]
long long circular_asd(NumericMatrix X, NumericVector x) {
    int n = X.nrow();
    std::vector<Point2D> P(n);
    for (int i = 0; i < n; i++) {
        P[i].x = X(i, 0);
        P[i].y = X(i, 1);
    }
    Point2D xx = {x[0], x[1]};
    return circular_asd_01(xx, P);
}

// An O(n log(n)) algorithm for angular simplicial depth based on the 01 sequence method
//
// For each point r in the sample, we want the same quantity returned by:
//   circular_asd_01(r, points \ {r})
// but without re-sorting for every r.
//
// Key observation: sorting by angle modulo PI is equivalent to sorting by (global) angle modulo PI,
// i.e. by alpha = atan2(p) mod PI. Relative angles for different rays induce only a cyclic shift.
std::vector<long long> circular_asd_01_all_points(const std::vector<Point2D>& points) {
    const int m = static_cast<int>(points.size());
    std::vector<long long> depths(m, 0);
    if (m <= 2) return depths;

    struct Entry {
        double alpha;   // angle modulo PI in [0, PI)
        int upper;      // 1 if atan2(p) in [0, PI), else 0
        int original_index;
    };

    std::vector<Entry> entries(m);
    for (int i = 0; i < m; ++i) {
        const double theta = atan2pi(points[i]);
        const int upper = (theta < M_PI) ? 1 : 0;
        double alpha = theta;
        if (alpha >= M_PI) alpha -= M_PI;
        entries[i] = {alpha, upper, i};
    }

    // Sort by alpha (equivalent to sorting by -x/y), tie-breaking deterministically.
    parallel_sort(entries, [](const Entry& a, const Entry& b) {
        if (a.alpha < b.alpha) return true;
        if (a.alpha > b.alpha) return false;
        if (a.upper != b.upper) return a.upper > b.upper; // put upper first
        return a.original_index < b.original_index;
    });

    const int n = m;
    std::vector<int> U(n);
    std::vector<int> orig(n);
    for (int i = 0; i < n; ++i) {
        U[i] = entries[i].upper;
        orig[i] = entries[i].original_index;
    }

    // Prefix counts and inversions on U.
    std::vector<long long> ones_prefix(n + 1, 0), zeros_prefix(n + 1, 0);
    std::vector<long long> inv10_prefix(n + 1, 0), inv01_prefix(n + 1, 0);
    for (int i = 0; i < n; ++i) {
        ones_prefix[i + 1] = ones_prefix[i] + U[i];
        zeros_prefix[i + 1] = zeros_prefix[i] + (1 - U[i]);

        // inv10: count of pairs (i<j) with U[i]=1, U[j]=0.
        inv10_prefix[i + 1] = inv10_prefix[i] + (U[i] == 0 ? ones_prefix[i] : 0);
        // inv01: count of pairs (i<j) with U[i]=0, U[j]=1.
        inv01_prefix[i + 1] = inv01_prefix[i] + (U[i] == 1 ? zeros_prefix[i] : 0);
    }

    // Suffix counts and inversions on U.
    std::vector<long long> ones_suffix(n + 1, 0), zeros_suffix(n + 1, 0);
    std::vector<long long> inv10_suffix(n + 1, 0), inv01_suffix(n + 1, 0);
    for (int i = n - 1; i >= 0; --i) {
        ones_suffix[i] = ones_suffix[i + 1] + U[i];
        zeros_suffix[i] = zeros_suffix[i + 1] + (1 - U[i]);

        inv10_suffix[i] = inv10_suffix[i + 1] + (U[i] == 1 ? zeros_suffix[i + 1] : 0);
        inv01_suffix[i] = inv01_suffix[i + 1] + (U[i] == 0 ? ones_suffix[i + 1] : 0);
    }

    // For ray at position k in the alpha-order, the angle_mod_PI order for the remaining points is:
    //   suffix (k+1..n-1) followed by prefix (0..k-1).
    // Labels (relative to the ray) are:
    //   - for alpha >= alpha_ray: label = (upper == upper_ray)
    //   - for alpha <  alpha_ray: label = (upper != upper_ray)
    // Under general position (no alpha ties), this yields closed-form expressions below.
    for (int k = 0; k < n; ++k) {
        const int upper_ray = U[k];

        const long long ones_pref = ones_prefix[k];
        const long long zeros_pref = zeros_prefix[k];
        const long long ones_suf = ones_suffix[k + 1];
        const long long zeros_suf = zeros_suffix[k + 1];

        long long N10 = 0;
        if (upper_ray == 1) {
            // suffix labels = U, prefix labels = 1-U
            N10 = inv10_suffix[k + 1] + inv01_prefix[k] + ones_suf * ones_pref;
        } else {
            // suffix labels = 1-U, prefix labels = U
            N10 = inv01_suffix[k + 1] + inv10_prefix[k] + zeros_suf * zeros_pref;
        }

        depths[orig[k]] = N10;
    }

    return depths;
}

// std::vector<long long> circular_asd_all_points(const std::vector<Point2D>& X) {
//     int n = X.size();
//     std::vector<long long> depths(n);
//     if (n <= 2) return depths;
//
//     depths[0] = circular_asd_01(X[0], std::vector<Point2D>(X.begin() + 1, X.end()));
//
//     return depths;
// }

// [[Rcpp::export]]
DataFrame circular_asd_all_arcs(NumericMatrix points_mat) {
    int n = points_mat.nrow();
    std::vector<Point2D> points_vec(n);

    // Convert matrix to std::vector<Point2D>
    for (int i = 0; i < n; i++) {
        points_vec[i].x = points_mat(i, 0);
        points_vec[i].y = points_mat(i, 1);
    }

    std::vector<long long> arcs = circular_asd_01_all_points(points_vec);

    // Convert std::vector<Arc> to DataFrame
    int m = arcs.size();
    NumericVector start_angles(m), end_angles(m), right_points_x(m), right_points_y(m);
    IntegerVector depths(m), right_point_depths(m);
    for (int i = 0; i < m; i++) {
        depths[i] = arcs[i];//.depth;
    }

    return DataFrame::create(
        Named("depth") = depths//,
    );
}
