#include <Rcpp.h>
#include "helpers.h"
#include "topological_sweep.h"
#include <cmath>
using namespace Rcpp;

namespace {

enum class PointColor {
    red,
    blue,
    query
};

inline double safe_project(double numerator, double denominator) {
    if (std::fabs(denominator) < EPS) {
        denominator = (denominator >= 0.0) ? EPS : -EPS;
    }
    return numerator / denominator;
}

long long count_origin_containing_triangles(const std::vector<Point2D>& ordered_points) {
    const int m = static_cast<int>(ordered_points.size());
    if (m < 3) return 0LL;

    std::vector<double> angles(static_cast<std::size_t>(m));
    for (int i = 0; i < m; i++) {
        angles[static_cast<std::size_t>(i)] = atan2pi(ordered_points[static_cast<std::size_t>(i)]);
    }

    long long non_containing_triangles = 0;
    int k = 1;
    for (int i = 0; i < m; i++) {
        if (k < i + 1) k = i + 1;
        while (k < i + m) {
            double angle_diff = angles[static_cast<std::size_t>(k % m)] - angles[static_cast<std::size_t>(i)];
            if (angle_diff < 0.0) angle_diff += 2 * M_PI;
            if (angle_diff >= M_PI - EPS) break;
            k++;
        }
        const int num_in_half_plane = k - i - 1;
        non_containing_triangles += combinations(num_in_half_plane, 2);
    }

    return combinations(m, 3) - non_containing_triangles;
}

long long spherical_asd_toposweep(Point3D ray, const std::vector<Point3D>& points) {
    if (magnitude_sq(ray) <= EPS) return 0LL;

    Point3D ray_unit = normalize(ray);
    Point2D ray_projected = {safe_project(ray_unit.x, ray_unit.z), safe_project(ray_unit.y, ray_unit.z)};
    const bool ray_positive_hemisphere = ray_unit.z > 0.0;

    std::vector<Point2D> projected_points;
    std::vector<PointColor> colors;
    std::vector<simplicialdepth::SweepLine> lines;
    projected_points.reserve(points.size() + 1);
    colors.reserve(points.size() + 1);
    lines.reserve(points.size() + 1);

    for (const auto& p_raw : points) {
        if (magnitude_sq(p_raw) <= EPS) continue;
        const Point3D p_unit = normalize(p_raw);
        const Point2D projected = {safe_project(p_unit.x, p_unit.z), safe_project(p_unit.y, p_unit.z)};
        projected_points.push_back(projected);
        lines.push_back({projected.y, projected.x});

        const bool same_hemisphere = (p_unit.z > 0.0) == ray_positive_hemisphere;
        colors.push_back(same_hemisphere ? PointColor::blue : PointColor::red);
    }

    const int n = static_cast<int>(projected_points.size());
    if (n < 3) return 0LL;

    const int query_index = n;
    projected_points.push_back(ray_projected);
    colors.push_back(PointColor::query);
    lines.push_back({ray_projected.y, ray_projected.x});

    simplicialdepth::TopologicalSweep topological_sweep;
    const std::vector<simplicialdepth::SweepIntersectionPair> swept_pairs = topological_sweep.sweep(lines);
    const std::vector<int> slope_orders = simplicialdepth::build_order_matrix(swept_pairs, n + 1);

    std::vector<int> angle_order_around_query;
    angle_order_around_query.reserve(static_cast<std::size_t>(n));
    for (int i = 0; i < n; i++) {
        const int s_index = slope_orders[static_cast<std::size_t>(query_index * n + i)];
        if (projected_points[static_cast<std::size_t>(s_index)].y - ray_projected.y >= 0.0) {
            angle_order_around_query.push_back(s_index);
        }
    }
    for (int i = 0; i < n; i++) {
        const int s_index = slope_orders[static_cast<std::size_t>(query_index * n + i)];
        if (projected_points[static_cast<std::size_t>(s_index)].y - ray_projected.y < 0.0) {
            angle_order_around_query.push_back(s_index);
        }
    }

    std::vector<Point2D> ordered_red_points;
    ordered_red_points.reserve(static_cast<std::size_t>(n));
    int n_red = 0;
    for (int i = 0; i < n; i++) {
        const int sorted_index = angle_order_around_query[static_cast<std::size_t>(i)];
        if (colors[static_cast<std::size_t>(sorted_index)] == PointColor::red) {
            ordered_red_points.push_back(projected_points[static_cast<std::size_t>(sorted_index)] - ray_projected);
            n_red++;
        }
    }
    const int n_blue = n - n_red;

    const long long N1 = count_origin_containing_triangles(ordered_red_points);

    long long N2 = 0;
    for (int i = 0; i < n; i++) {
        if (colors[static_cast<std::size_t>(i)] == PointColor::blue) continue;

        const bool reference_label = (ray_projected.y - projected_points[static_cast<std::size_t>(i)].y >= 0.0);
        std::vector<unsigned char> labels;
        labels.reserve(static_cast<std::size_t>(n_blue));
        int blue_counter = 0;
        int rank_of_reference_point = 0;

        for (int j = 0; j < n; j++) {
            const int s_index = slope_orders[static_cast<std::size_t>(i * n + j)];
            if (colors[static_cast<std::size_t>(s_index)] == PointColor::blue) {
                const bool label =
                    ((projected_points[static_cast<std::size_t>(s_index)].y -
                      projected_points[static_cast<std::size_t>(i)].y >= 0.0) == reference_label);
                labels.push_back(static_cast<unsigned char>(label ? 1 : 0));
                blue_counter += 1;
            } else if (colors[static_cast<std::size_t>(s_index)] == PointColor::query) {
                rank_of_reference_point = blue_counter;
            }
        }

        int zeros_left = 0;
        int zeros_right = 0;
        for (int j = 0; j < rank_of_reference_point; j++) {
            if (labels[static_cast<std::size_t>(j)] == 0U) zeros_left += 1;
        }
        for (int j = rank_of_reference_point; j < n_blue; j++) {
            if (labels[static_cast<std::size_t>(j)] == 0U) zeros_right += 1;
        }
        N2 += static_cast<long long>(zeros_left) * static_cast<long long>(zeros_right);

        long long a = 0;
        long long b = 0;
        for (int j = 0; j < rank_of_reference_point; j++) {
            if (labels[static_cast<std::size_t>(j)] == 1U) {
                a += 1;
            } else {
                b += a;
            }
        }
        N2 += b;

        a = 0;
        b = 0;
        for (int j = rank_of_reference_point; j < n_blue; j++) {
            if (labels[static_cast<std::size_t>(j)] == 0U) {
                a += 1;
            } else {
                b += a;
            }
        }
        N2 += b;
    }

    long long N3 = 0;
    for (int i = 0; i < n; i++) {
        if (colors[static_cast<std::size_t>(i)] == PointColor::red) continue;

        const Point2D blue_point_centered = projected_points[static_cast<std::size_t>(i)] - ray_projected;
        const Point2D label_vector = {blue_point_centered.y, -blue_point_centered.x};

        int rank_of_blue_point = 0;
        while (rank_of_blue_point < n &&
               slope_orders[static_cast<std::size_t>(query_index * n + rank_of_blue_point)] != i) {
            rank_of_blue_point += 1;
        }

        int rank_of_origin = 0;
        while (rank_of_origin < n &&
               slope_orders[static_cast<std::size_t>(i * n + rank_of_origin)] != query_index) {
            rank_of_origin += 1;
        }

        long long a = 0;
        long long b = 0;
        long long n_ones = 0;
        long long n_zeros = 0;
        for (int j = 0; j < n; j++) {
            const int s_index =
                slope_orders[static_cast<std::size_t>(query_index * n + (j + rank_of_blue_point) % n)];
            if (colors[static_cast<std::size_t>(s_index)] == PointColor::blue) continue;

            if (dot(projected_points[static_cast<std::size_t>(s_index)] - ray_projected, label_vector) >= 0.0) {
                a += 1;
                n_ones += 1;
            } else {
                b += a;
                n_zeros += 1;
            }
        }
        const long long trash_right = b;

        a = 0;
        b = 0;
        for (int j = 0; j < n; j++) {
            const int s_index = slope_orders[static_cast<std::size_t>(i * n + (j + rank_of_origin) % n)];
            if (colors[static_cast<std::size_t>(s_index)] == PointColor::blue) continue;

            if (dot(projected_points[static_cast<std::size_t>(s_index)] - ray_projected, label_vector) < 0.0) {
                a += 1;
            } else {
                b += a;
            }
        }
        const long long trash_left = b;
        N3 += n_zeros * n_ones - trash_left - trash_right;
    }

    return N1 + N2 + N3;
}

} // namespace

// [[Rcpp::export]]
long long spherical_asd(NumericMatrix X, NumericVector ray) {
    int n = X.nrow();
    // if (n < 3) return 0LL;
    if (n < 3) stop("n must be at least 3");
    if (X.ncol() != 3 || ray.size() != 3) {
        stop("P must be n x 3 and ray must be of length 3");
    }

    // Convert matrix to std::vector<Point3D>
    std::vector<Point3D> P(n);
    // for (int i = 0; i < n; i++)
    //     P.push_back({X(i,0), X(i,1), X(i,2)});
    for (int i = 0; i < n; i++) {
        P[i].x = X(i, 0);
        P[i].y = X(i, 1);
        P[i].z = X(i, 2);
    }

    // Convert vector to Point3D
    Point3D rayray = { ray[0], ray[1], ray[2] };

    return spherical_asd_toposweep(rayray, P);
}

// Compute leave-one-out spherical ASD for every row of X (each row is ray; X without that row is sample).
// This is the fast path for x = NULL wrappers: parallelize across i to avoid per-point thread startup.
// [[Rcpp::export]]
NumericVector spherical_asd_all_points(NumericMatrix X, int threads = 0) {
    const int n = X.nrow();
    if (n < 3) {
        return NumericVector(n);
    }
    if (X.ncol() != 3) {
        stop("X must be n x 3");
    }

    std::vector<Point3D> Pfull(static_cast<std::size_t>(n));
    for (int i = 0; i < n; i++) {
        Pfull[static_cast<std::size_t>(i)].x = X(i, 0);
        Pfull[static_cast<std::size_t>(i)].y = X(i, 1);
        Pfull[static_cast<std::size_t>(i)].z = X(i, 2);
    }

    NumericVector out(n);
    RcppParallel::RVector<double> out_view(out);

    struct Worker : public RcppParallel::Worker {
        const std::vector<Point3D>& P;
        const int n;
        RcppParallel::RVector<double> out;

        Worker(const std::vector<Point3D>& P_, int n_, RcppParallel::RVector<double> out_)
            : P(P_), n(n_), out(out_) {}

        void operator()(std::size_t begin, std::size_t end) {
            std::vector<Point3D> sample;
            sample.reserve(static_cast<std::size_t>(n - 1));

            for (std::size_t ii = begin; ii < end; ++ii) {
                const int i = static_cast<int>(ii);
                sample.clear();
                for (int j = 0; j < n; ++j) {
                    if (j == i) continue;
                    sample.push_back(P[static_cast<std::size_t>(j)]);
                }
                out[ii] = static_cast<double>(spherical_asd_toposweep(P[static_cast<std::size_t>(i)], sample));
            }
        }
    };

    Worker worker(Pfull, n, out_view);
    simplicialdepth::parallel_for(0, static_cast<std::size_t>(n), worker, threads);
    return out;
}
