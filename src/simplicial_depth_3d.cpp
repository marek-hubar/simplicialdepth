#include <Rcpp.h>
#include "helpers.h"
#include "topological_sweep.h"

#include <cmath>
#include <vector>

using namespace Rcpp;

namespace {

enum class PointColor { red, blue };

inline double safe_project(double numerator, double denominator) {
    if (std::fabs(denominator) < EPS) {
        denominator = (denominator >= 0.0) ? EPS : -EPS;
    }
    return numerator / denominator;
}

// Counts triangles (from points given in circular angular order around origin) that contain the origin.
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

inline void build_angle_order_around(int pivot,
                                    const std::vector<Point2D>& points,
                                    const std::vector<int>& slope_orders,
                                    std::vector<int>& out) {
    const int n = static_cast<int>(points.size());
    const int row_width = n - 1;
    out.clear();
    out.reserve(static_cast<std::size_t>(row_width));

    const std::size_t row_base = static_cast<std::size_t>(pivot * row_width);
    const Point2D pivot_pt = points[static_cast<std::size_t>(pivot)];

    for (int j = 0; j < row_width; j++) {
        const int idx = slope_orders[row_base + static_cast<std::size_t>(j)];
        if (points[static_cast<std::size_t>(idx)].y - pivot_pt.y >= 0.0) {
            out.push_back(idx);
        }
    }
    for (int j = 0; j < row_width; j++) {
        const int idx = slope_orders[row_base + static_cast<std::size_t>(j)];
        if (points[static_cast<std::size_t>(idx)].y - pivot_pt.y < 0.0) {
            out.push_back(idx);
        }
    }
}

long long compute_N1_red_inside_blue_triangles(const std::vector<Point2D>& points,
                                              const std::vector<PointColor>& colors,
                                              const std::vector<int>& slope_orders) {
    const int n = static_cast<int>(points.size());
    if (n < 4) return 0LL;

    long long N1 = 0;
    std::vector<int> angle_order;
    std::vector<Point2D> ordered_vectors;

    for (int pivot = 0; pivot < n; ++pivot) {
        if (colors[static_cast<std::size_t>(pivot)] != PointColor::red) continue;

        build_angle_order_around(pivot, points, slope_orders, angle_order);

        ordered_vectors.clear();
        ordered_vectors.reserve(angle_order.size());
        for (int idx : angle_order) {
            if (colors[static_cast<std::size_t>(idx)] == PointColor::blue) {
                ordered_vectors.push_back(points[static_cast<std::size_t>(idx)] -
                                          points[static_cast<std::size_t>(pivot)]);
            }
        }

        N1 += count_origin_containing_triangles(ordered_vectors);
    }

    return N1;
}

long long compute_N2_blue_inside_red_triangles(const std::vector<Point2D>& points,
                                              const std::vector<PointColor>& colors,
                                              const std::vector<int>& slope_orders) {
    const int n = static_cast<int>(points.size());
    if (n < 4) return 0LL;

    long long N2 = 0;
    std::vector<int> angle_order;
    std::vector<Point2D> ordered_vectors;

    for (int pivot = 0; pivot < n; ++pivot) {
        if (colors[static_cast<std::size_t>(pivot)] != PointColor::blue) continue;

        build_angle_order_around(pivot, points, slope_orders, angle_order);

        ordered_vectors.clear();
        ordered_vectors.reserve(angle_order.size());
        for (int idx : angle_order) {
            if (colors[static_cast<std::size_t>(idx)] == PointColor::red) {
                ordered_vectors.push_back(points[static_cast<std::size_t>(idx)] -
                                          points[static_cast<std::size_t>(pivot)]);
            }
        }

        N2 += count_origin_containing_triangles(ordered_vectors);
    }

    return N2;
}

// Implements the O(n^2) red-blue segment crossing count from src/sd_3d.md.
long long compute_N3_red_blue_segment_crossings(const std::vector<Point2D>& points,
                                               const std::vector<PointColor>& colors,
                                               const std::vector<int>& slope_orders) {
    const int n = static_cast<int>(points.size());
    if (n < 4) return 0LL;

    int r = 0, b = 0;
    for (int i = 0; i < n; ++i) {
        if (colors[static_cast<std::size_t>(i)] == PointColor::red)
            ++r;
        else
            ++b;
    }
    if (r < 2 || b < 2) return 0LL;

    const long long total_pairs = combinations(r, 2) * combinations(b, 2);

    long long S_non_intersecting = 0;

    const int m = n - 1;
    std::vector<int> angle_order;
    std::vector<int> idx2;
    std::vector<double> ang;
    std::vector<double> ang2;

    angle_order.reserve(static_cast<std::size_t>(m));
    idx2.resize(static_cast<std::size_t>(2 * m));
    ang.resize(static_cast<std::size_t>(m));
    ang2.resize(static_cast<std::size_t>(2 * m));

    for (int base = 0; base < n; ++base) {
        if (colors[static_cast<std::size_t>(base)] != PointColor::red) continue;

        build_angle_order_around(base, points, slope_orders, angle_order);

        const Point2D base_pt = points[static_cast<std::size_t>(base)];
        for (int j = 0; j < m; ++j) {
            const int idx = angle_order[static_cast<std::size_t>(j)];
            ang[static_cast<std::size_t>(j)] = atan2pi(points[static_cast<std::size_t>(idx)] - base_pt);
        }
        // Ensure angles are non-decreasing in this circular order by unwrapping.
        for (int j = 1; j < m; ++j) {
            while (ang[static_cast<std::size_t>(j)] < ang[static_cast<std::size_t>(j - 1)]) {
                ang[static_cast<std::size_t>(j)] += 2 * M_PI;
            }
        }

        for (int j = 0; j < m; ++j) {
            const int idx = angle_order[static_cast<std::size_t>(j)];
            idx2[static_cast<std::size_t>(j)] = idx;
            idx2[static_cast<std::size_t>(j + m)] = idx;

            ang2[static_cast<std::size_t>(j)] = ang[static_cast<std::size_t>(j)];
            ang2[static_cast<std::size_t>(j + m)] = ang[static_cast<std::size_t>(j)] + 2 * M_PI;
        }

        int k = 1;
        long long red_in_window = 0;
        long long blue_in_window = 0;

        for (int i = 0; i < m; ++i) {
            if (k < i + 1) {
                k = i + 1;
                red_in_window = 0;
                blue_in_window = 0;
            }

            while (k < i + m) {
                const double diff = ang2[static_cast<std::size_t>(k)] - ang2[static_cast<std::size_t>(i)];
                if (diff >= M_PI - EPS) break;
                const int idx_add = idx2[static_cast<std::size_t>(k)];
                if (colors[static_cast<std::size_t>(idx_add)] == PointColor::red) {
                    red_in_window += 1;
                } else {
                    blue_in_window += 1;
                }
                ++k;
            }

            const int idx_dir = idx2[static_cast<std::size_t>(i)];
            if (colors[static_cast<std::size_t>(idx_dir)] == PointColor::blue) {
                const long long r_L = red_in_window;
                const long long b_L = blue_in_window;
                const long long b_R = static_cast<long long>(b - 1) - b_L;
                S_non_intersecting += r_L * b_R;
            }

            const int remove_pos = i + 1;
            if (remove_pos < k) {
                const int idx_remove = idx2[static_cast<std::size_t>(remove_pos)];
                if (colors[static_cast<std::size_t>(idx_remove)] == PointColor::red) {
                    red_in_window -= 1;
                } else {
                    blue_in_window -= 1;
                }
            }
        }
    }

    return total_pairs - S_non_intersecting;
}

long long simplicial_depth_3d_toposweep_count(const std::vector<Point3D>& X,
                                              const Point3D& q,
                                              int threads) {
    (void)threads; // currently unused (single-query work is O(n^2) but mostly memory-bound)

    std::vector<Point2D> projected_points;
    std::vector<PointColor> colors;
    std::vector<simplicialdepth::SweepLine> lines;

    projected_points.reserve(X.size());
    colors.reserve(X.size());
    lines.reserve(X.size());

    for (const auto& p : X) {
        const Point3D v = p - q;
        if (magnitude_sq(v) <= EPS) continue;

        Point3D u = normalize(v);

        // Avoid z==0 in gnomonic projection; keep sign for hemisphere classification.
        if (std::fabs(u.z) < EPS) {
            u.z = (u.z >= 0.0) ? EPS : -EPS;
        }

        const PointColor c = (u.z > 0.0) ? PointColor::blue : PointColor::red;
        const Point2D proj = {safe_project(u.x, u.z), safe_project(u.y, u.z)};

        const int idx = static_cast<int>(projected_points.size());
        projected_points.push_back(proj);
        colors.push_back(c);

        // Dual mapping point (x,y) -> line y = m x + b with (m=y, b=x).
        // Add tiny deterministic jitter to m to avoid parallel lines (m ties).
        const double jitter = static_cast<double>(idx + 1) * 1e-15;
        lines.push_back({proj.y + jitter, proj.x});
    }

    const int n = static_cast<int>(projected_points.size());
    if (n < 4) return 0LL;

    simplicialdepth::TopologicalSweep topological_sweep;
    const std::vector<simplicialdepth::SweepIntersectionPair> swept_pairs = topological_sweep.sweep(lines);
    const std::vector<int> slope_orders = simplicialdepth::build_order_matrix(swept_pairs, n);

    const long long N1 = compute_N1_red_inside_blue_triangles(projected_points, colors, slope_orders);
    const long long N2 = compute_N2_blue_inside_red_triangles(projected_points, colors, slope_orders);
    const long long N3 = compute_N3_red_blue_segment_crossings(projected_points, colors, slope_orders);

    return N1 + N2 + N3;
}

} // namespace

// [[Rcpp::export]]
long long simplicial_depth_3d(NumericMatrix X, NumericVector q, int threads = 0) {
    const int n = X.nrow();
    if (n < 4) return 0LL;
    if (X.ncol() != 3 || q.size() != 3) {
        stop("X must be n x 3 and q must be of length 3");
    }

    std::vector<Point3D> P(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        P[static_cast<std::size_t>(i)].x = X(i, 0);
        P[static_cast<std::size_t>(i)].y = X(i, 1);
        P[static_cast<std::size_t>(i)].z = X(i, 2);
    }

    const Point3D qq = {q[0], q[1], q[2]};
    return simplicial_depth_3d_toposweep_count(P, qq, threads);
}
