#include <Rcpp.h>
#include "helpers.h"
#include "simplicial_depth_2d.h"
#include <cmath>
#include <cassert>
using namespace Rcpp;

// The main algorithm for fast 3D angular simplicial depth
long long spherical_asd(Point3D ray, const std::vector<Point3D>& P) {
    int n = P.size();
    if (n < 3) return 0LL;

    // Step 1: Project points to a unit sphere.
    std::vector<Point3D> P_prime;
    P_prime.reserve(n);
    for (const auto& pi : P) {
        if (magnitude_sq(pi) > EPS) {
            P_prime.push_back(normalize(pi));
        }
    }
    n = P_prime.size();
    if (n < 3) return 0;

    bool ray_color = (ray.z > 0);
    // Project ray onto the z=1 plane
    Point2D ray_projected = {ray.x / ray.z, ray.y / ray.z};
    // Project all the other points onto this plane
    std::vector<Point2D> P_projected(n);
    for (int i=0; i<n; i++) P_projected[i] = {P_prime[i].x / P_prime[i].z, P_prime[i].y / P_prime[i].z};
    // Find the color of each point (is it in the same half space as ray?)
    std::vector<bool> P_colors(n);
    for (int i=0; i<n; i++) P_colors[i] = ((P_prime[i].z > 0) == ray_color);
    // Shift all the points so that the projected ray becomes the origin
    for (int i=0; i<n; i++) P_projected[i] = P_projected[i] - ray_projected;

    std::vector<Point2D> P_red;
    std::vector<Point2D> P_blue;
    P_red.reserve(n);
    P_blue.reserve(n);
    for (int i=0; i<n; i++) {
        if (P_colors[i] == 0) P_blue.push_back(P_projected[i]);
        else P_red.push_back(P_projected[i]);
    }
    int n_blue_points = P_blue.size();
    int n_red_points = P_red.size();

    
    long long N1 = simplicial_depth_2d({0.0, 0.0}, P_red); // Number of triangles of the opposite color of ray that contain the origin

    // Helper struct for the rotational sweep algorithm
    struct ProjPoint {
        Point2D p_2d;
        double slope; // Angle on the projection plane
        bool label;
        int index;
    };

    long long N2 = 0;
    for (const auto& red_point : P_red) {
        Point2D center_shifted = {-red_point.x, -red_point.y};
        std::vector<ProjPoint> P_blue_shifted(n_blue_points);
        bool reference_label = center_shifted.y > 0 ? 1 : 0;
        double reference_slope = -center_shifted.x / center_shifted.y;
        for (int j=0; j<n_blue_points; j++) {
            P_blue_shifted[j].p_2d = P_blue[j] - red_point;
            P_blue_shifted[j].label = ((P_blue_shifted[j].p_2d.y > 0) == reference_label);
            P_blue_shifted[j].slope = -P_blue_shifted[j].p_2d.x / P_blue_shifted[j].p_2d.y;
            P_blue_shifted[j].index = j;
        }

        std::sort(P_blue_shifted.begin(), P_blue_shifted.end(), [](const ProjPoint& a, const ProjPoint& b){return a.slope < b.slope;});

        int reference_order = 0;
        while (reference_order < n_blue_points) {
            if (P_blue_shifted[reference_order].slope < reference_slope) reference_order += 1;
            else break;
        }

        // Case 1: 0 *1* 0
        // Count zeros to the left of reference order, to the right of ref. order and then multiply them
        int s_left = 0;
        int s_right = 0;
        for (int j=0; j<reference_order; j++) {
            if (P_blue_shifted[j].label == 0) s_left++;
        }
        for (int j=reference_order; j<n_blue_points; j++) {
            if (P_blue_shifted[j].label == 0) s_right++;
        }
        long long case_1 = s_left * s_right;


        // Case 2: 1 0 *1*
        // Counting number of 1 0 subsequences to the left of reference_order
        long long case_2 = 0, ones = 0;
        for (int j=0; j<reference_order; j++) {
            if (P_blue_shifted[j].label == 1) ones++;
            else case_2 += ones;
        }

        // Case 2: *1* 0 1
        // Counting number of 1 0 subsequences to the left of reference_order
        long long case_3 = 0, zeros = 0;
        for (int j=reference_order; j<n_blue_points; j++) {
            if (P_blue_shifted[j].label == 0) zeros++;
            else case_3 += zeros;
        }

        N2 = N2 + case_1 + case_2 + case_3;
    }

    // Number of intersections between all red segments with blue segments whose one endpoint is the origin
    long long N3 = 0;
    for (const auto& blue_point : P_blue) {
        // Rotate all red points counterclockwise so that the blue point Blue[b,] ends up on the negative x-axis
        double theta = atan2(blue_point.y, blue_point.x);
        double alpha = M_PI - theta;
        Matrix2D rotation_matrix = {{cos(alpha), -sin(alpha)}, {sin(alpha), cos(alpha)}};

        std::vector<ProjPoint> P_red_rotated(n_red_points);
        int red_labels_ones = 0, red_labels_zeros = 0;
        for (int j=0; j<n_red_points; j++) {
            P_red_rotated[j].p_2d = rotation_matrix * P_red[j];
            P_red_rotated[j].slope = -P_red_rotated[j].p_2d.x / P_red_rotated[j].p_2d.y;
            P_red_rotated[j].label = (P_red_rotated[j].p_2d.y > 0);
            red_labels_ones += (P_red_rotated[j].label == 1);
            red_labels_zeros += (P_red_rotated[j].label == 0);
        }
        Point2D blue_point_rotated = rotation_matrix * blue_point;
        std::sort(P_red_rotated.begin(), P_red_rotated.end(), [](const ProjPoint& a, const ProjPoint& b) {return a.slope < b.slope;});
        
        long long total_right = 0, ones = 0;
        for (const auto& red_point : P_red_rotated) {
            if (red_point.label == 1) ones++;
            else total_right += ones;
        }

        for (int j=0; j<n_red_points; j++) {
            P_red_rotated[j].p_2d = P_red_rotated[j].p_2d - blue_point_rotated;
            P_red_rotated[j].slope = -P_red_rotated[j].p_2d.x / P_red_rotated[j].p_2d.y;
        }
        std::sort(P_red_rotated.begin(), P_red_rotated.end(), [](const ProjPoint& a, const ProjPoint& b) {return a.slope < b.slope;});

        long long total_left = 0, zeros = 0;
        for (const auto& red_point : P_red_rotated) {
            if (red_point.label == 0) zeros++;
            else total_left += zeros;
        }
        N3 = N3 + red_labels_ones * red_labels_zeros - total_right - total_left;
    }
    return N1 + N2 + N3;
}

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

    return spherical_asd(rayray, P);
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
                out[ii] = static_cast<double>(spherical_asd(P[static_cast<std::size_t>(i)], sample));
            }
        }
    };

    Worker worker(Pfull, n, out_view);
    simplicialdepth::parallel_for(0, static_cast<std::size_t>(n), worker, threads);
    return out;
}
