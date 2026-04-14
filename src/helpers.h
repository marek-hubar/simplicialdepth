#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "rcpp_parallel_utils.h"

#if __has_include(<tbb/parallel_sort.h>)
    #include <tbb/parallel_sort.h>
    #define SIMPLICIALDEPTH_HAS_TBB_PARALLEL_SORT 1
#else
    #define SIMPLICIALDEPTH_HAS_TBB_PARALLEL_SORT 0
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const double EPS = 1e-12;



template<typename T, typename Compare>
void parallel_sort(std::vector<T>& vec, Compare comp) {
    const std::size_t n = vec.size();
    if (n <= 1000) {
        std::sort(vec.begin(), vec.end(), comp);
        return;
    }

#if SIMPLICIALDEPTH_HAS_TBB_PARALLEL_SORT
    tbb::parallel_sort(vec.begin(), vec.end(), comp);
#else
    std::sort(vec.begin(), vec.end(), comp);
#endif
}

inline int mod(int k, int n) {
    if (k >= 0) return k % n;
    return mod(n + k, n);
}

// Represents a point or a vector in 2D space
struct Point2D {
    double x, y;
    
    Point2D operator+(const Point2D& other) const {
        return {x + other.x, y + other.y};
    }
    Point2D operator-(const Point2D& other) const {
        return {x - other.x, y - other.y};
    }
};

// Overload the << operator for std::ostream
inline std::ostream& operator<<(std::ostream& os, const Point2D& point) {
    os << "(" << point.x << ", " << point.y << ")";
    return os;
}

// --- Vector Operations ---
inline double dot(const Point2D& v1, const Point2D& v2) { return v1.x * v2.x + v1.y * v2.y; }
inline double magnitude_sq(const Point2D& v) { return dot(v, v); }
inline double magnitude(const Point2D& v) { return std::sqrt(magnitude_sq(v)); }
inline Point2D normalize(const Point2D& v) {
    double mag = magnitude(v);
    if (mag < EPS) return {0, 0};
    return {v.x / mag, v.y / mag};
}

inline double atan2pi(double y, double x) {
    double angle = atan2(y, x);
    angle = angle < 0 ? angle + 2 * M_PI : angle;
    return angle;
}

inline double atan2pi(Point2D p) {
    return atan2pi(p.y, p.x);
}

inline double oposite_angle(double angle) {
    return angle < M_PI ? angle + M_PI : angle - M_PI;
}

struct Point3D {
    double x, y, z;
    
    Point3D operator+(const Point3D& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }
    Point3D operator-(const Point3D& other) const {
        return {x - other.x, y - other.y, z - other.z};
    }
    Point3D operator*(const Point3D& other) const {
        return {x - other.x, y - other.y, z - other.z};
    }
};
inline Point3D operator*(double scalar, const Point3D& point) {
    return {scalar * point.x, scalar * point.y, scalar * point.z};
}
// Overload the << operator for std::ostream
inline std::ostream& operator<<(std::ostream& os, const Point3D& point) {
    os << "(" << point.x << ", " << point.y << ", " << point.z << ")";
    return os;
}

// --- Vector Operations ---
inline double dot(const Point3D& v1, const Point3D& v2) { return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z; }
inline double magnitude_sq(const Point3D& v) { return dot(v, v); }
inline double magnitude(const Point3D& v) { return std::sqrt(magnitude_sq(v)); }

inline Point3D normalize(const Point3D& v) {
    double mag = magnitude(v);
    if (mag < EPS) return {0, 0};
    return {v.x / mag, v.y / mag, v.z / mag};
}
inline Point3D cross(const Point3D& v1, const Point3D& v2) {
    return { v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x };
}

// Safe sign helper returning -1, 0, or +1
static inline int sgn(double v) {
    return (v > 0) - (v < 0);
}

// Computes n choose k for k <= 3
inline long long combinations(int n, int k) {
    if (k < 0 || k > n) return 0;
    if (k == 0 || k == n) return 1;
    if (k > n / 2) k = n - k; // use symmetry

    long long result = 1;
    for (int i = 1; i <= k; ++i) {
        result = result * (n - k + i) / i;
    }
    return result;
}

struct Arc {
    double left_angle;
    double right_angle;
    long long depth;
    long long right_point_depth;
    Point2D right_point;

    double mid_point() const {
        double mid = (left_angle + right_angle) / 2.0;
        if (std::fabs(left_angle - right_angle) > M_PI)
            mid = oposite_angle(mid);
        return mid;
    }
};

struct Matrix2D {
    Point2D first_row, second_row;
};

inline Point2D operator*(const Matrix2D& matrix, const Point2D& point) {
    return {dot(matrix.first_row, point), dot(matrix.second_row, point)};
}
