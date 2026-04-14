#include <Rcpp.h>
using namespace Rcpp;
#include "helpers.h"
#include <RcppParallel.h>
#include <cassert>

long long SDk(const std::vector<Point3D>& Xinput, const Point3D& x, int k) {
    int n = Xinput.size();
    std::vector<Point3D> X;
    X.reserve(n);
    for (const auto& p : Xinput) {
        Point3D q{p.x - x.x, p.y - x.y, p.z - x.z};
        if (fabs(q.y) < EPS) q.y = EPS;
        if (fabs(q.x) > EPS) X.push_back(q);
    }
    n = (int)X.size();
    if (n == 0 || k > n) return 0LL;

    struct HelperStruct {
        double x, y;
        int color, label;
        double sort_key;
    };

    struct PointY {
        double y1, y2;
        int sign;
    };
    std::vector<PointY> Y; Y.reserve(n);
    int n_up = 0;
    for (const auto& p : X) {
        double y1 = p.x / p.z;
        double y2 = p.y / p.z;
        int s = sgn(p.z);
        Y.push_back({y1, y2, s});
        if (s == 1) n_up++;
    }
    int n_down = n - n_up;

    long long SDk = combinations(n, k) - combinations(n_up, k) - combinations(n_down, k);
    if (n == 1) return SDk;

    std::vector<long long> counts(n - 1, 0LL);

    std::vector<HelperStruct> sorted_points(n-1);

    for (int idx = 0; idx < n; idx++) {
        int j = 0;
        for (int i = 0; i < n; i++) {
            if (i != idx) {
                HelperStruct* sorted_point = &sorted_points[j];
                PointY* pointy_i = &Y[i], *pointy_idx = &Y[idx];
                sorted_point->x = pointy_i->y1 - pointy_idx->y1,
                sorted_point->y = pointy_i->y2 - pointy_idx->y2,
                sorted_point->color = (pointy_i->sign != pointy_idx->sign) ? 1 : 0,
                sorted_point->label = (sorted_point->y * pointy_i->sign > 0.0) ? 1 : 0,
                sorted_point->sort_key = -(sorted_point->x / sorted_point->y);
                j++;
            }
        }
        int m = (int)sorted_points.size();
        if (m == 0) continue;

        std::sort(sorted_points.begin(), sorted_points.end(), [](const auto& a, const auto& b) {
            return a.sort_key < b.sort_key;
        });

        int EL = 0;
        for (int j = 1; j < m; j++) EL += sorted_points[j].label;
        counts[EL] += sorted_points[0].color;

        for (int j = 0; j < m - 1; j++) {
            if (sorted_points[j].label == 0) EL++;
            if (sorted_points[j+1].label == 1) EL--;
            counts[EL] += sorted_points[j+1].color;
        }
    }

    long long SUM = 0LL;
    for (int l = 1; l <= n - 1; l++) {
        long long cnt = counts[l - 1];
        if (cnt == 0) continue;
        SUM += cnt * (combinations(l - 1, k - 2) + combinations(n - l - 1, k - 2));
    }

    SDk -= SUM / 4;
    return SDk;
}

// [[Rcpp::export]]
long long SDk(NumericMatrix Xinput_mat, NumericVector x_vec, int k) {
    int n = Xinput_mat.nrow();
    if (Xinput_mat.ncol() != 3 || x_vec.size() != 3) {
        stop("Xinput must be n x 3 and x must be length 3");
    }

    // Convert matrix to std::vector<Point3D>
    std::vector<Point3D> Xinput(n);
    for (int i = 0; i < n; i++) {
        Xinput[i].x = Xinput_mat(i, 0);
        Xinput[i].y = Xinput_mat(i, 1);
        Xinput[i].z = Xinput_mat(i, 2);
    }

    // Convert vector to Point3D
    Point3D x = { x_vec[0], x_vec[1], x_vec[2] };

    return SDk(Xinput, x, k);
}

long long SDk_parallel(const std::vector<Point3D>& Xinput, const Point3D& x, int k, int threads = 0) {
    int n = Xinput.size();
    if (n < k) return 0LL;
    std::vector<Point3D> X;
    X.reserve(n);
    for (const auto& p : Xinput) {
        Point3D q{p.x - x.x, p.y - x.y, p.z - x.z};
        if (fabs(q.y) < EPS) q.y = EPS; // For sorting in the projection plane we must divide by q.y
        if (fabs(q.x) > EPS) X.push_back(q);
    }
    n = (int)X.size();
    if (n == 0 || k > n) return 0LL;


    struct HelperStruct {
        double x, y;
        int color, label;
        double sort_key;
    };

    struct PointY {
        double y1, y2;
        int sign;
    };
    std::vector<PointY> Y; Y.reserve(n);
    int n_up = 0;
    for (const auto& p : X) {
        double y1 = p.x / p.z;
        double y2 = p.y / p.z;
        int s = sgn(p.z);
        Y.push_back({y1, y2, s});
        if (s == 1) n_up++;
    }
    int n_down = n - n_up;

    long long SDk = combinations(n, k) - combinations(n_up, k) - combinations(n_down, k);
    if (n == 1) return SDk;

    struct SDkReducer : public RcppParallel::Worker {
        const std::vector<PointY>& Y;
        int n;
        std::vector<long long> counts;
        std::vector<HelperStruct> sorted_points;

        SDkReducer(const std::vector<PointY>& Y_, int n_)
            : Y(Y_), n(n_), counts(static_cast<std::size_t>(n_ - 1), 0LL),
              sorted_points(static_cast<std::size_t>(n_ - 1)) {}

        SDkReducer(SDkReducer& other, RcppParallel::Split)
            : Y(other.Y), n(other.n), counts(static_cast<std::size_t>(other.n - 1), 0LL),
              sorted_points(static_cast<std::size_t>(other.n - 1)) {}

        void operator()(std::size_t begin, std::size_t end) {
            for (std::size_t idx_ = begin; idx_ < end; ++idx_) {
                const int idx = static_cast<int>(idx_);
                int j = 0;
                for (int i = 0; i < n; i++) {
                    if (i == idx) continue;
                    HelperStruct* sorted_point = &sorted_points[j];
                    const PointY* pointy_i = &Y[i];
                    const PointY* pointy_idx = &Y[idx];
                    sorted_point->x = pointy_i->y1 - pointy_idx->y1;
                    sorted_point->y = pointy_i->y2 - pointy_idx->y2;
                    sorted_point->color = (pointy_i->sign != pointy_idx->sign) ? 1 : 0;
                    sorted_point->label = (sorted_point->y * pointy_i->sign > 0.0) ? 1 : 0;
                    sorted_point->sort_key = -(sorted_point->x / sorted_point->y);
                    j++;
                }

                const int m = n - 1;
                if (m == 0) continue;

                std::sort(sorted_points.begin(), sorted_points.end(),
                          [](const auto& a, const auto& b) {
                              return a.sort_key < b.sort_key;
                          });

                int EL = 0;
                for (int jj = 1; jj < m; jj++) EL += sorted_points[jj].label;
                counts[static_cast<std::size_t>(EL)] += sorted_points[0].color;

                for (int jj = 0; jj < m - 1; jj++) {
                    if (sorted_points[jj].label == 0) EL++;
                    if (sorted_points[jj + 1].label == 1) EL--;
                    counts[static_cast<std::size_t>(EL)] += sorted_points[jj + 1].color;
                }
            }
        }

        void join(SDkReducer& rhs) {
            for (std::size_t i = 0; i < counts.size(); ++i) {
                counts[i] += rhs.counts[i];
            }
        }
    };

    SDkReducer reducer(Y, n);
    simplicialdepth::parallel_reduce(0, static_cast<std::size_t>(n), reducer, 1, threads);

    const std::vector<long long>& counts = reducer.counts;
    long long SUM = 0LL;
    for (int l = 1; l <= n - 1; l++) {
        long long cnt = counts[l - 1];
        if (cnt == 0) continue;
        SUM += cnt * (combinations(l - 1, k - 2) + combinations(n - l - 1, k - 2));
    }

    SDk -= SUM / 4;
    // std::cout << "SDk: " << SDk << std::endl;
    return SDk;
}

// [[Rcpp::export]]
long long SDk_parallel(NumericMatrix Xinput_mat, NumericVector x_vec, int k) {
    int n = Xinput_mat.nrow();
    if (n < k) return 0LL;

    if (Xinput_mat.ncol() != 3 || x_vec.size() != 3) {
        stop("Xinput must be n x 3 and x must be length of 3");
    }

    // Convert matrix to std::vector<Point3D>
    std::vector<Point3D> Xinput(n);
    // for (int i = 0; i < n; i++)
    //     Xinput.push_back({Xinput_mat(i,0), Xinput_mat(i,1), Xinput_mat(i,2)});
    for (int i = 0; i < n; i++) {
        Xinput[i].x = Xinput_mat(i, 0);
        Xinput[i].y = Xinput_mat(i, 1);
        Xinput[i].z = Xinput_mat(i, 2);
    }

    // Convert vector to Point3D
    Point3D x = { x_vec[0], x_vec[1], x_vec[2] };

    return SDk_parallel(Xinput, x, k, 0);
}

// Like SDk_parallel(), but allows per-call thread control.
// [[Rcpp::export]]
long long SDk_parallel_threads(NumericMatrix Xinput_mat, NumericVector x_vec, int k, int threads = 0) {
    int n = Xinput_mat.nrow();
    if (n < k) return 0LL;

    if (Xinput_mat.ncol() != 3 || x_vec.size() != 3) {
        stop("Xinput must be n x 3 and x must be length of 3");
    }

    std::vector<Point3D> Xinput(static_cast<std::size_t>(n));
    for (int i = 0; i < n; i++) {
        Xinput[static_cast<std::size_t>(i)].x = Xinput_mat(i, 0);
        Xinput[static_cast<std::size_t>(i)].y = Xinput_mat(i, 1);
        Xinput[static_cast<std::size_t>(i)].z = Xinput_mat(i, 2);
    }

    Point3D x = { x_vec[0], x_vec[1], x_vec[2] };
    return SDk_parallel(Xinput, x, k, threads);
}

// Compute leave-one-out SDk for every row of X (each row is x; X without that row is sample).
// This is the fast path for x = NULL wrappers.
// [[Rcpp::export]]
NumericVector SDk_all_points(NumericMatrix Xinput_mat, int k, int threads = 0) {
    const int n = Xinput_mat.nrow();
    if (n == 0) return NumericVector(0);
    if (Xinput_mat.ncol() != 3) stop("Xinput must be n x 3");

    std::vector<Point3D> Xfull(static_cast<std::size_t>(n));
    for (int i = 0; i < n; i++) {
        Xfull[static_cast<std::size_t>(i)].x = Xinput_mat(i, 0);
        Xfull[static_cast<std::size_t>(i)].y = Xinput_mat(i, 1);
        Xfull[static_cast<std::size_t>(i)].z = Xinput_mat(i, 2);
    }

    NumericVector out(n);
    RcppParallel::RVector<double> out_view(out);

    struct Worker : public RcppParallel::Worker {
        const std::vector<Point3D>& X;
        const int n;
        const int k;
        RcppParallel::RVector<double> out;

        Worker(const std::vector<Point3D>& X_, int n_, int k_, RcppParallel::RVector<double> out_)
            : X(X_), n(n_), k(k_), out(out_) {}

        void operator()(std::size_t begin, std::size_t end) {
            std::vector<Point3D> sample;
            sample.reserve(static_cast<std::size_t>(n - 1));

            for (std::size_t ii = begin; ii < end; ++ii) {
                const int i = static_cast<int>(ii);
                sample.clear();
                for (int j = 0; j < n; ++j) {
                    if (j == i) continue;
                    sample.push_back(X[static_cast<std::size_t>(j)]);
                }

                const Point3D x = X[static_cast<std::size_t>(i)];
                out[ii] = static_cast<double>(SDk(sample, x, k));
            }
        }
    };

    Worker worker(Xfull, n, k, out_view);
    simplicialdepth::parallel_for(0, static_cast<std::size_t>(n), worker, threads);
    return out;
}



