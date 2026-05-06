point_in_simplex_2d <- function(q, a, b, c, eps = 1e-10) {
    m <- rbind(
        c(a[1], b[1], c[1]),
        c(a[2], b[2], c[2]),
        c(1, 1, 1)
    )
    if (abs(det(m)) <= eps) return(FALSE)
    lambdas <- solve(m, c(q[1], q[2], 1))
    all(lambdas >= -eps)
}

point_in_simplex_3d <- function(q, a, b, c, d, eps = 1e-10) {
    m <- rbind(
        c(a[1], b[1], c[1], d[1]),
        c(a[2], b[2], c[2], d[2]),
        c(a[3], b[3], c[3], d[3]),
        c(1, 1, 1, 1)
    )
    if (abs(det(m)) <= eps) return(FALSE)
    lambdas <- solve(m, c(q[1], q[2], q[3], 1))
    all(lambdas >= -eps)
}

bruteforce_sd2d_count <- function(X, q) {
    n <- nrow(X)
    if (n < 3) return(0)
    total <- 0L
    for (i in 1:(n - 2)) {
        for (j in (i + 1):(n - 1)) {
            for (k in (j + 1):n) {
                if (point_in_simplex_2d(q, X[i, ], X[j, ], X[k, ])) {
                    total <- total + 1L
                }
            }
        }
    }
    total
}

bruteforce_sd3d_count <- function(X, q) {
    n <- nrow(X)
    if (n < 4) return(0)
    total <- 0L
    for (i in 1:(n - 3)) {
        for (j in (i + 1):(n - 2)) {
            for (k in (j + 1):(n - 1)) {
                for (l in (k + 1):n) {
                    if (point_in_simplex_3d(q, X[i, ], X[j, ], X[k, ], X[l, ])) {
                        total <- total + 1L
                    }
                }
            }
        }
    }
    total
}

testthat::test_that("2D simplicialdepth single query matches brute force", {
    set.seed(101)
    X <- matrix(rnorm(2 * 14), ncol = 2)
    q <- c(0.2, -0.1)

    res <- simplicialdepth::simplicialdepth(X, q)
    expected <- bruteforce_sd2d_count(X, q) / choose(nrow(X), 3)

    testthat::expect_equal(as.numeric(res$depth), as.numeric(expected))
})

testthat::test_that("2D simplicialdepth all points matches brute force leave-one-out", {
    set.seed(102)
    X <- matrix(rnorm(2 * 12), ncol = 2)

    res <- simplicialdepth::simplicialdepth(X)
    expected <- vapply(seq_len(nrow(X)), function(i) {
        bruteforce_sd2d_count(X[-i, , drop = FALSE], X[i, ]) / choose(nrow(X) - 1, 3)
    }, numeric(1))

    testthat::expect_equal(as.numeric(res$depth), expected)
})

testthat::test_that("3D simplicialdepth single query matches brute force", {
    set.seed(201)
    X <- matrix(rnorm(3 * 11), ncol = 3)
    q <- c(0.1, -0.2, 0.3)

    res <- simplicialdepth::simplicialdepth(X, q, threads = 1L)
    expected <- bruteforce_sd3d_count(X, q) / choose(nrow(X), 4)

    testthat::expect_equal(as.numeric(res$depth), as.numeric(expected))
})

testthat::test_that("3D simplicialdepth all points matches brute force leave-one-out", {
    set.seed(202)
    X <- matrix(rnorm(3 * 10), ncol = 3)

    res <- simplicialdepth::simplicialdepth(X, threads = 1L)
    expected <- vapply(seq_len(nrow(X)), function(i) {
        bruteforce_sd3d_count(X[-i, , drop = FALSE], X[i, ]) / choose(nrow(X) - 1, 4)
    }, numeric(1))

    testthat::expect_equal(as.numeric(res$depth), expected)
})
