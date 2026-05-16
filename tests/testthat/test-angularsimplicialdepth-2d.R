triangle_contains_origin <- function(p1, p2, p3) {
    cross2 <- function(a, b) a[1] * b[2] - a[2] * b[1]
    orient <- function(a, b, c) cross2(b - a, c - a)
    eps <- 1e-12

    area2 <- orient(p1, p2, p3)
    if (abs(area2) <= eps) {
        return(FALSE)
    }

    o1 <- orient(p1, p2, c(0, 0))
    o2 <- orient(p2, p3, c(0, 0))
    o3 <- orient(p3, p1, c(0, 0))

    has_pos <- (o1 > eps) || (o2 > eps) || (o3 > eps)
    has_neg <- (o1 < -eps) || (o2 < -eps) || (o3 < -eps)

    !(has_pos && has_neg)
}

testthat::test_that("2D angularsimplicialdepth for single point", {
    set.seed(456)
    n <- 30
    X <- matrix(rnorm(2 * n), ncol = 2)
    X <- X / sqrt(rowSums(X^2))
    x <- rnorm(2)

    res <- simplicialdepth::angularsimplicialdepth(X, x)

    expected <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            if (triangle_contains_origin(X[i, ], X[j, ], -x)) {
                expected <- expected + 1
            }
        }
    }
    expected <- expected / choose(n, 2)
    testthat::expect_equal(as.numeric(res$depth), expected)
})

testthat::test_that("2D all-points ASD matches per-point circular_asd", {
    set.seed(456)
    n <- 50
    X <- matrix(rnorm(2 * n), ncol = 2)
    X <- X / sqrt(rowSums(X^2))

    circular_asd_single <- get("circular_asd", envir = asNamespace("simplicialdepth"))
    circular_asd_all <- get("circular_asd_all_points", envir = asNamespace("simplicialdepth"))

    res_counts <- circular_asd_all(X)$depth
    expected_counts <- vapply(seq_len(n), function(k) {
        circular_asd_single(X[-k, , drop = FALSE], X[k, ])
    }, numeric(1))

    testthat::expect_equal(as.numeric(res_counts), as.numeric(expected_counts))
})
