testthat::test_that("3D angularsimplicialdepth uses spherical_asd normalization", {
    set.seed(456)
    X <- matrix(rnorm(60), ncol = 3)
    X <- X / sqrt(rowSums(X^2))
    spherical_asd_all_points_cpp <- get("spherical_asd_all_points", envir = asNamespace("simplicialdepth"))

    res_all <- simplicialdepth::angularsimplicialdepth(X, threads = 1L)
    raw_all <- as.numeric(spherical_asd_all_points_cpp(X, threads = 1L))
    expected_all <- raw_all / choose(nrow(X) - 1, 3)

    testthat::expect_equal(as.numeric(res_all$depth), expected_all)
})

testthat::test_that("3D angularsimplicialdepth single query matches spherical_asd", {
    set.seed(789)
    X <- matrix(rnorm(75), ncol = 3)
    X <- X / sqrt(rowSums(X^2))
    q <- X[3, ]
    spherical_asd_cpp <- get("spherical_asd", envir = asNamespace("simplicialdepth"))

    res_q <- simplicialdepth::angularsimplicialdepth(X, q, threads = 1L)
    X_new <- X[-3, , drop = FALSE]
    raw_q <- as.numeric(spherical_asd_cpp(X_new, -q))
    expected_q <- raw_q / choose(nrow(X_new), 3)

    testthat::expect_equal(as.numeric(res_q$depth), expected_q)
})

testthat::test_that("3D angularsimplicialdepth matrix queries match vector queries", {
    set.seed(999)
    X <- matrix(rnorm(90), ncol = 3)
    X <- X / sqrt(rowSums(X^2))
    Q <- X[c(2, 7, 15), , drop = FALSE]

    res_mat <- simplicialdepth::angularsimplicialdepth(X, Q, threads = 1L)
    res_vec <- vapply(
        seq_len(nrow(Q)),
        function(i) as.numeric(simplicialdepth::angularsimplicialdepth(X, Q[i, ], threads = 1L)$depth),
        numeric(1)
    )

    testthat::expect_equal(as.numeric(res_mat$depth), res_vec)
})
