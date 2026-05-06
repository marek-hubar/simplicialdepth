testthat::test_that("spherical_asd_all_points matches leave-one-out spherical_asd", {
    set.seed(123)
    X <- matrix(rnorm(45), ncol = 3)
    X <- X / sqrt(rowSums(X^2))

    spherical_asd_all_points_cpp <- get("spherical_asd_all_points", envir = asNamespace("simplicialdepth"))
    spherical_asd_cpp <- get("spherical_asd", envir = asNamespace("simplicialdepth"))

    all_counts <- as.numeric(spherical_asd_all_points_cpp(X, threads = 1L))
    loo_counts <- vapply(
        seq_len(nrow(X)),
        function(i) as.numeric(spherical_asd_cpp(X[-i, , drop = FALSE], X[i, ])),
        numeric(1)
    )

    testthat::expect_equal(all_counts, loo_counts)
})
