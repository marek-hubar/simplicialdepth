if (!requireNamespace("testthat", quietly = TRUE)) {
    stop("Please install testthat to run qa-tests.")
}

testthat::test_dir("qa-tests", reporter = "summary", stop_on_failure = TRUE)
