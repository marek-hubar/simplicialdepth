if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE, export_all = FALSE)
} else if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE, export_all = FALSE)
} else if (requireNamespace("simplicialdepth", quietly = TRUE)) {
    library(simplicialdepth)
} else {
    stop(
        "Cannot load simplicialdepth for qa-tests. ",
        "Install pkgload/devtools, or install the package first."
    )
}
