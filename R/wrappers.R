input_check <- function(X, x=NULL) {
    if (!is.matrix(X) || !is.numeric(X))
        stop("X must be a numeric matrix")

    n <- nrow(X)
    d <- ncol(X)

    if (!is.null(x)) {
        if ((!is.vector(x) && !is.matrix(x)) || !is.numeric(x))
            stop("x must be a numeric matrix or a numeric vector.")
        if (is.matrix(x)) {
            if (ncol(x) != d) {
                stop("x must have the same number of columns as X (i.e. ncol(x) = ncol(X)).")
            }
        } else {
            if (length(x) != d) {
                stop("x must have length equal to the number of columns of X (i.e. length(x) = ncol(X)).")
            }
        }
    }
}

remove_x_from_X <- function(X, x) {
    keep <- rowSums(X == matrix(x, nrow(X), ncol(X), byrow = TRUE)) != ncol(X)
    return(X[keep, ])
}

#' Constructor for depth_result objects
#'
#' @param depth Numeric vector of depths.
#' @param max_depth Maximum depth (scalar).
#' @param max_point Coordinates of the point with maximum depth.
#' @param max_index Index of the point with maximum depth.
#'
#' @return An object of class \code{depth_result}.
#' @export
depth_result <- function(depth, max_depth, max_point, max_index) {
  structure(
    list(
      depth     = depth,
      max_depth = max_depth,
      max_point = max_point,
      max_index = max_index
    ),
    class = "depth_result"
  )
}

#' @method print depth_result
#' @export Print method
print.depth_result <- function(x, ...) {
  cat("Depth result\n")
  cat("  Number of points:", length(x$depth), "\n")
  cat("  Maximum depth   :", x$max_depth, "\n")
  cat("  At point        :", paste(x$max_point, collapse = ", "), "\n")
  cat("  With index      :", paste(x$max_index, collapse = ", "), "\n")
  invisible(x)
}

spherical_asd_on_one_point <- function(X, i) {
    new_X <- X[-i,]
    x <- X[i,]
    return(spherical_asd(new_X, x))
}

#' Angular Simplicial Depth
#'
#' Computes the angular simplicial depth for 2D or 3D data, either for all
#' points in a sample or for user-specified query points.
#'
#' @param X A numeric matrix of size \eqn{n \times d} containing the sample points, 
#'   where each row is an observation and each column is a coordinate. 
#'   All points are projected onto a unit circle (2D) or sphere (3D) before depth computation.
#' @param x Optional. A query point or set of points:
#'   \itemize{
#'     \item If `NULL` (default): computes depths for all rows of `X`.
#'     \item If a numeric vector of length \eqn{d}: computes the depth of this single point.
#'     \item If a numeric matrix with \eqn{d} columns: computes depths for each row.
#'   }
#'
#' @details
#' Depths are normalized by appropriate binomial coefficients:
#' \itemize{
#'   \item In 2D: divide by \eqn{choose(n, 2)} for query points, \eqn{choose(n-1, 2)} for sample points.
#'   \item In 3D: divide by \eqn{choose(n, 3)} or \eqn{choose(n-1, 3)}.
#' }
#' The function identifies the point with maximum depth.
#'
#' @return An object of class \code{depth_result} with components:
#' \describe{
#'   \item{depth}{Numeric vector of depths for each query point or sample point.}
#'   \item{max_depth}{Maximum depth observed.}
#'   \item{max_point}{Coordinates of the point with maximum depth.}
#'   \item{max_index}{Index of the point with maximum depth.}
#' }
#'
#' @examples
#' X <- matrix(rnorm(20), ncol = 2)
#' result <- angularsimplicialdepth(X)
#' result$depth       # depths of each point
#' result$max_depth   # maximum depth
#' result$max_point   # point with maximum depth
#' result$max_index   # index of max depth
#'
#' @export
angularsimplicialdepth <- function(X, x = NULL, threads = NULL) {
    input_check(X, x)

    threads_i <- if (is.null(threads)) 0L else as.integer(threads)

    n <- nrow(X)
    d <- ncol(X)
    if (is.null(x)) {
        if (n < d+1)
            stop("For x=NULL, X must have more rows than columns.")
    } else {
        if (n < d)
            stop("X must have at least as many rows as columns.")
    }

    if (d == 2) {
        if (is.null(x)) {
            result <- circular_asd_all_arcs(X)$depth
            result <- result / choose(n - 1, 2)
            max_depth_index <- which.max(result)
            max_depth <- result[max_depth_index]
            return(depth_result(depth = result,
                                max_depth = max_depth,
                                max_point = X[max_depth_index,],
                                max_index = max_depth_index))
        } else {
            if (is.matrix(x)) {
                result <- apply(x, 1, function(row) {
                    X_new <- remove_x_from_X(X, row)
                    if (nrow(X_new) >= 2) {
                        circular_asd(X_new, row) / choose(nrow(X_new), 2)
                    } else {
                        return(0)
                    }})
                max_depth_index <- which.max(result)
                return(depth_result(depth=result,
                                    max_depth=result[max_depth_index],
                                    max_point=x[max_depth_index, ],
                                    max_index=max_depth_index))
            }
            # The only option left is that x is a numeric vector of proper length
            X_new <- remove_x_from_X(X, x)
            if (nrow(X_new) >= 2) {
                result <- circular_asd(X_new, x) / choose(nrow(X_new), 2)
            } else {
                result <- 0
            }
            return(depth_result(depth=c(result),
                                max_depth=result,
                                max_point=x,
                                max_index=1))
        }
    }
    if (d == 3) {
        if (is.null(x)) {
            result <- spherical_asd_all_points(X, threads_i) / choose(n-1, 3)
            max_depth_index <- which.max(result)
            return(depth_result(depth=result,
                                max_depth=result[max_depth_index],
                                max_point=X[max_depth_index, ],
                                max_index=max_depth_index))
        }
        if (is.matrix(x)) {
            # result <- apply(x, 1, function(row) spherical_asd(X, row)) / choose(n, 3)
            result <- apply(x, 1, function(row) {
                X_new <- remove_x_from_X(X, row)
                if (nrow(X_new) >= 3) {
                    spherical_asd(X_new, -row) / choose(nrow(X_new), 3)
                } else {
                    return(0)
                }})
            max_depth_index <- which.max(result)
            return(depth_result(depth=result,
                                max_depth=result[max_depth_index],
                                max_point=x[max_depth_index, ],
                                max_index=max_depth_index))
        }


        X_new <- remove_x_from_X(X, x)

        if (nrow(X_new) >= 3) {
            result <- spherical_asd(X_new, -x) / choose(nrow(X_new), 3)
        } else {
            result <- 0
        }
        return(depth_result(depth=c(result),
                            max_depth=result,
                            max_point=x,
                            max_index=1))
    }
    stop("Not implemented: only 2D or 3D simplicial depth is supported")
}

sd2d_on_one_point <- function(X, i) {
    new_X <- X[-i,]
    x <- X[i,]
    return(simplicial_depth_2d(new_X, x))
}

sdk_on_one_point <- function(X, i, k) {
    new_X <- X[-i,]
    x <- X[i,]
    return(SDk_parallel(new_X, x, k))
}

#' Simplicial Depth
#'
#' Computes the simplicial depth for 2D or 3D data, either for all points in a sample
#' or for user-specified query points.
#'
#' @param X Numeric matrix of size \eqn{n \times d} containing sample points.
#' @param x Optional. A query point or set of points (vector or matrix) with same column dimension as `X`.
#'
#' @details
#' Depths are normalized by:
#' \itemize{
#'   \item 2D: \eqn{choose(n-1, 3)} for sample points, \eqn{choose(n, 3)} for query points.
#'   \item 3D: \eqn{choose(n-1, 4)} for sample points, \eqn{choose(n, 4)} for query points.
#' }
#'
#' @return An object of class \code{depth_result} with components:
#' \describe{
#'   \item{depth}{Numeric vector of depths.}
#'   \item{max_depth}{Maximum depth.}
#'   \item{max_point}{Coordinates of the point with maximum depth.}
#'   \item{max_index}{Index of the point with maximum depth.}
#' }
#'
#' @examples
#' X <- matrix(rnorm(20), ncol = 2)
#' result <- simplicialdepth(X)
#' result$depth
#' result$max_depth
#' result$max_point
#' result$max_index
#'
#' @export
simplicialdepth <- function(X, x=NULL, threads = NULL) {
    input_check(X, x)

    threads_i <- if (is.null(threads)) 0L else as.integer(threads)

    # X <- remove_x_from_X(X, x)

    n <- nrow(X)
    d <- ncol(X)

    if (is.null(x)) {
        if (n < d+2)
            stop("For x=NULL, X must have at least 2 more rows than columns.")
    } else {
        if (n < d+1)
            stop("X must have more rows than columns.")
    }

    if (d == 2) {
        if (is.null(x)) {
            result <- vapply(seq_len(nrow(X)), function(i) sd2d_on_one_point(X, i), numeric(1)) / choose(n-1, 3)
            max_depth_index <- which.max(result)
            return(depth_result(depth=result,
                                max_depth=result[max_depth_index],
                                max_point=X[max_depth_index, ],
                                max_index=max_depth_index))
        }
        if (is.matrix(x)) {
            result <- apply(x, 1, function(row) {
                X_new <- remove_x_from_X(X, row)
                if (nrow(X_new) >= 3) {
                    simplicial_depth_2d(X_new, row) / choose(nrow(X_new), 3)
                } else {
                    return(0)
                }})
            max_depth_index <- which.max(result)
            return(depth_result(depth=result,
                                max_depth=result[max_depth_index],
                                max_point=x[max_depth_index, ],
                                max_index=max_depth_index))
        }
        # The only remaining option is that x is a numeric vector
        X_new <- remove_x_from_X(X, x)
        if (nrow(X_new) >= 3) {
            result <- simplicial_depth_2d(X,x) / choose(n, 3)
        } else {
            result <- 0
        }
        return(depth_result(depth=c(result),
                            max_depth=result,
                            max_point=x,
                            max_index=1))
    }
    if (d == 3) {
        if (is.null(x)) {
            result <- SDk_all_points(X, 4L, threads_i) / choose(n-1, 4)
            max_depth_index <- which.max(result)
            return(depth_result(depth=result,
                                max_depth=result[max_depth_index],
                                max_point=X[max_depth_index, ],
                                max_index=max_depth_index))
        }
        if (is.matrix(x)) {
            result <- apply(x, 1, function(row) {
                X_new <- remove_x_from_X(X, row)
                if (nrow(X_new) >= 4) {
                    SDk_parallel_threads(X_new, row, 4L, threads_i) / choose(nrow(X_new), 4)
                } else {
                    return(0)
                }})
            max_depth_index <- which.max(result)
            return(depth_result(depth=result,
                                max_depth=result[max_depth_index],
                                max_point=x[max_depth_index, ],
                                max_index=max_depth_index))
        }
        X_new <- remove_x_from_X(X, x)
        if (nrow(X_new) >= 4) {
            result <- SDk_parallel_threads(X_new, x, 4L, threads_i) / choose(nrow(X_new), 4)
        } else {
            return(0)
        }
        return(depth_result(depth=c(result),
                            max_depth=result,
                            max_point=x,
                            max_index=1))
    }
}

#' k-Hull Depth
#'
#' Computes the \eqn{k}-hull depth for 2D or 3D data, either for all points in a sample
#' or for user-specified query points.
#'
#' @param X Numeric matrix of size \eqn{n \times d} containing sample points.
#' @param x Optional. Query point(s) as a vector or matrix with same number of columns as `X`.
#' @param k Integer, hull size parameter. Must satisfy \eqn{k > d} and \eqn{k \le n}.
#'
#' @details
#' Depths are normalized by \eqn{choose(n-1, k)} for sample points and
#' \eqn{choose(n, k)} for query points. If \eqn{k \le d} or \eqn{k > n}, zero depths are returned.
#'
#' @return An object of class \code{depth_result} with components:
#' \describe{
#'   \item{depth}{Numeric vector of depths.}
#'   \item{max_depth}{Maximum depth.}
#'   \item{max_point}{Coordinates of the point with maximum depth.}
#'   \item{max_index}{Index of the point with maximum depth.}
#' }
#'
#' @examples
#' X <- matrix(rnorm(20), ncol = 2)
#' result <- khulldepth(X, k = 4)
#' result$depth
#' result$max_depth
#' result$max_point
#' result$max_index
#'
#' @export
khulldepth <- function(X, x=NULL, k, threads = NULL) {
    input_check(X, x)

    threads_i <- if (is.null(threads)) 0L else as.integer(threads)

    n <- nrow(X)
    d <- ncol(X)

    if (d != 3) {
        stop("Not implemented: khulldepth currently supports only 3D data (ncol(X) = 3).")
    }
    
    if (is.null(x)) {
        if (n < k+1)
            stop("For x=NULL, X must have more than k rows.")
    } else {
        if (n < k)
            stop("X must have at least k rows.")
    }
    if (k <= d)
        stop("k must be bigger than the number of columns of X.")

    if (is.null(x)) {
        result <- SDk_all_points(X, as.integer(k), threads_i) / choose(n-1, k)
        max_depth_index <- which.max(result)
        return(depth_result(depth=result,
                            max_depth=result[max_depth_index],
                            max_point=X[max_depth_index, ],
                            max_index=max_depth_index))
    }
    if (is.matrix(x)) {
        result <- apply(x, 1, function(row) {
            X_new <- remove_x_from_X(X, row)
            if (nrow(X_new) >= k) {
                SDk_parallel_threads(X_new, row, as.integer(k), threads_i) / choose(nrow(X_new), k)
            } else {
                return(0)
            }})
        max_depth_index <- which.max(result)
        return(depth_result(depth=result,
                            max_depth=result[max_depth_index],
                            max_point=x[max_depth_index, ],
                            max_index=max_depth_index))
    }
    # The only option left is that x is a numeric vector of length equal to d
    X_new <- remove_x_from_X(X, x)
    if (nrow(X_new) >= k) {
        result <- SDk_parallel_threads(X_new, x, as.integer(k), threads_i) / choose(nrow(X_new), k)
    } else {
        result <- 0
    }
    return(depth_result(depth=c(result),
                        max_depth=result,
                        max_point=x,
                        max_index=1))
}
