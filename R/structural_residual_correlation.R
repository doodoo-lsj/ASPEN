# Correlated structural errors preserving the fitted conditional outcome variance.
structural_rho_version <- function() "structural_rho_conditional_variance_v1"

structural_rho_adjustment <- function(Sigma_M, sigma_y0, rho, tolerance = 1e-10) {
  Sigma_M <- as.matrix(Sigma_M)
  k <- nrow(Sigma_M)
  rho <- if (is.null(dim(rho))) matrix(rho, nrow = 1L) else as.matrix(rho)
  stopifnot(ncol(Sigma_M) == k, ncol(rho) == k,
            length(sigma_y0) == 1L, is.finite(tolerance), tolerance > 0, tolerance < 1)
  invalid <- matrix(NA_real_, nrow(rho), k, dimnames = list(NULL, colnames(Sigma_M)))
  result <- list(adjustment = invalid, covariance = invalid,
                 sigma_y = rep(NA_real_, nrow(rho)),
                 q = rep(NA_real_, nrow(rho)), valid = rep(FALSE, nrow(rho)))
  if (!is.finite(sigma_y0) || sigma_y0 <= 0 || any(!is.finite(Sigma_M)) ||
      !isTRUE(isSymmetric(Sigma_M)) ||
      inherits(try(chol(Sigma_M), silent = TRUE), "try-error")) return(result)
  sd_m <- sqrt(diag(Sigma_M))
  R <- Sigma_M / outer(sd_m, sd_m)
  solved <- tryCatch(t(solve(R, t(rho))), error = function(e) invalid)
  q <- rowSums(rho * solved)
  valid <- is.finite(q) & q >= -tolerance & q < 1 - tolerance &
    apply(is.finite(rho) & abs(rho) < 1, 1L, all) &
    apply(is.finite(solved), 1L, all)
  sigma_y <- rep(NA_real_, nrow(rho))
  sigma_y[valid] <- sigma_y0 / sqrt(1 - pmax(0, q[valid]))
  result$adjustment <- sweep(solved * sigma_y, 2L, sd_m, "/")
  result$covariance <- sweep(rho * sigma_y, 2L, sd_m, "*")
  result$sigma_y <- sigma_y
  result$q <- q
  result$valid <- valid & is.finite(sigma_y) &
    apply(is.finite(result$adjustment), 1L, all) &
    apply(is.finite(result$covariance), 1L, all)
  result$adjustment[!result$valid, ] <- NA_real_
  result$covariance[!result$valid, ] <- NA_real_
  result$sigma_y[!result$valid] <- NA_real_
  result
}
