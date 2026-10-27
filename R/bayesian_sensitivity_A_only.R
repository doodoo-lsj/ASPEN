# Bayesian Sensitivity A runtime helpers for the Shiny app.
#
# This file intentionally contains only residual-correlation Sensitivity A.
# It does not contain B1/B2 latent-U Stan code and does not fit Bayesian models.

source(local({
  source_paths <- Filter(function(x) is.character(x) && length(x) == 1L,
                         lapply(sys.frames(), function(frame) frame$ofile))
  file.path(dirname(tail(source_paths, 1L)[[1L]]), "structural_residual_correlation.R")
}), local = TRUE)

bayes_A_bundle_mediators <- function(residual_draws) {
  residual_draws$mediators %||% c("PC1_R", "PC2_R", "PC3")
}

bayes_A_bundle_exposures <- function(residual_draws) {
  residual_draws$exposures %||% names(residual_draws$alpha) %||% c("X_adjacent", "X_far")
}

bayes_A_contrast_label <- function(exposure, residual_draws = NULL) {
  labels <- residual_draws$contrast_labels %||% NULL
  if (!is.null(labels)) {
    out <- unname(labels[exposure])
    out <- ifelse(is.na(out), exposure, out)
    return(out)
  }

  dplyr::case_when(
    exposure == "X_adjacent" ~ "adjacent vs inside",
    exposure == "X_far" ~ "far vs inside",
    TRUE ~ exposure
  )
}

bayes_A_residual_correlation_array <- function(residual_draws, mediators = bayes_A_bundle_mediators(residual_draws)) {
  n_draws <- residual_draws$n_draws
  k <- length(mediators)

  if (!is.null(residual_draws$R)) {
    R <- residual_draws$R
    if (length(dim(R)) != 3 || dim(R)[[1]] != n_draws || dim(R)[[2]] != k || dim(R)[[3]] != k) {
      stop("residual_draws$R must be an n_draws x K x K residual-correlation array.", call. = FALSE)
    }
    return(R)
  }

  required_old <- c("r12", "r13", "r23")
  missing_old <- setdiff(required_old, names(residual_draws))
  if (k != 3 || length(missing_old) > 0) {
    stop("Residual correlation draws must be supplied as R array for dynamic mediator sets.", call. = FALSE)
  }

  R <- array(0, dim = c(n_draws, k, k), dimnames = list(NULL, mediators, mediators))
  for (s in seq_len(n_draws)) {
    R[s, , ] <- diag(1, k)
  }
  R[, 1, 2] <- R[, 2, 1] <- residual_draws$r12
  R[, 1, 3] <- R[, 3, 1] <- residual_draws$r13
  R[, 2, 3] <- R[, 3, 2] <- residual_draws$r23
  R
}

validate_bayesian_A_residual_draws <- function(residual_draws) {
  required <- c("n_draws", "sig", "sigma_y", "beta", "alpha", "total_effect", "outcome_direct_coef")
  missing <- setdiff(required, names(residual_draws))
  if (length(missing) > 0) {
    stop("Bayesian Sensitivity A residual draw bundle is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  mediators <- bayes_A_bundle_mediators(residual_draws)
  exposures <- bayes_A_bundle_exposures(residual_draws)

  if (length(mediators) == 0 || length(exposures) == 0) {
    stop("Bayesian Sensitivity A residual draw bundle must define at least one mediator and exposure.", call. = FALSE)
  }
  if (!all(mediators %in% colnames(residual_draws$sig))) {
    stop("residual_draws$sig must contain mediator columns: ", paste(mediators, collapse = ", "), call. = FALSE)
  }
  if (!all(mediators %in% colnames(residual_draws$beta))) {
    stop("residual_draws$beta must contain mediator columns: ", paste(mediators, collapse = ", "), call. = FALSE)
  }
  if (!all(exposures %in% names(residual_draws$alpha))) {
    stop("residual_draws$alpha must contain exposure entries: ", paste(exposures, collapse = ", "), call. = FALSE)
  }
  if (!all(exposures %in% colnames(residual_draws$total_effect))) {
    stop("residual_draws$total_effect must contain exposure columns: ", paste(exposures, collapse = ", "), call. = FALSE)
  }
  if (!all(exposures %in% colnames(residual_draws$outcome_direct_coef))) {
    stop("residual_draws$outcome_direct_coef must contain exposure columns: ", paste(exposures, collapse = ", "), call. = FALSE)
  }

  n_draws <- residual_draws$n_draws
  lengths <- c(
    nrow(residual_draws$sig),
    length(residual_draws$sigma_y),
    nrow(residual_draws$beta),
    nrow(residual_draws$total_effect),
    nrow(residual_draws$outcome_direct_coef),
    vapply(exposures, function(exposure) nrow(residual_draws$alpha[[exposure]]), integer(1))
  )
  if (!all(lengths == n_draws)) {
    stop("Inconsistent posterior draw counts in Bayesian Sensitivity A residual draw bundle.", call. = FALSE)
  }

  invisible(bayes_A_residual_correlation_array(residual_draws, mediators))
}

bayesian_A_rho_grid <- function(rho_grid = seq(-0.5, 0.5, by = 0.05),
                                mode = c("common", "one_at_a_time"),
                                mediators = c("PC1_R", "PC2_R", "PC3"),
                                mediators_varied = mediators) {
  mode <- match.arg(mode)
  rho_cols <- paste0("rho_", seq_along(mediators))

  make_row <- function(scenario, mediator_varied, rho, rho_vec) {
    out <- tibble::tibble(
      scenario = scenario,
      mediator_varied = mediator_varied,
      rho = rho
    )
    for (i in seq_along(rho_cols)) {
      out[[rho_cols[[i]]]] <- rho_vec[[i]]
    }
    out
  }

  if (mode == "common") {
    return(purrr::map_dfr(rho_grid, function(rho) {
      make_row("common", "all", rho, rep(rho, length(mediators)))
    }))
  }

  varied_idx <- match(intersect(mediators_varied, mediators), mediators)
  if (length(varied_idx) == 0) {
    stop("mediators_varied must include at least one supplied mediator.", call. = FALSE)
  }

  purrr::map_dfr(varied_idx, function(idx) {
    purrr::map_dfr(rho_grid, function(rho) {
      rho_vec <- rep(0, length(mediators))
      rho_vec[[idx]] <- rho
      make_row("one_at_a_time", mediators[[idx]], rho, rho_vec)
    })
  })
}

summarize_bayesian_A_draws <- function(draws, value_col = "value", group_cols) {
  draws |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      mean = mean(.data[[value_col]], na.rm = TRUE),
      median = stats::median(.data[[value_col]], na.rm = TRUE),
      sd = stats::sd(.data[[value_col]], na.rm = TRUE),
      q025 = stats::quantile(.data[[value_col]], 0.025, na.rm = TRUE),
      q975 = stats::quantile(.data[[value_col]], 0.975, na.rm = TRUE),
      Pr_gt_0 = mean(.data[[value_col]] > 0, na.rm = TRUE),
      Pr_lt_0 = mean(.data[[value_col]] < 0, na.rm = TRUE),
      .groups = "drop"
    )
}

bayes_A_wide_for_exposure <- function(residual_draws, exposure, beta = residual_draws$beta) {
  mediators <- bayes_A_bundle_mediators(residual_draws)
  n_draws <- residual_draws$n_draws
  alpha_mat <- residual_draws$alpha[[exposure]][, mediators, drop = FALSE]
  beta_mat <- beta[, mediators, drop = FALSE]
  TE <- residual_draws$total_effect[, exposure]
  outcome_direct_coef <- residual_draws$outcome_direct_coef[, exposure]
  IIE <- alpha_mat * beta_mat
  total_IIE <- rowSums(IIE)
  direct_vec <- TE - total_IIE
  PM <- total_IIE / TE

  out <- tibble::tibble(
    draw = seq_len(n_draws),
    exposure = exposure,
    contrast = bayes_A_contrast_label(exposure, residual_draws),
    direct_effect = direct_vec,
    outcome_direct_coef = outcome_direct_coef,
    total_IIE = total_IIE,
    TE = TE,
    PM = PM
  )
  for (m in mediators) {
    out[[paste0("IIE_", m)]] <- IIE[[m]]
    out[[paste0("alpha_", m)]] <- alpha_mat[[m]]
    out[[paste0("beta_", m)]] <- beta_mat[[m]]
  }
  out
}

bayes_A_path_from_wide <- function(wide) {
  wide |>
    dplyr::select(
      dplyr::any_of(c("scenario", "mediator_varied", "rho")),
      dplyr::matches("^rho_[0-9]+$"),
      draw,
      dplyr::any_of(c("valid_parameter", "feasibility")),
      exposure,
      contrast,
      dplyr::starts_with("alpha_"),
      dplyr::starts_with("beta_"),
      dplyr::starts_with("IIE_")
    ) |>
    tidyr::pivot_longer(
      cols = dplyr::matches("^(alpha|beta|IIE)_"),
      names_to = c(".value", "mediator"),
      names_pattern = "^(alpha|beta|IIE)_(.*)$"
    )
}

compute_bayesian_A_baseline_from_residual_draws <- function(residual_draws) {
  validate_bayesian_A_residual_draws(residual_draws)
  exposures <- bayes_A_bundle_exposures(residual_draws)

  wide <- dplyr::bind_rows(lapply(exposures, function(exposure) {
    bayes_A_wide_for_exposure(residual_draws, exposure)
  }))

  path_draws <- bayes_A_path_from_wide(wide)
  decomposition_draws <- wide |>
    dplyr::select(draw, exposure, contrast, direct_effect, outcome_direct_coef, total_IIE, TE, PM)

  list(
    wide = wide,
    path_draws = path_draws,
    decomposition_draws = decomposition_draws,
    path_summary = path_draws |>
      dplyr::group_by(contrast, exposure, mediator) |>
      dplyr::summarise(
        alpha_mean = mean(alpha, na.rm = TRUE),
        alpha_q025 = stats::quantile(alpha, 0.025, na.rm = TRUE),
        alpha_q975 = stats::quantile(alpha, 0.975, na.rm = TRUE),
        beta_mean = mean(beta, na.rm = TRUE),
        beta_q025 = stats::quantile(beta, 0.025, na.rm = TRUE),
        beta_q975 = stats::quantile(beta, 0.975, na.rm = TRUE),
        IIE_mean = mean(IIE, na.rm = TRUE),
        IIE_median = stats::median(IIE, na.rm = TRUE),
        IIE_q025 = stats::quantile(IIE, 0.025, na.rm = TRUE),
        IIE_q975 = stats::quantile(IIE, 0.975, na.rm = TRUE),
        Pr_IIE_gt_0 = mean(IIE > 0, na.rm = TRUE),
        Pr_IIE_lt_0 = mean(IIE < 0, na.rm = TRUE),
        .groups = "drop"
      ),
    decomposition_summary = decomposition_draws |>
      tidyr::pivot_longer(cols = c(direct_effect, outcome_direct_coef, total_IIE, TE, PM), names_to = "quantity", values_to = "value") |>
      summarize_bayesian_A_draws(group_cols = c("contrast", "exposure", "quantity")) |>
      tidyr::pivot_wider(
        names_from = quantity,
        values_from = c(mean, median, sd, q025, q975, Pr_gt_0, Pr_lt_0),
        names_glue = "{quantity}_{.value}"
      ),
    residual_draws = residual_draws,
    sampled_indices = residual_draws$sampled_indices,
    n_draws = residual_draws$n_draws
  )
}

bayes_A_sigma_M <- function(residual_draws, draw, mediators, R) {
  sig <- as.numeric(unlist(residual_draws$sig[draw, mediators, drop = TRUE], use.names = FALSE))
  Sigma_M <- diag(sig^2, length(mediators))
  if (length(mediators) > 1) {
    for (i in seq_along(mediators)) {
      for (j in seq_along(mediators)) {
        if (i != j) {
          Sigma_M[i, j] <- R[draw, i, j] * sig[[i]] * sig[[j]]
        }
      }
    }
  }
  dimnames(Sigma_M) <- list(mediators, mediators)
  Sigma_M
}

compute_bayesian_sensitivity_A_from_draws <- function(residual_draws,
                                                      rho_grid = seq(-0.5, 0.5, by = 0.05),
                                                      mode = c("common", "one_at_a_time"),
                                                      mediators = bayes_A_bundle_mediators(residual_draws),
                                                      mediators_varied = mediators,
                                                      feasibility_tol = 1e-10) {
  mode <- match.arg(mode)
  validate_bayesian_A_residual_draws(residual_draws)
  mediators <- intersect(mediators, bayes_A_bundle_mediators(residual_draws))
  if (length(mediators) == 0) {
    stop("No valid mediators were supplied for Bayesian Sensitivity A.", call. = FALSE)
  }

  exposures <- bayes_A_bundle_exposures(residual_draws)
  R <- bayes_A_residual_correlation_array(residual_draws, mediators = bayes_A_bundle_mediators(residual_draws))
  all_mediators <- bayes_A_bundle_mediators(residual_draws)
  mediator_idx <- match(mediators, all_mediators)
  R <- R[, mediator_idx, mediator_idx, drop = FALSE]
  grid <- bayesian_A_rho_grid(rho_grid, mode = mode, mediators = mediators, mediators_varied = mediators_varied)
  rho_cols <- paste0("rho_", seq_along(mediators))

  rho_matrix <- as.matrix(grid[, rho_cols, drop = FALSE])

  rows <- dplyr::bind_rows(lapply(seq_len(residual_draws$n_draws), function(s) {
    sig_vec <- as.numeric(unlist(residual_draws$sig[s, mediators, drop = TRUE], use.names = FALSE))
    beta_vec <- as.numeric(unlist(residual_draws$beta[s, mediators, drop = TRUE], use.names = FALSE))
    sigma_y <- residual_draws$sigma_y[[s]]

    Sigma_M <- matrix(R[s, , ], length(mediators)) * outer(sig_vec, sig_vec)
    adjustment <- structural_rho_adjustment(Sigma_M, sigma_y, rho_matrix, feasibility_tol)
    solve_result <- adjustment$adjustment
    feasibility <- adjustment$q
    valid <- adjustment$valid

    beta_adj <- sweep(solve_result, 2, beta_vec, function(adjustment, beta) beta - adjustment)
    beta_adj[!valid, ] <- NA_real_
    colnames(beta_adj) <- paste0("beta_", mediators)

    dplyr::bind_cols(
      grid,
      tibble::tibble(
        draw = s,
        valid_parameter = valid,
        feasibility = feasibility,
        sigma_y_structural = adjustment$sigma_y
      ),
      tibble::as_tibble(beta_adj)
    )
  }))

  make_contrast <- function(exposure) {
    alpha <- residual_draws$alpha[[exposure]][, mediators, drop = FALSE]
    TE <- residual_draws$total_effect[, exposure]
    outcome_direct_coef <- residual_draws$outcome_direct_coef[, exposure]

    out <- rows |>
      dplyr::mutate(
        exposure = exposure,
        contrast = bayes_A_contrast_label(exposure, residual_draws),
        TE = TE[draw],
        outcome_direct_coef = outcome_direct_coef[draw]
      )
    total_IIE <- rep(0, nrow(out))
    for (m in mediators) {
      out[[paste0("alpha_", m)]] <- alpha[out$draw, m]
      out[[paste0("IIE_", m)]] <- out[[paste0("alpha_", m)]] * out[[paste0("beta_", m)]]
      total_IIE <- total_IIE + out[[paste0("IIE_", m)]]
    }
    out$total_IIE <- total_IIE
    out$direct_effect <- out$TE - out$total_IIE
    out$PM <- out$total_IIE / out$TE
    out
  }

  wide <- dplyr::bind_rows(lapply(exposures, make_contrast))
  path_draws <- bayes_A_path_from_wide(wide)
  summary_draws <- wide |>
    dplyr::select(
      scenario, mediator_varied, rho, dplyr::all_of(rho_cols),
      draw, valid_parameter, feasibility, exposure, contrast,
      direct_effect, outcome_direct_coef, total_IIE, TE, PM
    )

  group_cols <- c("scenario", "mediator_varied", "rho", rho_cols, "contrast", "exposure")
  list(
    parameterization = structural_rho_version(),
    grid = grid,
    residual_draws = residual_draws,
    draws_wide = wide,
    path_draws = path_draws,
    summary_draws = summary_draws,
    path = path_draws |>
      dplyr::filter(valid_parameter) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c(group_cols, "mediator")))) |>
      dplyr::summarise(
        alpha_mean = mean(alpha, na.rm = TRUE),
        alpha_q025 = stats::quantile(alpha, 0.025, na.rm = TRUE),
        alpha_q975 = stats::quantile(alpha, 0.975, na.rm = TRUE),
        beta_mean = mean(beta, na.rm = TRUE),
        beta_q025 = stats::quantile(beta, 0.025, na.rm = TRUE),
        beta_q975 = stats::quantile(beta, 0.975, na.rm = TRUE),
        IIE_mean = mean(IIE, na.rm = TRUE),
        IIE_median = stats::median(IIE, na.rm = TRUE),
        IIE_q025 = stats::quantile(IIE, 0.025, na.rm = TRUE),
        IIE_q975 = stats::quantile(IIE, 0.975, na.rm = TRUE),
        IIE_Pr_gt_0 = mean(IIE > 0, na.rm = TRUE),
        IIE_Pr_lt_0 = mean(IIE < 0, na.rm = TRUE),
        .groups = "drop"
      ),
    summary = summary_draws |>
      dplyr::filter(valid_parameter) |>
      tidyr::pivot_longer(cols = c(direct_effect, outcome_direct_coef, total_IIE, TE, PM), names_to = "quantity", values_to = "value") |>
      summarize_bayesian_A_draws(group_cols = c(group_cols, "quantity")) |>
      tidyr::pivot_wider(
        names_from = quantity,
        values_from = c(mean, median, sd, q025, q975, Pr_gt_0, Pr_lt_0),
        names_glue = "{quantity}_{.value}"
      ),
    feasibility = wide |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c("scenario", "mediator_varied", "rho", rho_cols)))) |>
      dplyr::summarise(
        valid_draw_fraction = mean(valid_parameter),
        max_feasibility = max(feasibility, na.rm = TRUE),
        .groups = "drop"
      ),
    validation = list(
      TE_fixed_decomposition_max_abs = max(abs(wide$TE - wide$direct_effect - wide$total_IIE), na.rm = TRUE),
      rho0_available = any(rowSums(abs(grid[, rho_cols, drop = FALSE])) == 0)
    )
  )
}

bayesian_A_identity_columns <- function(mediators = c("PC1_R", "PC2_R", "PC3")) {
  c(
    "direct_effect", "outcome_direct_coef", "total_IIE", "TE", "PM",
    paste0("alpha_", mediators),
    paste0("beta_", mediators),
    paste0("IIE_", mediators)
  )
}

validate_bayesian_A_rho0_draw_identity <- function(A_result,
                                                   baseline_quantities,
                                                   tolerance = 1e-12) {
  mediators <- bayes_A_bundle_mediators(A_result$residual_draws %||% baseline_quantities$residual_draws %||% list())
  if (is.null(mediators)) {
    mediators <- sub("^alpha_", "", grep("^alpha_", names(baseline_quantities$wide), value = TRUE))
  }
  cols <- bayesian_A_identity_columns(mediators)
  rho_cols <- grep("^rho_[0-9]+$", names(A_result$draws_wide), value = TRUE)

  A0 <- A_result$draws_wide |>
    dplyr::filter(rowSums(abs(dplyr::pick(dplyr::all_of(rho_cols)))) == 0) |>
    dplyr::select(
      scenario, mediator_varied, rho, dplyr::all_of(rho_cols),
      draw, exposure, contrast, dplyr::all_of(cols)
    )
  B0 <- baseline_quantities$wide |>
    dplyr::select(draw, exposure, contrast, dplyr::all_of(cols)) |>
    dplyr::rename_with(~ paste0("baseline_", .x), dplyr::all_of(cols))

  joined <- A0 |>
    dplyr::left_join(B0, by = c("draw", "exposure", "contrast"))

  detail <- purrr::map_dfr(cols, function(col) {
    diff <- joined[[col]] - joined[[paste0("baseline_", col)]]
    tibble::tibble(
      scenario = joined$scenario,
      mediator_varied = joined$mediator_varied,
      draw = joined$draw,
      exposure = joined$exposure,
      contrast = joined$contrast,
      quantity = col,
      A_value = joined[[col]],
      baseline_value = joined[[paste0("baseline_", col)]],
      diff = diff,
      abs_diff = abs(diff)
    )
  })

  max_abs_diff <- max(detail$abs_diff, na.rm = TRUE)
  list(
    passed = is.finite(max_abs_diff) && max_abs_diff <= tolerance,
    tolerance = tolerance,
    max_abs_draw_diff = max_abs_diff,
    max_discrepancy = detail |>
      dplyr::arrange(dplyr::desc(abs_diff)) |>
      dplyr::slice(1),
    detail = detail
  )
}

validate_bayesian_A_rho0_scenarios_identical <- function(A_common,
                                                         A_one_at_a_time,
                                                         tolerance = 1e-12) {
  mediators <- bayes_A_bundle_mediators(A_common$residual_draws)
  cols <- bayesian_A_identity_columns(mediators)
  rho_cols <- grep("^rho_[0-9]+$", names(A_one_at_a_time$draws_wide), value = TRUE)

  common0 <- A_common$draws_wide |>
    dplyr::filter(rowSums(abs(dplyr::pick(dplyr::all_of(grep("^rho_[0-9]+$", names(A_common$draws_wide), value = TRUE))))) == 0) |>
    dplyr::select(draw, exposure, contrast, dplyr::all_of(cols)) |>
    dplyr::rename_with(~ paste0("common_", .x), dplyr::all_of(cols))
  one0 <- A_one_at_a_time$draws_wide |>
    dplyr::filter(rowSums(abs(dplyr::pick(dplyr::all_of(rho_cols)))) == 0) |>
    dplyr::select(
      scenario, mediator_varied, rho, dplyr::all_of(rho_cols),
      draw, exposure, contrast, dplyr::all_of(cols)
    )
  joined <- one0 |>
    dplyr::left_join(common0, by = c("draw", "exposure", "contrast"))

  detail <- purrr::map_dfr(cols, function(col) {
    diff <- joined[[col]] - joined[[paste0("common_", col)]]
    tibble::tibble(
      mediator_varied = joined$mediator_varied,
      draw = joined$draw,
      exposure = joined$exposure,
      contrast = joined$contrast,
      quantity = col,
      one_at_a_time_value = joined[[col]],
      common_value = joined[[paste0("common_", col)]],
      diff = diff,
      abs_diff = abs(diff)
    )
  })
  max_abs_diff <- max(detail$abs_diff, na.rm = TRUE)

  list(
    passed = is.finite(max_abs_diff) && max_abs_diff <= tolerance,
    tolerance = tolerance,
    max_abs_draw_diff = max_abs_diff,
    max_discrepancy = detail |>
      dplyr::arrange(dplyr::desc(abs_diff)) |>
      dplyr::slice(1),
    detail = detail
  )
}
