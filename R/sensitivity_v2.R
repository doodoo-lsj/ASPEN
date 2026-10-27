# Sensitivity helpers for the current v2 spatial mediation model.
#
# This file intentionally does not modify or source sensitivity_0611.R.
# It ports the reusable sensitivity logic to the current analysis defaults:
# - data: pt16_emory_GEX_immune_FULL_v2.rds
# - exposure: X_adjacent and X_far, inside as reference
# - mediators: PC1_R, PC2_R, PC3
# - covariates: x_coord and y_coord

source(local({
  source_paths <- Filter(function(x) is.character(x) && length(x) == 1L,
                         lapply(sys.frames(), function(frame) frame$ofile))
  file.path(dirname(tail(source_paths, 1L)[[1L]]), "structural_residual_correlation.R")
}), local = TRUE)

prepare_v2_mediation_data <- function(rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
                                      coord_x = "imagecol",
                                      coord_y = "imagerow") {
  gex <- readRDS(rds_path)

  required <- c("X_cat", "Y", "meta", "M_expr")
  missing_required <- setdiff(required, names(gex))
  if (length(missing_required) > 0) {
    stop("RDS is missing required elements: ", paste(missing_required, collapse = ", "))
  }

  meta <- as.data.frame(gex$meta)
  missing_coords <- setdiff(c(coord_x, coord_y), names(meta))
  if (length(missing_coords) > 0) {
    stop("Metadata is missing coordinate columns: ", paste(missing_coords, collapse = ", "))
  }

  eda_df <- meta |>
    dplyr::mutate(
      X_cat = factor(gex$X_cat, levels = c("inside", "adjacent", "far")),
      Y = as.numeric(gex$Y),
      x_coord = .data[[coord_x]],
      y_coord = .data[[coord_y]]
    )

  m_scaled <- scale(gex$M_expr)
  pca_fit <- prcomp(m_scaled, center = FALSE, scale. = FALSE)

  pca_df <- eda_df |>
    dplyr::mutate(
      PC1_original = pca_fit$x[, 1],
      PC2_original = pca_fit$x[, 2],
      PC3_original = pca_fit$x[, 3],
      PC1_R = -PC1_original,
      PC2_R = -PC2_original,
      PC3 = PC3_original
    )

  dat_full_allpc <- pca_df |>
    dplyr::mutate(
      X_cat = factor(X_cat, levels = c("inside", "adjacent", "far")),
      X_adjacent = as.integer(X_cat == "adjacent"),
      X_far = as.integer(X_cat == "far")
    ) |>
    dplyr::select(
      Y, X_cat, X_adjacent, X_far,
      PC1_R, PC2_R, PC3,
      x_coord, y_coord
    ) |>
    tidyr::drop_na()

  list(
    gex = gex,
    pca_fit = pca_fit,
    pca_df = pca_df,
    dat_full_allpc = dat_full_allpc
  )
}

decompose_linear_multix <- function(data,
                                    exposures = c("X_adjacent", "X_far"),
                                    outcome = "Y",
                                    mediators = c("PC1_R", "PC2_R", "PC3"),
                                    covariates = c("x_coord", "y_coord")) {
  cov_part <- if (!is.null(covariates) && length(covariates) > 0) {
    paste(covariates, collapse = " + ")
  } else {
    NULL
  }

  exposure_part <- paste(exposures, collapse = " + ")
  mediator_part <- paste(mediators, collapse = " + ")

  te_formula <- stats::as.formula(
    paste(
      outcome, "~", exposure_part,
      if (!is.null(cov_part)) paste("+", cov_part) else ""
    )
  )
  fit_te <- stats::lm(te_formula, data = data)
  te_coef <- stats::coef(fit_te)[exposures]

  alpha_tbl <- purrr::map_dfr(mediators, function(med) {
    m_formula <- stats::as.formula(
      paste(
        med, "~", exposure_part,
        if (!is.null(cov_part)) paste("+", cov_part) else ""
      )
    )

    fit_m <- stats::lm(m_formula, data = data)

    tibble::tibble(
      mediator = med,
      exposure = exposures,
      alpha_X_to_M = stats::coef(fit_m)[exposures]
    )
  })

  y_formula <- stats::as.formula(
    paste(
      outcome, "~", exposure_part, "+", mediator_part,
      if (!is.null(cov_part)) paste("+", cov_part) else ""
    )
  )
  fit_y <- stats::lm(y_formula, data = data)

  beta_tbl <- tibble::tibble(
    mediator = mediators,
    beta_M_to_Y = stats::coef(fit_y)[mediators]
  )

  path_tbl <- alpha_tbl |>
    dplyr::left_join(beta_tbl, by = "mediator") |>
    dplyr::mutate(indirect_component = alpha_X_to_M * beta_M_to_Y)

  summary_tbl <- path_tbl |>
    dplyr::group_by(exposure) |>
    dplyr::summarise(NIE = sum(indirect_component), .groups = "drop") |>
    dplyr::mutate(
      TE = te_coef[exposure],
      NDE = TE - NIE,
      outcome_direct_coef = stats::coef(fit_y)[exposure],
      PM = NIE / TE,
      mediators = paste(mediators, collapse = " + "),
      covariates = ifelse(
        is.null(covariates) || length(covariates) == 0,
        "none",
        paste(covariates, collapse = " + ")
      )
    ) |>
    dplyr::select(exposure, mediators, covariates, TE, NDE, NIE, PM, outcome_direct_coef)

  list(
    summary = summary_tbl,
    path = path_tbl,
    fit_total = fit_te,
    fit_outcome = fit_y
  )
}

residual_corr_sensitivity <- function(data,
                                      exposures = c("X_adjacent", "X_far"),
                                      outcome = "Y",
                                      mediators = c("PC1_R", "PC2_R", "PC3"),
                                      covariates = c("x_coord", "y_coord"),
                                      rho_grid = seq(-0.5, 0.5, by = 0.01),
                                      mode = c("common", "one_at_a_time"),
                                      feasibility_tol = 1e-10) {
  mode <- match.arg(mode)

  cov_part <- if (!is.null(covariates) && length(covariates) > 0) {
    paste(covariates, collapse = " + ")
  } else {
    NULL
  }

  exposure_part <- paste(exposures, collapse = " + ")
  mediator_part <- paste(mediators, collapse = " + ")

  te_formula <- stats::as.formula(
    paste(
      outcome, "~", exposure_part,
      if (!is.null(cov_part)) paste("+", cov_part) else ""
    )
  )
  fit_te <- stats::lm(te_formula, data = data)
  te <- stats::coef(fit_te)[exposures]

  mediator_fits <- lapply(mediators, function(med) {
    m_formula <- stats::as.formula(
      paste(
        med, "~", exposure_part,
        if (!is.null(cov_part)) paste("+", cov_part) else ""
      )
    )
    stats::lm(m_formula, data = data)
  })
  names(mediator_fits) <- mediators

  alpha_mat <- sapply(mediators, function(med) {
    stats::coef(mediator_fits[[med]])[exposures]
  })
  alpha_mat <- matrix(
    alpha_mat,
    nrow = length(exposures),
    ncol = length(mediators),
    dimnames = list(exposures, mediators)
  )
  rownames(alpha_mat) <- exposures
  colnames(alpha_mat) <- mediators

  mediator_resid <- sapply(mediators, function(med) stats::resid(mediator_fits[[med]]))
  mediator_resid <- matrix(
    mediator_resid,
    ncol = length(mediators),
    dimnames = list(NULL, mediators)
  )
  colnames(mediator_resid) <- mediators

  sigma_m <- stats::cov(mediator_resid)
  sigma_m <- matrix(
    sigma_m,
    nrow = length(mediators),
    ncol = length(mediators),
    dimnames = list(mediators, mediators)
  )
  sd_m <- stats::setNames(apply(mediator_resid, 2, stats::sd), mediators)

  y_formula <- stats::as.formula(
    paste(
      outcome, "~", exposure_part, "+", mediator_part,
      if (!is.null(cov_part)) paste("+", cov_part) else ""
    )
  )
  fit_y <- stats::lm(y_formula, data = data)
  beta_hat <- stats::coef(fit_y)[mediators]
  nde_hat <- stats::coef(fit_y)[exposures]
  sd_y_resid <- stats::sd(stats::resid(fit_y))

  nie_hat <- as.vector(alpha_mat %*% beta_hat)
  names(nie_hat) <- exposures

  observed_summary <- tibble::tibble(
    exposure = exposures,
    TE = te[exposures],
    NDE = te[exposures] - nie_hat[exposures],
    NIE = nie_hat[exposures],
    PM = nie_hat[exposures] / te[exposures],
    outcome_direct_coef = nde_hat[exposures]
  )

  observed_path <- tibble::as_tibble(alpha_mat, rownames = "exposure") |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(mediators),
      names_to = "mediator",
      values_to = "alpha_X_to_M"
    ) |>
    dplyr::left_join(
      tibble::tibble(mediator = mediators, beta_M_to_Y = beta_hat[mediators]),
      by = "mediator"
    ) |>
    dplyr::mutate(indirect_component = alpha_X_to_M * beta_M_to_Y)

  make_result_for_rho <- function(rho_vec, label, scenario, varied_mediator = "all_mediators") {
    adjustment <- structural_rho_adjustment(sigma_m, sd_y_resid, rho_vec, feasibility_tol)
    cov_my <- stats::setNames(as.numeric(adjustment$covariance), mediators)
    feasibility_ratio <- adjustment$q[[1]]
    valid_parameter <- adjustment$valid[[1]]

    beta_adj <- if (valid_parameter) {
      as.vector(beta_hat - adjustment$adjustment[1, ])
    } else {
      rep(NA_real_, length(mediators))
    }
    names(beta_adj) <- mediators

    nie_adj <- if (valid_parameter) {
      as.vector(alpha_mat %*% beta_adj)
    } else {
      rep(NA_real_, length(exposures))
    }
    names(nie_adj) <- exposures

    rho_tbl <- tibble::as_tibble_row(stats::setNames(as.list(rho_vec), paste0("rho_", mediators)))

    summary_tbl <- tibble::tibble(
      scenario = scenario,
      sensitivity_type = label,
      varied_mediator = varied_mediator,
      exposure = exposures,
      TE = te[exposures],
      NIE_adjusted = nie_adj[exposures],
      NDE_adjusted = te[exposures] - nie_adj[exposures],
      PM_adjusted = nie_adj[exposures] / te[exposures],
      feasibility_ratio = feasibility_ratio,
      valid_parameter = valid_parameter
    ) |>
      dplyr::bind_cols(rho_tbl[rep(1, length(exposures)), , drop = FALSE])

    beta_tbl <- tibble::tibble(
      scenario = scenario,
      sensitivity_type = label,
      varied_mediator = varied_mediator,
      mediator = mediators,
      beta_observed = beta_hat[mediators],
      beta_adjusted = beta_adj[mediators],
      residual_covariance = cov_my[mediators],
      feasibility_ratio = feasibility_ratio,
      valid_parameter = valid_parameter
    ) |>
      dplyr::bind_cols(rho_tbl[rep(1, length(mediators)), , drop = FALSE])

    path_base <- tibble::as_tibble(alpha_mat, rownames = "exposure") |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(mediators),
        names_to = "mediator",
        values_to = "alpha_X_to_M"
      ) |>
      dplyr::left_join(
        tibble::tibble(
          mediator = mediators,
          beta_observed = beta_hat[mediators],
          beta_adjusted = beta_adj[mediators],
          residual_covariance = cov_my[mediators]
        ),
        by = "mediator"
      ) |>
      dplyr::mutate(
        scenario = scenario,
        sensitivity_type = label,
        varied_mediator = varied_mediator,
        indirect_observed = alpha_X_to_M * beta_observed,
        indirect_adjusted = alpha_X_to_M * beta_adjusted,
        feasibility_ratio = feasibility_ratio,
        valid_parameter = valid_parameter
      )

    path_tbl <- dplyr::bind_cols(
      path_base,
      rho_tbl[rep(1, nrow(path_base)), , drop = FALSE]
    ) |>
      dplyr::select(
        scenario,
        sensitivity_type,
        varied_mediator,
        exposure,
        mediator,
        dplyr::starts_with("rho_"),
        alpha_X_to_M,
        beta_observed,
        beta_adjusted,
        residual_covariance,
        indirect_observed,
        indirect_adjusted,
        feasibility_ratio,
        valid_parameter
      )

    list(summary = summary_tbl, beta = beta_tbl, path = path_tbl)
  }

  results <- list()
  beta_results <- list()
  path_results <- list()

  if (mode == "common") {
    for (rho in rho_grid) {
      rho_vec <- stats::setNames(rep(rho, length(mediators)), mediators)
      tmp <- make_result_for_rho(
        rho_vec,
        "common rho for all mediators",
        scenario = "common",
        varied_mediator = "all_mediators"
      )
      results[[as.character(rho)]] <- dplyr::mutate(tmp$summary, rho = rho)
      beta_results[[as.character(rho)]] <- dplyr::mutate(tmp$beta, rho = rho)
      path_results[[as.character(rho)]] <- dplyr::mutate(tmp$path, rho = rho)
    }
  } else {
    counter <- 1
    for (med in mediators) {
      for (rho in rho_grid) {
        rho_vec <- stats::setNames(rep(0, length(mediators)), mediators)
        rho_vec[med] <- rho
        tmp <- make_result_for_rho(
          rho_vec,
          paste0("rho varied for ", med, " only"),
          scenario = "one_at_a_time",
          varied_mediator = med
        )
        results[[counter]] <- dplyr::mutate(tmp$summary, varied_mediator = med, rho = rho)
        beta_results[[counter]] <- dplyr::mutate(tmp$beta, varied_mediator = med, rho = rho)
        path_results[[counter]] <- dplyr::mutate(tmp$path, varied_mediator = med, rho = rho)
        counter <- counter + 1
      }
    }
  }

  list(
    sensitivity_parameterization = structural_rho_version(),
    observed_summary = observed_summary,
    observed_path = observed_path,
    sensitivity_summary = dplyr::bind_rows(results),
    sensitivity_beta = dplyr::bind_rows(beta_results),
    sensitivity_path = dplyr::bind_rows(path_results),
    alpha_mat = alpha_mat,
    beta_hat = beta_hat,
    Sigma_M = sigma_m,
    sd_M = sd_m,
    sd_Y_resid = sd_y_resid,
    fit_total = fit_te,
    fit_outcome = fit_y,
    mediator_fits = mediator_fits
  )
}

validate_rho0_matches_observed <- function(sensitivity_result, tolerance = 1e-8) {
  summary_rho0 <- sensitivity_result$sensitivity_summary |>
    dplyr::filter(abs(rho) < tolerance, valid_parameter) |>
    dplyr::select(exposure, NIE_adjusted) |>
    dplyr::left_join(
      dplyr::select(sensitivity_result$observed_summary, exposure, NIE),
      by = "exposure"
    ) |>
    dplyr::mutate(abs_diff = abs(NIE_adjusted - NIE))

  path_rho0 <- sensitivity_result$sensitivity_path |>
    dplyr::filter(abs(rho) < tolerance, valid_parameter) |>
    dplyr::select(exposure, mediator, beta_adjusted, indirect_adjusted) |>
    dplyr::left_join(
      dplyr::select(
        sensitivity_result$observed_path,
        exposure,
        mediator,
        beta_M_to_Y,
        indirect_component
      ),
      by = c("exposure", "mediator")
    ) |>
    dplyr::mutate(
      beta_abs_diff = abs(beta_adjusted - beta_M_to_Y),
      indirect_abs_diff = abs(indirect_adjusted - indirect_component)
    )

  list(
    summary = summary_rho0,
    path = path_rho0,
    passed = all(summary_rho0$abs_diff <= tolerance, na.rm = TRUE) &&
      all(path_rho0$beta_abs_diff <= tolerance, na.rm = TRUE) &&
      all(path_rho0$indirect_abs_diff <= tolerance, na.rm = TRUE)
  )
}

find_rho_zero_crossings <- function(data,
                                    value_col,
                                    group_cols,
                                    rho_col = "rho") {
  split_key <- interaction(data[, group_cols, drop = FALSE], drop = TRUE)

  purrr::imap_dfr(split(data, split_key), function(df, key) {
    df <- df |>
      dplyr::filter(valid_parameter, !is.na(.data[[value_col]]), !is.na(.data[[rho_col]])) |>
      dplyr::arrange(.data[[rho_col]])

    group_values <- if (nrow(df) > 0) {
      df[1, group_cols, drop = FALSE]
    } else {
      tibble::as_tibble_row(stats::setNames(rep(NA_character_, length(group_cols)), group_cols))
    }

    if (nrow(df) == 0) {
      return(dplyr::bind_cols(
        group_values,
        tibble::tibble(
          value = value_col,
          rho_zero_crossing = NA_real_,
          crossing_in_grid = FALSE
        )
      ))
    }

    rho <- df[[rho_col]]
    value <- df[[value_col]]
    exact_zero <- which(abs(value) < .Machine$double.eps^0.5)

    if (length(exact_zero) > 0) {
      crossing <- rho[exact_zero[1]]
      crossing_found <- TRUE
    } else {
      cross_idx <- which(value[-length(value)] * value[-1] < 0)
      crossing_found <- length(cross_idx) > 0
      crossing <- if (crossing_found) {
        i <- cross_idx[1]
        rho[i] - value[i] * (rho[i + 1] - rho[i]) / (value[i + 1] - value[i])
      } else {
        NA_real_
      }
    }

    dplyr::bind_cols(
      group_values,
      tibble::tibble(
        value = value_col,
        rho_zero_crossing = crossing,
        crossing_in_grid = crossing_found
      )
    )
  })
}

summarize_rho_tipping_points <- function(sensitivity_result) {
  total_tipping <- find_rho_zero_crossings(
    sensitivity_result$sensitivity_summary,
    value_col = "NIE_adjusted",
    group_cols = c("scenario", "sensitivity_type", "varied_mediator", "exposure")
  )

  mediator_tipping <- find_rho_zero_crossings(
    sensitivity_result$sensitivity_path,
    value_col = "indirect_adjusted",
    group_cols = c("scenario", "sensitivity_type", "varied_mediator", "exposure", "mediator")
  )

  list(total = total_tipping, mediator_specific = mediator_tipping)
}

build_rho_sensitivity_plots <- function(common_result,
                                        one_at_a_time_result = NULL) {
  common_total <- common_result$sensitivity_summary |>
    dplyr::filter(valid_parameter) |>
    dplyr::mutate(
      contrast = dplyr::case_when(
        exposure == "X_adjacent" ~ "adjacent vs inside",
        exposure == "X_far" ~ "far vs inside",
        TRUE ~ exposure
      )
    )

  p_total_common <- ggplot2::ggplot(
    common_total,
    ggplot2::aes(x = rho, y = NIE_adjusted, color = contrast)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::labs(
      x = "common rho",
      y = "Adjusted total indirect effect",
      color = "Contrast"
    ) +
    ggplot2::theme_bw(base_size = 12)

  common_path <- common_result$sensitivity_path |>
    dplyr::filter(valid_parameter) |>
    dplyr::mutate(
      contrast = dplyr::case_when(
        exposure == "X_adjacent" ~ "adjacent vs inside",
        exposure == "X_far" ~ "far vs inside",
        TRUE ~ exposure
      )
    )

  p_path_common <- ggplot2::ggplot(
    common_path,
    ggplot2::aes(x = rho, y = indirect_adjusted, color = mediator)
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~ contrast, scales = "free_y") +
    ggplot2::labs(
      x = "common rho",
      y = "Adjusted mediator-specific indirect effect",
      color = "Mediator"
    ) +
    ggplot2::theme_bw(base_size = 12)

  plots <- list(
    total_common = p_total_common,
    mediator_specific_common = p_path_common
  )

  if (!is.null(one_at_a_time_result)) {
    one_total <- one_at_a_time_result$sensitivity_summary |>
      dplyr::filter(valid_parameter) |>
      dplyr::mutate(
        contrast = dplyr::case_when(
          exposure == "X_adjacent" ~ "adjacent vs inside",
          exposure == "X_far" ~ "far vs inside",
          TRUE ~ exposure
        )
      )

    plots$total_one_at_a_time <- ggplot2::ggplot(
      one_total,
      ggplot2::aes(x = rho, y = NIE_adjusted, color = varied_mediator)
    ) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::facet_wrap(~ contrast, scales = "free_y") +
      ggplot2::labs(
        x = "rho for one mediator",
        y = "Adjusted total indirect effect",
        color = "Varied mediator"
      ) +
      ggplot2::theme_bw(base_size = 12)

    one_path <- one_at_a_time_result$sensitivity_path |>
      dplyr::filter(valid_parameter, mediator == varied_mediator) |>
      dplyr::mutate(
        contrast = dplyr::case_when(
          exposure == "X_adjacent" ~ "adjacent vs inside",
          exposure == "X_far" ~ "far vs inside",
          TRUE ~ exposure
        )
      )

    plots$mediator_specific_one_at_a_time <- ggplot2::ggplot(
      one_path,
      ggplot2::aes(x = rho, y = indirect_adjusted, color = mediator)
    ) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::facet_wrap(~ contrast, scales = "free_y") +
      ggplot2::labs(
        x = "rho for one mediator",
        y = "Adjusted mediator-specific indirect effect",
        color = "Mediator"
      ) +
      ggplot2::theme_bw(base_size = 12)
  }

  plots
}

partial_r2_sensitivity <- function(data,
                                   exposures = c("X_adjacent", "X_far"),
                                   outcome = "Y",
                                   mediators = c("PC1_R", "PC2_R", "PC3"),
                                   covariates = c("x_coord", "y_coord"),
                                   r2_m_grid = seq(0, 0.30, by = 0.01),
                                   r2_y_grid = seq(0, 0.30, by = 0.01),
                                   mode = c("common", "one_at_a_time"),
                                   target_contrast = c("X_adjacent", "X_far")) {
  mode <- match.arg(mode)
  target_contrast <- match.arg(target_contrast)

  cov_part <- if (!is.null(covariates) && length(covariates) > 0) {
    paste(covariates, collapse = " + ")
  } else {
    NULL
  }

  exposure_part <- paste(exposures, collapse = " + ")
  mediator_part <- paste(mediators, collapse = " + ")

  te_formula <- stats::as.formula(
    paste(
      outcome, "~", exposure_part,
      if (!is.null(cov_part)) paste("+", cov_part) else ""
    )
  )
  fit_te <- stats::lm(te_formula, data = data)
  te <- stats::coef(fit_te)[exposures]

  mediator_fits <- lapply(mediators, function(med) {
    m_formula <- stats::as.formula(
      paste(
        med, "~", exposure_part,
        if (!is.null(cov_part)) paste("+", cov_part) else ""
      )
    )
    stats::lm(m_formula, data = data)
  })
  names(mediator_fits) <- mediators

  alpha_mat <- sapply(mediators, function(med) {
    stats::coef(mediator_fits[[med]])[exposures]
  })
  rownames(alpha_mat) <- exposures
  colnames(alpha_mat) <- mediators

  y_formula <- stats::as.formula(
    paste(
      outcome, "~", exposure_part, "+", mediator_part,
      if (!is.null(cov_part)) paste("+", cov_part) else ""
    )
  )
  fit_y <- stats::lm(y_formula, data = data)
  fit_y_sum <- summary(fit_y)

  beta_hat <- stats::coef(fit_y)[mediators]
  se_beta <- stats::coef(fit_y_sum)[mediators, "Std. Error"]
  df_resid <- stats::df.residual(fit_y)

  nie_hat <- as.vector(alpha_mat %*% beta_hat)
  names(nie_hat) <- exposures

  observed_summary <- tibble::tibble(
    exposure = exposures,
    TE = te[exposures],
    NDE = te[exposures] - nie_hat[exposures],
    NIE = nie_hat[exposures],
    PM = nie_hat[exposures] / te[exposures],
    outcome_direct_coef = stats::coef(fit_y)[exposures]
  )

  grid <- expand.grid(r2_m = r2_m_grid, r2_y = r2_y_grid) |>
    tibble::as_tibble()

  make_result <- function(r2_m, r2_y, varied_mediator = NA_character_) {
    beta_adj <- beta_hat
    affected_mediators <- if (mode == "common") mediators else varied_mediator

    for (med in affected_mediators) {
      bias_mag <- se_beta[med] * sqrt(df_resid * r2_y * r2_m / (1 - r2_m))
      alpha_target <- alpha_mat[target_contrast, med]

      if (alpha_target * beta_hat[med] >= 0) {
        beta_adj[med] <- beta_hat[med] - sign(beta_hat[med]) * bias_mag
      } else {
        beta_adj[med] <- beta_hat[med] + sign(beta_hat[med]) * bias_mag
      }
    }

    nie_adj <- as.vector(alpha_mat %*% beta_adj)
    names(nie_adj) <- exposures

    tibble::tibble(
      varied_mediator = ifelse(is.na(varied_mediator), "common mediators", varied_mediator),
      r2_m = r2_m,
      r2_y = r2_y,
      exposure = exposures,
      TE = te[exposures],
      NIE_adjusted = nie_adj[exposures],
      NDE_adjusted = te[exposures] - nie_adj[exposures],
      PM_adjusted = nie_adj[exposures] / te[exposures]
    )
  }

  sensitivity_summary <- if (mode == "common") {
    purrr::pmap_dfr(list(grid$r2_m, grid$r2_y), ~ make_result(..1, ..2))
  } else {
    purrr::map_dfr(mediators, function(med) {
      purrr::pmap_dfr(list(grid$r2_m, grid$r2_y), ~ make_result(..1, ..2, med))
    })
  }

  list(
    observed_summary = observed_summary,
    sensitivity_summary = sensitivity_summary,
    alpha_mat = alpha_mat,
    beta_hat = beta_hat,
    se_beta = se_beta,
    df_resid = df_resid,
    fit_total = fit_te,
    fit_outcome = fit_y,
    mediator_fits = mediator_fits
  )
}

partial_r2_added_covariates <- function(data, outcome, base_vars, added_vars) {
  base_formula <- stats::as.formula(paste(outcome, "~", paste(base_vars, collapse = " + ")))
  full_formula <- stats::as.formula(
    paste(outcome, "~", paste(c(base_vars, added_vars), collapse = " + "))
  )

  fit_base <- stats::lm(base_formula, data = data)
  fit_full <- stats::lm(full_formula, data = data)

  r2_base <- summary(fit_base)$r.squared
  r2_full <- summary(fit_full)$r.squared

  tibble::tibble(
    outcome = outcome,
    base_vars = paste(base_vars, collapse = " + "),
    added_vars = paste(added_vars, collapse = " + "),
    r2_base = r2_base,
    r2_full = r2_full,
    partial_r2_added = (r2_full - r2_base) / (1 - r2_base)
  )
}

find_tipping_points <- function(sens_table) {
  sens_table |>
    dplyr::group_by(model, contrast) |>
    dplyr::arrange(rho, .by_group = TRUE) |>
    dplyr::mutate(
      NIE_at_rho0 = NIE_adjusted[rho == 0][1],
      abs_NIE_ratio = abs(NIE_adjusted) / abs(NIE_at_rho0),
      sign_change = sign(NIE_adjusted) != sign(NIE_at_rho0),
      below_half = abs_NIE_ratio <= 0.5
    ) |>
    dplyr::summarise(
      observed_NIE = NIE_at_rho0[1],
      rho_zero_crossing = ifelse(any(sign_change, na.rm = TRUE), rho[which(sign_change)[1]], NA_real_),
      rho_half_effect = ifelse(any(below_half, na.rm = TRUE), rho[which(below_half)[1]], NA_real_),
      min_NIE = min(NIE_adjusted, na.rm = TRUE),
      max_NIE = max(NIE_adjusted, na.rm = TRUE),
      .groups = "drop"
    )
}

find_pr2_tipping <- function(pr2_table) {
  pr2_table2 <- dplyr::mutate(pr2_table, r2_sum = r2_m + r2_y)

  split(pr2_table2, interaction(pr2_table2$model, pr2_table2$target, drop = TRUE)) |>
    purrr::map_dfr(function(df) {
      df <- dplyr::arrange(df, r2_m, r2_y)

      observed_nie <- df |>
        dplyr::filter(abs(r2_m) < 1e-12, abs(r2_y) < 1e-12) |>
        dplyr::slice(1) |>
        dplyr::pull(NIE_adjusted)

      if (length(observed_nie) == 0) {
        observed_nie <- NA_real_
      }

      df <- df |>
        dplyr::mutate(
          abs_ratio = abs(NIE_adjusted) / abs(observed_nie),
          below_half = abs_ratio <= 0.5,
          sign_change = sign(NIE_adjusted) != sign(observed_nie)
        )

      sym_df <- df |>
        dplyr::filter(abs(r2_m - r2_y) < 1e-12) |>
        dplyr::arrange(r2_m)

      sym_half <- dplyr::slice(dplyr::filter(sym_df, below_half), 1)
      sym_zero <- dplyr::slice(dplyr::filter(sym_df, sign_change), 1)
      min_half <- df |>
        dplyr::filter(below_half) |>
        dplyr::arrange(r2_sum, r2_m, r2_y) |>
        dplyr::slice(1)
      min_zero <- df |>
        dplyr::filter(sign_change) |>
        dplyr::arrange(r2_sum, r2_m, r2_y) |>
        dplyr::slice(1)

      tibble::tibble(
        model = df$model[1],
        target = df$target[1],
        observed_NIE = observed_nie,
        symmetric_r2_half = if (nrow(sym_half) == 0) NA_real_ else sym_half$r2_m,
        symmetric_r2_zero = if (nrow(sym_zero) == 0) NA_real_ else sym_zero$r2_m,
        min_sum_r2_half = if (nrow(min_half) == 0) NA_real_ else min_half$r2_sum,
        r2_m_half_min_sum = if (nrow(min_half) == 0) NA_real_ else min_half$r2_m,
        r2_y_half_min_sum = if (nrow(min_half) == 0) NA_real_ else min_half$r2_y,
        min_sum_r2_zero = if (nrow(min_zero) == 0) NA_real_ else min_zero$r2_sum,
        r2_m_zero_min_sum = if (nrow(min_zero) == 0) NA_real_ else min_zero$r2_m,
        r2_y_zero_min_sum = if (nrow(min_zero) == 0) NA_real_ else min_zero$r2_y
      )
    }) |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4))) |>
    dplyr::arrange(model, target)
}
