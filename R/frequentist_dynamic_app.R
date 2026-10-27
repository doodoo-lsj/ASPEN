# Dynamic frequentist mediation runtime helpers for prepared Shiny data.
#
# This layer mirrors the Shiny Bayesian runtime data specification, but uses the
# existing lm-based frequentist decomposition and residual-rho sensitivity.

source(local({
  source_paths <- Filter(function(x) is.character(x) && length(x) == 1L,
                         lapply(sys.frames(), function(frame) frame$ofile))
  file.path(dirname(tail(source_paths, 1L)[[1L]]), "structural_residual_correlation.R")
}), local = TRUE)

frequentist_dynamic_model_version <- function() {
  "frequentist_mediation_dynamic_v3_original_pca_structural_rho"
}

frequentist_dynamic_prepare_data <- function(analysis_obj, selected_pcs, pc_signs) {
  prepared <- bayesian_dynamic_prepare_data(
    analysis_obj = analysis_obj,
    selected_pcs = selected_pcs,
    pc_signs = pc_signs,
    standardize_numeric_covariates = TRUE
  )
  prepared$model_version <- frequentist_dynamic_model_version()
  prepared
}

map_frequentist_mediator_names <- function(df, mediator_map, columns = c("mediator", "varied_mediator")) {
  name_map <- stats::setNames(mediator_map$mediator, mediator_map$model_name)
  out <- df
  for (column in intersect(columns, names(out))) {
    mapped <- unname(name_map[out[[column]]])
    out[[column]] <- ifelse(is.na(mapped), out[[column]], mapped)
  }
  out
}

make_frequentist_dynamic_fit <- function(analysis_obj, selected_pcs, pc_signs) {
  prepared <- frequentist_dynamic_prepare_data(analysis_obj, selected_pcs, pc_signs)
  fit <- decompose_linear_multix(
    data = prepared$data,
    exposures = prepared$exposures,
    outcome = "Y",
    mediators = prepared$mediator_model_names,
    covariates = prepared$covariate_terms
  )

  fit$summary <- fit$summary |>
    dplyr::mutate(contrast = unname(prepared$exposure_labels[exposure]))
  fit$summary$contrast <- ifelse(is.na(fit$summary$contrast), fit$summary$exposure, fit$summary$contrast)
  fit$path <- fit$path |>
    map_frequentist_mediator_names(prepared$mediator_map) |>
    dplyr::mutate(contrast = unname(prepared$exposure_labels[exposure]))
  fit$path$contrast <- ifelse(is.na(fit$path$contrast), fit$path$exposure, fit$path$contrast)

  list(
    source = "dynamic_fit",
    model_id = frequentist_dynamic_model_version(),
    prepared = prepared,
    fit = fit,
    selected_pcs = selected_pcs,
    created_at = Sys.time()
  )
}

make_frequentist_example_fit <- function(analysis_obj, pc_signs) {
  result <- make_frequentist_dynamic_fit(
    analysis_obj, selected_pcs = c("PC1", "PC2", "PC3"), pc_signs = pc_signs
  )
  result$source <- "validated_example"
  result$model_id <- "frequentist_example_current"
  result
}

frequentist_dynamic_cache_key <- function(analysis_obj, selected_pcs, pc_signs) {
  digest::digest(
    list(
      model_version = frequentist_dynamic_model_version(),
      data_source = analysis_obj$spec$data_source,
      metadata = analysis_obj$spec$analysis_data$metadata,
      feature_matrix = analysis_obj$spec$analysis_data$feature_matrix,
      exposures = analysis_obj$exposures,
      covariates = analysis_obj$covariates,
      selected_pcs = selected_pcs,
      pc_signs = unname(pc_signs[selected_pcs])
    ),
    algo = "xxhash64"
  )
}

bootstrap_frequentist_dynamic <- function(frequentist_result,
                                          B = 1000,
                                          seed = 123,
                                          rho_grid = seq(-0.3, 0.3, by = 0.05),
                                          progress_callback = NULL) {
  prepared <- frequentist_result$prepared
  set.seed(seed)

  n <- nrow(prepared$data)
  boot_summary <- vector("list", B)
  boot_path <- vector("list", B)
  sensitivity_components_list <- vector("list", B)

  sensitivity_components <- function(dat) {
    cov_part <- if (!is.null(prepared$covariate_terms) && length(prepared$covariate_terms) > 0) {
      paste(prepared$covariate_terms, collapse = " + ")
    } else {
      NULL
    }
    exposure_part <- paste(prepared$exposures, collapse = " + ")
    mediator_part <- paste(prepared$mediator_model_names, collapse = " + ")

    fit_te <- stats::lm(
      stats::as.formula(paste("Y ~", exposure_part, if (!is.null(cov_part)) paste("+", cov_part) else "")),
      data = dat
    )
    te <- stats::coef(fit_te)[prepared$exposures]

    mediator_fits <- lapply(prepared$mediator_model_names, function(med) {
      stats::lm(
        stats::as.formula(paste(med, "~", exposure_part, if (!is.null(cov_part)) paste("+", cov_part) else "")),
        data = dat
      )
    })
    names(mediator_fits) <- prepared$mediator_model_names

    alpha_mat <- sapply(prepared$mediator_model_names, function(med) {
      stats::coef(mediator_fits[[med]])[prepared$exposures]
    })
    alpha_mat <- matrix(
      alpha_mat,
      nrow = length(prepared$exposures),
      ncol = length(prepared$mediator_model_names),
      dimnames = list(prepared$exposures, prepared$mediator_model_names)
    )

    mediator_resid <- sapply(prepared$mediator_model_names, function(med) stats::resid(mediator_fits[[med]]))
    mediator_resid <- matrix(
      mediator_resid,
      ncol = length(prepared$mediator_model_names),
      dimnames = list(NULL, prepared$mediator_model_names)
    )
    sigma_m <- stats::cov(mediator_resid)
    sigma_m <- matrix(
      sigma_m,
      nrow = length(prepared$mediator_model_names),
      ncol = length(prepared$mediator_model_names),
      dimnames = list(prepared$mediator_model_names, prepared$mediator_model_names)
    )
    sd_m <- stats::setNames(apply(mediator_resid, 2, stats::sd), prepared$mediator_model_names)

    fit_y <- stats::lm(
      stats::as.formula(paste("Y ~", exposure_part, "+", mediator_part, if (!is.null(cov_part)) paste("+", cov_part) else "")),
      data = dat
    )
    beta_hat <- stats::coef(fit_y)[prepared$mediator_model_names]
    nde_hat <- stats::coef(fit_y)[prepared$exposures]
    sd_y_resid <- stats::sd(stats::resid(fit_y))
    nie_hat <- as.vector(alpha_mat %*% beta_hat)
    names(nie_hat) <- prepared$exposures

    observed_summary <- tibble::tibble(
      exposure = prepared$exposures,
      TE = te[prepared$exposures],
      NDE = te[prepared$exposures] - nie_hat[prepared$exposures],
      NIE = nie_hat[prepared$exposures],
      PM = nie_hat[prepared$exposures] / te[prepared$exposures],
      outcome_direct_coef = nde_hat[prepared$exposures]
    )

    observed_path <- tibble::as_tibble(alpha_mat, rownames = "exposure") |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(prepared$mediator_model_names),
        names_to = "mediator",
        values_to = "alpha_X_to_M"
      ) |>
      dplyr::left_join(
        tibble::tibble(mediator = prepared$mediator_model_names, beta_M_to_Y = beta_hat[prepared$mediator_model_names]),
        by = "mediator"
      ) |>
      dplyr::mutate(indirect_component = alpha_X_to_M * beta_M_to_Y)

    list(
      observed_summary = observed_summary,
      observed_path = observed_path,
      alpha_mat = alpha_mat,
      beta_hat = beta_hat,
      te = te,
      sigma_m = sigma_m,
      sd_m = sd_m,
      sd_y_resid = sd_y_resid
    )
  }

  for (b in seq_len(B)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    dat_b <- prepared$data[idx, , drop = FALSE]

    comp <- sensitivity_components(dat_b)
    boot_summary[[b]] <- comp$observed_summary
    boot_path[[b]] <- comp$observed_path
    sensitivity_components_list[[b]] <- comp[c("alpha_mat", "beta_hat", "te", "sigma_m", "sd_m", "sd_y_resid")]

    if (is.function(progress_callback)) {
      progress_callback(b, B)
    }
  }

  boot <- list(
    boot_summary = boot_summary,
    boot_path = boot_path,
    boot_summary_df = dplyr::bind_rows(boot_summary, .id = "boot_id"),
    boot_path_df = dplyr::bind_rows(boot_path, .id = "boot_id")
  )
  ci <- summarize_bootstrap_ci(boot$boot_summary_df, boot$boot_path_df)
  ci$summary_ci <- ci$summary_ci |>
    dplyr::mutate(contrast = unname(prepared$exposure_labels[exposure]))
  ci$summary_ci$contrast <- ifelse(is.na(ci$summary_ci$contrast), ci$summary_ci$exposure, ci$summary_ci$contrast)
  ci$path_ci <- ci$path_ci |>
    map_frequentist_mediator_names(prepared$mediator_map) |>
    dplyr::mutate(contrast = unname(prepared$exposure_labels[exposure]))
  ci$path_ci$contrast <- ifelse(is.na(ci$path_ci$contrast), ci$path_ci$exposure, ci$path_ci$contrast)

  sensitivity_bundle <- list(
    components = sensitivity_components_list,
    exposures = prepared$exposures,
    mediators = prepared$mediator_model_names,
    mediator_map = prepared$mediator_map,
    exposure_labels = prepared$exposure_labels
  )

  list(
    bootstrap = boot,
    ci = ci,
    sensitivity_bundle = sensitivity_bundle
  )
}

frequentist_sensitivity_bootstrap_from_bundle <- function(bundle,
                                                          rho_grid = seq(-0.3, 0.3, by = 0.05),
                                                          feasibility_tol = 1e-10) {
  exposures <- bundle$exposures
  mediators <- bundle$mediators
  exposure_labels <- bundle$exposure_labels

  map_mediator_names <- function(df) {
    map_frequentist_mediator_names(df, bundle$mediator_map)
  }

  make_result_for_component <- function(comp, boot_id, rho_vec, rho, scenario, varied_mediator) {
    adjustment <- structural_rho_adjustment(comp$sigma_m, comp$sd_y_resid, rho_vec, feasibility_tol)
    cov_my <- stats::setNames(as.numeric(adjustment$covariance), mediators)
    feasibility_ratio <- adjustment$q[[1]]
    valid_parameter <- adjustment$valid[[1]]

    beta_adj <- if (valid_parameter) {
      as.vector(comp$beta_hat - adjustment$adjustment[1, ])
    } else {
      rep(NA_real_, length(mediators))
    }
    names(beta_adj) <- mediators

    nie_adj <- if (valid_parameter) {
      as.vector(comp$alpha_mat %*% beta_adj)
    } else {
      rep(NA_real_, length(exposures))
    }
    names(nie_adj) <- exposures

    summary <- tibble::tibble(
      boot_id = boot_id,
      scenario = scenario,
      varied_mediator = varied_mediator,
      exposure = exposures,
      contrast = unname(exposure_labels[exposures]),
      rho = rho,
      TE = comp$te[exposures],
      NIE_adjusted = nie_adj[exposures],
      NDE_adjusted = comp$te[exposures] - nie_adj[exposures],
      PM_adjusted = nie_adj[exposures] / comp$te[exposures],
      feasibility_ratio = feasibility_ratio,
      valid_parameter = valid_parameter
    )
    summary$contrast <- ifelse(is.na(summary$contrast), summary$exposure, summary$contrast)

    path <- tibble::as_tibble(comp$alpha_mat, rownames = "exposure") |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(mediators),
        names_to = "mediator",
        values_to = "alpha_X_to_M"
      ) |>
      dplyr::left_join(
        tibble::tibble(
          mediator = mediators,
          beta_observed = comp$beta_hat[mediators],
          beta_adjusted = beta_adj[mediators],
          residual_covariance = cov_my[mediators]
        ),
        by = "mediator"
      ) |>
      dplyr::mutate(
        boot_id = boot_id,
        scenario = scenario,
        varied_mediator = varied_mediator,
        contrast = unname(exposure_labels[exposure]),
        rho = rho,
        indirect_observed = alpha_X_to_M * beta_observed,
        indirect_adjusted = alpha_X_to_M * beta_adjusted,
        feasibility_ratio = feasibility_ratio,
        valid_parameter = valid_parameter
      )
    path$contrast <- ifelse(is.na(path$contrast), path$exposure, path$contrast)

    list(
      summary = summary,
      path = path
    )
  }

  common_summary <- list()
  common_path <- list()
  one_summary <- list()
  one_path <- list()
  common_i <- 1
  one_i <- 1

  for (b in seq_along(bundle$components)) {
    comp <- bundle$components[[b]]
    boot_id <- as.character(b)

    for (rho in rho_grid) {
      rho_vec <- stats::setNames(rep(rho, length(mediators)), mediators)
      common <- make_result_for_component(
        comp = comp,
        boot_id = boot_id,
        rho_vec = rho_vec,
        rho = rho,
        scenario = "common",
        varied_mediator = "all_mediators"
      )
      common_summary[[common_i]] <- common$summary
      common_path[[common_i]] <- common$path
      common_i <- common_i + 1
    }

    for (med in mediators) {
      for (rho in rho_grid) {
        rho_vec <- stats::setNames(rep(0, length(mediators)), mediators)
        rho_vec[med] <- rho
        one <- make_result_for_component(
          comp = comp,
          boot_id = boot_id,
          rho_vec = rho_vec,
          rho = rho,
          scenario = "one_at_a_time",
          varied_mediator = med
        )
        one_summary[[one_i]] <- one$summary
        one_path[[one_i]] <- one$path
        one_i <- one_i + 1
      }
    }
  }

  list(
    common_summary_df = dplyr::bind_rows(common_summary),
    common_path_df = dplyr::bind_rows(common_path) |> map_mediator_names(),
    one_summary_df = dplyr::bind_rows(one_summary) |>
      map_frequentist_mediator_names(bundle$mediator_map, columns = "varied_mediator"),
    one_path_df = dplyr::bind_rows(one_path) |> map_mediator_names()
  )
}

summarize_frequentist_sensitivity_bundle_ci <- function(bundle,
                                                        rho_grid = seq(-0.3, 0.3, by = 0.05),
                                                        feasibility_tol = 1e-10) {
  components <- bundle$components
  exposures <- bundle$exposures
  mediators <- bundle$mediators
  exposure_labels <- bundle$exposure_labels
  mediator_name_map <- stats::setNames(bundle$mediator_map$mediator, bundle$mediator_map$model_name)
  mediator_labels <- unname(mediator_name_map[mediators])
  mediator_labels <- ifelse(is.na(mediator_labels), mediators, mediator_labels)

  exposure_label <- function(x) {
    out <- unname(exposure_labels[x])
    ifelse(is.na(out), x, out)
  }
  q025 <- function(x) stats::quantile(x, 0.025, na.rm = TRUE)
  q975 <- function(x) stats::quantile(x, 0.975, na.rm = TRUE)
  summarize_vec <- function(x) {
    c(mean = mean(x, na.rm = TRUE), lwr = unname(q025(x)), upr = unname(q975(x)))
  }

  common_summary_rows <- list()
  common_path_rows <- list()
  one_summary_rows <- list()
  one_path_rows <- list()
  common_summary_i <- 1
  common_path_i <- 1
  one_summary_i <- 1
  one_path_i <- 1

  for (rho in rho_grid) {
    common_direct <- matrix(NA_real_, nrow = length(components), ncol = length(exposures))
    common_indirect <- matrix(NA_real_, nrow = length(components), ncol = length(exposures))
    common_total <- matrix(NA_real_, nrow = length(components), ncol = length(exposures))
    common_alpha <- array(NA_real_, dim = c(length(components), length(exposures), length(mediators)))
    common_beta <- matrix(NA_real_, nrow = length(components), ncol = length(mediators))
    common_path_indirect <- array(NA_real_, dim = c(length(components), length(exposures), length(mediators)))

    for (b in seq_along(components)) {
      comp <- components[[b]]
      adjustment <- structural_rho_adjustment(comp$sigma_m, comp$sd_y_resid,
                                              rep(rho, length(mediators)), feasibility_tol)
      solved <- as.numeric(adjustment$adjustment)
      valid <- adjustment$valid[[1]]
      if (valid) {
        beta_adj <- as.numeric(comp$beta_hat - solved)
        names(beta_adj) <- mediators
        indirect_components <- sweep(comp$alpha_mat, 2, beta_adj, `*`)
        indirect_total <- rowSums(indirect_components)
        common_direct[b, ] <- comp$te[exposures] - indirect_total[exposures]
        common_indirect[b, ] <- indirect_total[exposures]
        common_total[b, ] <- comp$te[exposures]
        common_alpha[b, , ] <- comp$alpha_mat[exposures, mediators, drop = FALSE]
        common_beta[b, ] <- beta_adj[mediators]
        common_path_indirect[b, , ] <- indirect_components[exposures, mediators, drop = FALSE]
      }
    }

    for (e in seq_along(exposures)) {
      d <- summarize_vec(common_direct[, e])
      ind <- summarize_vec(common_indirect[, e])
      tot <- summarize_vec(common_total[, e])
      common_summary_rows[[common_summary_i]] <- tibble::tibble(
        scenario = "common",
        varied_mediator = "all_mediators",
        exposure = exposures[[e]],
        contrast = exposure_label(exposures[[e]]),
        rho = rho,
        adjusted_direct_mean = d[["mean"]],
        adjusted_direct_lwr = d[["lwr"]],
        adjusted_direct_upr = d[["upr"]],
        adjusted_indirect_mean = ind[["mean"]],
        adjusted_indirect_lwr = ind[["lwr"]],
        adjusted_indirect_upr = ind[["upr"]],
        adjusted_total_mean = tot[["mean"]],
        adjusted_total_lwr = tot[["lwr"]],
        adjusted_total_upr = tot[["upr"]]
      )
      common_summary_i <- common_summary_i + 1

      for (m in seq_along(mediators)) {
        a <- summarize_vec(common_alpha[, e, m])
        beta <- summarize_vec(common_beta[, m])
        indirect <- summarize_vec(common_path_indirect[, e, m])
        common_path_rows[[common_path_i]] <- tibble::tibble(
          scenario = "common",
          varied_mediator = "all_mediators",
          exposure = exposures[[e]],
          contrast = exposure_label(exposures[[e]]),
          mediator = mediator_labels[[m]],
          rho = rho,
          alpha_mean = a[["mean"]],
          alpha_lwr = a[["lwr"]],
          alpha_upr = a[["upr"]],
          beta_mean = beta[["mean"]],
          beta_lwr = beta[["lwr"]],
          beta_upr = beta[["upr"]],
          indirect_mean = indirect[["mean"]],
          indirect_lwr = indirect[["lwr"]],
          indirect_upr = indirect[["upr"]]
        )
        common_path_i <- common_path_i + 1
      }
    }

    for (v in seq_along(mediators)) {
      one_direct <- matrix(NA_real_, nrow = length(components), ncol = length(exposures))
      one_indirect <- matrix(NA_real_, nrow = length(components), ncol = length(exposures))
      one_total <- matrix(NA_real_, nrow = length(components), ncol = length(exposures))
      one_alpha <- array(NA_real_, dim = c(length(components), length(exposures), length(mediators)))
      one_beta <- matrix(NA_real_, nrow = length(components), ncol = length(mediators))
      one_path_indirect <- array(NA_real_, dim = c(length(components), length(exposures), length(mediators)))

      for (b in seq_along(components)) {
        comp <- components[[b]]
        rho_vec <- stats::setNames(rep(0, length(mediators)), mediators)
        rho_vec[[v]] <- rho
        adjustment <- structural_rho_adjustment(comp$sigma_m, comp$sd_y_resid, rho_vec, feasibility_tol)
        solved <- as.numeric(adjustment$adjustment)
        valid <- adjustment$valid[[1]]
        if (valid) {
          beta_adj <- as.numeric(comp$beta_hat - solved)
          names(beta_adj) <- mediators
          indirect_components <- sweep(comp$alpha_mat, 2, beta_adj, `*`)
          indirect_total <- rowSums(indirect_components)
          one_direct[b, ] <- comp$te[exposures] - indirect_total[exposures]
          one_indirect[b, ] <- indirect_total[exposures]
          one_total[b, ] <- comp$te[exposures]
          one_alpha[b, , ] <- comp$alpha_mat[exposures, mediators, drop = FALSE]
          one_beta[b, ] <- beta_adj[mediators]
          one_path_indirect[b, , ] <- indirect_components[exposures, mediators, drop = FALSE]
        }
      }

      for (e in seq_along(exposures)) {
        d <- summarize_vec(one_direct[, e])
        ind <- summarize_vec(one_indirect[, e])
        tot <- summarize_vec(one_total[, e])
        one_summary_rows[[one_summary_i]] <- tibble::tibble(
          scenario = "one_at_a_time",
          varied_mediator = mediator_labels[[v]],
          exposure = exposures[[e]],
          contrast = exposure_label(exposures[[e]]),
          rho = rho,
          adjusted_direct_mean = d[["mean"]],
          adjusted_direct_lwr = d[["lwr"]],
          adjusted_direct_upr = d[["upr"]],
          adjusted_indirect_mean = ind[["mean"]],
          adjusted_indirect_lwr = ind[["lwr"]],
          adjusted_indirect_upr = ind[["upr"]],
          adjusted_total_mean = tot[["mean"]],
          adjusted_total_lwr = tot[["lwr"]],
          adjusted_total_upr = tot[["upr"]]
        )
        one_summary_i <- one_summary_i + 1

        for (m in seq_along(mediators)) {
          a <- summarize_vec(one_alpha[, e, m])
          beta <- summarize_vec(one_beta[, m])
          indirect <- summarize_vec(one_path_indirect[, e, m])
          one_path_rows[[one_path_i]] <- tibble::tibble(
            scenario = "one_at_a_time",
            varied_mediator = mediator_labels[[v]],
            exposure = exposures[[e]],
            contrast = exposure_label(exposures[[e]]),
            mediator = mediator_labels[[m]],
            rho = rho,
            alpha_mean = a[["mean"]],
            alpha_lwr = a[["lwr"]],
            alpha_upr = a[["upr"]],
            beta_mean = beta[["mean"]],
            beta_lwr = beta[["lwr"]],
            beta_upr = beta[["upr"]],
            indirect_mean = indirect[["mean"]],
            indirect_lwr = indirect[["lwr"]],
            indirect_upr = indirect[["upr"]]
          )
          one_path_i <- one_path_i + 1
        }
      }
    }
  }

  list(
    common_summary_ci = dplyr::bind_rows(common_summary_rows),
    common_path_ci = dplyr::bind_rows(common_path_rows),
    one_summary_ci = dplyr::bind_rows(one_summary_rows),
    one_path_ci = dplyr::bind_rows(one_path_rows)
  )
}

summarize_frequentist_sensitivity_bootstrap_ci <- function(boot_sensitivity,
                                                           rho_grid = seq(-0.3, 0.3, by = 0.05)) {
  if ("components" %in% names(boot_sensitivity)) {
    return(summarize_frequentist_sensitivity_bundle_ci(
      bundle = boot_sensitivity,
      rho_grid = rho_grid
    ))
  }

  summarize_summary_ci <- function(df) {
    df |>
      dplyr::filter(valid_parameter) |>
      dplyr::mutate(total_adjusted = NDE_adjusted + NIE_adjusted) |>
      dplyr::group_by(scenario, varied_mediator, exposure, contrast, rho) |>
      dplyr::summarise(
        adjusted_direct_mean = mean(NDE_adjusted, na.rm = TRUE),
        adjusted_direct_lwr = stats::quantile(NDE_adjusted, 0.025, na.rm = TRUE),
        adjusted_direct_upr = stats::quantile(NDE_adjusted, 0.975, na.rm = TRUE),
        adjusted_indirect_mean = mean(NIE_adjusted, na.rm = TRUE),
        adjusted_indirect_lwr = stats::quantile(NIE_adjusted, 0.025, na.rm = TRUE),
        adjusted_indirect_upr = stats::quantile(NIE_adjusted, 0.975, na.rm = TRUE),
        adjusted_total_mean = mean(total_adjusted, na.rm = TRUE),
        adjusted_total_lwr = stats::quantile(total_adjusted, 0.025, na.rm = TRUE),
        adjusted_total_upr = stats::quantile(total_adjusted, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
  }

  summarize_path_ci <- function(df) {
    df |>
      dplyr::filter(valid_parameter) |>
      dplyr::group_by(scenario, varied_mediator, exposure, contrast, mediator, rho) |>
      dplyr::summarise(
        alpha_mean = mean(alpha_X_to_M, na.rm = TRUE),
        alpha_lwr = stats::quantile(alpha_X_to_M, 0.025, na.rm = TRUE),
        alpha_upr = stats::quantile(alpha_X_to_M, 0.975, na.rm = TRUE),
        beta_mean = mean(beta_adjusted, na.rm = TRUE),
        beta_lwr = stats::quantile(beta_adjusted, 0.025, na.rm = TRUE),
        beta_upr = stats::quantile(beta_adjusted, 0.975, na.rm = TRUE),
        indirect_mean = mean(indirect_adjusted, na.rm = TRUE),
        indirect_lwr = stats::quantile(indirect_adjusted, 0.025, na.rm = TRUE),
        indirect_upr = stats::quantile(indirect_adjusted, 0.975, na.rm = TRUE),
        .groups = "drop"
      )
  }

  list(
    common_summary_ci = summarize_summary_ci(boot_sensitivity$common_summary_df),
    common_path_ci = summarize_path_ci(boot_sensitivity$common_path_df),
    one_summary_ci = summarize_summary_ci(boot_sensitivity$one_summary_df),
    one_path_ci = summarize_path_ci(boot_sensitivity$one_path_df)
  )
}

compute_frequentist_dynamic_sensitivity <- function(frequentist_result,
                                                    rho_grid = seq(-0.3, 0.3, by = 0.05),
                                                    progress_callback = NULL) {
  prepared <- frequentist_result$prepared
  notify <- function(value, detail) {
    if (is.function(progress_callback)) {
      progress_callback(value, detail)
    }
  }

  notify(0.10, "Computing common residual-correlation curve")
  common <- residual_corr_sensitivity(
    data = prepared$data,
    exposures = prepared$exposures,
    outcome = "Y",
    mediators = prepared$mediator_model_names,
    covariates = prepared$covariate_terms,
    rho_grid = rho_grid,
    mode = "common"
  )
  notify(0.45, "Computing one-mediator-at-a-time curves")
  one_at_a_time <- residual_corr_sensitivity(
    data = prepared$data,
    exposures = prepared$exposures,
    outcome = "Y",
    mediators = prepared$mediator_model_names,
    covariates = prepared$covariate_terms,
    rho_grid = rho_grid,
    mode = "one_at_a_time"
  )

  notify(0.70, "Checking rho = 0 baseline agreement")
  rho0_common <- validate_rho0_matches_observed(common)
  rho0_one_at_a_time <- validate_rho0_matches_observed(one_at_a_time)

  map_result <- function(x) {
    x$observed_path <- x$observed_path |>
      map_frequentist_mediator_names(prepared$mediator_map)
    x$sensitivity_summary <- x$sensitivity_summary |>
      map_frequentist_mediator_names(prepared$mediator_map, columns = "varied_mediator")
    x$sensitivity_path <- x$sensitivity_path |>
      map_frequentist_mediator_names(prepared$mediator_map)
    x$sensitivity_beta <- x$sensitivity_beta |>
      map_frequentist_mediator_names(prepared$mediator_map)
    x
  }

  notify(0.82, "Preparing mediator labels")
  common <- map_result(common)
  one_at_a_time <- map_result(one_at_a_time)

  notify(0.92, "Summarizing zero crossings")
  list(
    common = common,
    one_at_a_time = one_at_a_time,
    rho0_common = rho0_common,
    rho0_one_at_a_time = rho0_one_at_a_time,
    tipping_common = summarize_rho_tipping_points(common),
    tipping_one_at_a_time = summarize_rho_tipping_points(one_at_a_time)
  )
}

validate_frequentist_sign_reversal <- function(analysis_obj, selected_pcs, pc_signs, pc_to_reverse) {
  if (!pc_to_reverse %in% selected_pcs) {
    stop("pc_to_reverse must be one of the selected PCs.", call. = FALSE)
  }
  base <- make_frequentist_dynamic_fit(analysis_obj, selected_pcs, pc_signs)
  flipped_signs <- pc_signs
  flipped_signs[[pc_to_reverse]] <- -flipped_signs[[pc_to_reverse]]
  flipped <- make_frequentist_dynamic_fit(analysis_obj, selected_pcs, flipped_signs)

  base_path <- base$fit$path
  flipped_path <- flipped$fit$path
  joined_path <- dplyr::left_join(
    base_path,
    flipped_path,
    by = c("exposure", "contrast", "mediator"),
    suffix = c("_base", "_flipped")
  )
  joined_summary <- dplyr::left_join(
    base$fit$summary,
    flipped$fit$summary,
    by = c("exposure", "contrast"),
    suffix = c("_base", "_flipped")
  )

  reversed_rows <- joined_path$mediator == pc_to_reverse
  list(
    alpha_flips = all(abs(joined_path$alpha_X_to_M_base[reversed_rows] + joined_path$alpha_X_to_M_flipped[reversed_rows]) < 1e-8),
    beta_flips = all(abs(joined_path$beta_M_to_Y_base[reversed_rows] + joined_path$beta_M_to_Y_flipped[reversed_rows]) < 1e-8),
    indirect_invariant = all(abs(joined_path$indirect_component_base - joined_path$indirect_component_flipped) < 1e-8),
    decomposition_invariant = all(abs(joined_summary$TE_base - joined_summary$TE_flipped) < 1e-8) &&
      all(abs(joined_summary$NDE_base - joined_summary$NDE_flipped) < 1e-8) &&
      all(abs(joined_summary$NIE_base - joined_summary$NIE_flipped) < 1e-8) &&
      all(abs(joined_summary$PM_base - joined_summary$PM_flipped) < 1e-8)
  )
}
