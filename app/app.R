suppressPackageStartupMessages({
  library(shiny)
  library(DT)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

find_project_root <- function() {
  candidates <- unique(normalizePath(
    c(getwd(), file.path(getwd(), "..")),
    winslash = "/",
    mustWork = FALSE
  ))

	  for (candidate in candidates) {
	    if (file.exists(file.path(candidate, "R", "pca_mediation_pipeline.R")) &&
	        file.exists(file.path(candidate, "R", "data_input_app.R")) &&
	        file.exists(file.path(candidate, "R", "bayesian_mediation_v1.R")) &&
	        file.exists(file.path(candidate, "R", "bayesian_sensitivity_A_only.R")) &&
	        file.exists(file.path(candidate, "R", "analysis_bayesian_app.R")) &&
	        file.exists(file.path(candidate, "R", "bayesian_dynamic_app.R")) &&
	        file.exists(file.path(candidate, "R", "frequentist_dynamic_app.R")) &&
	        file.exists(file.path(candidate, "R", "sensitivity_v2.R")) &&
	        file.exists(file.path(candidate, "data", "raw", "pt16_emory_GEX_immune_FULL_v2.rds")) &&
	        file.exists(file.path(candidate, "results", "bayesian", "current", "emory_bayesian_v1_sensitivity_A_app_artifact.rds"))) {
	      return(candidate)
	    }
	  }

  stop("Could not locate the project root from the Shiny app directory.")
}

project_root <- find_project_root()
source(file.path(project_root, "R", "pca_mediation_pipeline.R"))
source(file.path(project_root, "R", "data_input_app.R"))
source(file.path(project_root, "R", "bayesian_mediation_v1.R"))
source(file.path(project_root, "R", "bayesian_sensitivity_A_only.R"))
source(file.path(project_root, "R", "analysis_bayesian_app.R"))
source(file.path(project_root, "R", "bayesian_dynamic_app.R"))
source(file.path(project_root, "R", "frequentist_dynamic_app.R"))
source(file.path(project_root, "R", "sensitivity_v2.R"))
source(file.path(project_root, "R", "app_result_export.R"))

rds_path <- file.path(project_root, "data", "raw", "pt16_emory_GEX_immune_FULL_v2.rds")
bayesian_artifact_path <- file.path(project_root, "results", "bayesian", "current",
  "emory_bayesian_v1_sensitivity_A_app_artifact.rds")
bayesian_app_artifact <- load_bayesian_app_artifact(bayesian_artifact_path)
bayesian_app_artifact_pc2_dynamic_path <- file.path(
  project_root,
  "results",
  "bayesian",
  "current",
  "emory_bayesian_dynamic_PC1_PC2_default_app_artifact.rds"
)
bayesian_app_artifact_pc2_dynamic <- if (file.exists(bayesian_app_artifact_pc2_dynamic_path)) {
  load_bayesian_app_artifact(bayesian_app_artifact_pc2_dynamic_path)
} else {
  NULL
}
bayesian_validated_artifacts <- Filter(
  function(x) !is.null(x$artifact),
  list(
    PC1_PC2 = list(pc_ids = c("PC1", "PC2"), pc_signs = c(PC1 = 1, PC2 = 1),
                  artifact = bayesian_app_artifact_pc2_dynamic),
    PC1_PC2_PC3 = list(pc_ids = c("PC1", "PC2", "PC3"),
                      pc_signs = c(PC1 = -1, PC2 = -1, PC3 = 1),
                      artifact = bayesian_app_artifact)
  )
)
for (entry in bayesian_validated_artifacts) {
  orientation <- entry$artifact$settings$pc_orientation
  artifact_mediators <- bayes_A_bundle_mediators(entry$artifact$sensitivity_A$residual_draws)
  # The publication artifact predates explicit orientation metadata.
  if (is.null(orientation) && identical(artifact_mediators, c("PC1_R", "PC2_R", "PC3"))) {
    orientation <- c(PC1 = -1, PC2 = -1, PC3 = 1)
  }
  if (is.null(orientation) || !setequal(names(orientation), entry$pc_ids) ||
      !isTRUE(all(orientation[entry$pc_ids] == entry$pc_signs))) {
    stop("Example artifact orientation does not match its selected-PC convention.")
  }
}
example_inputs <- example_data_inputs(rds_path)
example_spec <- build_analysis_spec(
  data_source = "example",
  main_data = example_inputs$main_data,
  feature_data = example_inputs$feature_data,
  roles = example_inputs$defaults,
  feature_id = example_inputs$defaults$feature_observation_id,
  excluded_feature_columns = example_inputs$defaults$excluded_feature_columns
)
example_analysis <- compute_analysis_from_spec(example_spec)
coord_x <- "imagecol"
coord_y <- "imagerow"
exposures <- c("X_adjacent", "X_far")
mediators <- c("PC1_R", "PC2_R", "PC3")
covariates <- c("x_coord", "y_coord")
rho_grid <- seq(-0.5, 0.5, by = 0.05)
bayesian_sensitivity_rho_grid <- seq(-0.3, 0.3, by = 0.05)
validated_pc_ids <- c("PC1", "PC2", "PC3")
pc_display_names <- c(PC1_R = "PC1", PC2_R = "PC2", PC3 = "PC3")
pc_mediator_ids <- c(PC1 = "PC1_R", PC2 = "PC2_R", PC3 = "PC3")

contrast_label <- function(x) {
  dplyr::case_when(
    x == "X_adjacent" ~ "adjacent vs inside",
    x == "X_far" ~ "far vs inside",
    TRUE ~ x
  )
}

round_numeric <- function(df, digits = 4) {
  df |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, digits)))
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

mediator_label <- function(x) {
  out <- unname(pc_display_names[x])
  out[x == "all"] <- "All mediators"
  out[grepl("^PC[0-9]+$", x)] <- x[grepl("^PC[0-9]+$", x)]
  ifelse(is.na(out), x, out)
}

mediator_filter_choices <- function(mediator_ids = mediators) {
  c("All mediators" = "all", stats::setNames(mediator_ids, mediator_label(mediator_ids)))
}

contrast_filter_choices <- function() {
  c("All contrasts" = "all", "adjacent vs inside" = "adjacent vs inside", "far vs inside" = "far vs inside")
}

pc_signs_for_selection <- function(selected, analysis_obj = NULL) {
  signs <- stats::setNames(rep(1, 6), paste0("PC", seq_len(6)))
  # This is an interpretation convention for the Emory publication example only.
  publication_example <- !is.null(analysis_obj) &&
    isTRUE(analysis_obj$bayesian_artifact_allowed) &&
    validate_example_bayesian_config(analysis_obj$spec) &&
    identical(analysis_obj$spec$analysis_data, example_analysis$spec$analysis_data)
  if (publication_example && setequal(selected, c("PC1", "PC2", "PC3"))) {
    signs[c("PC1", "PC2")] <- -1
  }
  signs
}

pc_sign_vector <- function(input, analysis_obj = NULL) {
  pc_signs_for_selection(selected_pc_ids(input, analysis_obj), analysis_obj)
}

mediation_sign_vector <- function(input) {
  pc_signs <- pc_sign_vector(input)
  c(
    PC1_R = pc_signs[["PC1"]],
    PC2_R = pc_signs[["PC2"]],
    PC3 = pc_signs[["PC3"]]
  )
}

pc_candidate_ids <- function(analysis_obj, max_pcs = 6) {
  score_cols <- ncol(analysis_obj$pca$pca_fit$x)
  feature_cols <- analysis_obj$spec$mediator_features$n_features %||% score_cols
  max_available <- min(max_pcs, score_cols, feature_cols)
  if (!is.finite(max_available) || max_available < 1) {
    return(character(0))
  }
  paste0("PC", seq_len(max_available))
}

selected_pc_ids <- function(input, analysis_obj = NULL) {
  candidates <- if (is.null(analysis_obj)) {
    paste0("PC", seq_len(6))
  } else {
    pc_candidate_ids(analysis_obj)
  }
  default_selected <- utils::head(candidates, 2)

  selected <- NULL
  pc_flags <- vapply(seq_along(candidates), function(i) {
    isTRUE(input[[paste0("pca_pc", i, "_selected")]])
  }, logical(1))
  if (any(!vapply(seq_along(candidates), function(i) {
    is.null(input[[paste0("pca_pc", i, "_selected")]])
  }, logical(1)))) {
    selected <- candidates[pc_flags]
  }
  selected <- selected %||% input$pca_mediation_pcs %||% default_selected
  if (!is.null(analysis_obj)) {
    selected <- intersect(selected, candidates)
  }
  selected
}

validated_selection_selected <- function(selected) {
  setequal(selected, validated_pc_ids)
}

bayesian_validated_artifact_for_selection <- function(selected) {
  selected <- sort(selected)
  for (key in names(bayesian_validated_artifacts)) {
    entry <- bayesian_validated_artifacts[[key]]
    if (setequal(selected, entry$pc_ids)) {
      return(list(key = key, pc_ids = entry$pc_ids, pc_signs = entry$pc_signs,
                  artifact = entry$artifact))
    }
  }
  NULL
}

orient_spatial_scores <- function(plot_df, signs, pca_fit = NULL) {
  out <- plot_df |>
    mutate(
      PC1 = PC1_R * signs[["PC1"]],
      PC2 = PC2_R * signs[["PC2"]],
      PC3 = PC3 * signs[["PC3"]]
    )

  if (!is.null(pca_fit)) {
    for (i in 4:min(6, ncol(pca_fit$x))) {
      pc <- paste0("PC", i)
      out[[pc]] <- pca_fit$x[, i] * signs[[pc]]
    }
  }

  out
}

orient_loading_data <- function(loading_df, signs, pca_fit = NULL) {
  out <- loading_df |>
    mutate(
      PC1 = PC1_R * signs[["PC1"]],
      PC2 = PC2_R * signs[["PC2"]],
      PC3 = PC3 * signs[["PC3"]]
    )

  if (!is.null(pca_fit)) {
    rotation <- pca_fit$rotation
    for (i in 4:min(6, ncol(rotation))) {
      pc <- paste0("PC", i)
      out[[pc]] <- rotation[, i] * signs[[pc]]
    }
  }

  out
}

orient_interval_columns <- function(df, prefix, signs) {
  sign_vec <- signs[df$mediator]
  old_lwr <- df[[paste0(prefix, "_lwr")]]
  old_upr <- df[[paste0(prefix, "_upr")]]
  reversed <- sign_vec < 0

  df[[paste0(prefix, "_est")]] <- df[[paste0(prefix, "_est")]] * sign_vec
  df[[paste0(prefix, "_lwr")]] <- ifelse(reversed, -old_upr, old_lwr)
  df[[paste0(prefix, "_upr")]] <- ifelse(reversed, -old_lwr, old_upr)
  df
}

orient_bayesian_path <- function(df, signs) {
  if (!all(df$mediator %in% names(signs))) {
    df$mediator_id <- df$mediator
    df$mediator <- mediator_label(df$mediator)
    return(df)
  }
  df <- orient_interval_columns(df, "alpha", signs)
  df <- orient_interval_columns(df, "beta", signs)
  df$mediator_id <- df$mediator
  df$mediator <- mediator_label(df$mediator)
  df
}

orient_frequentist_path <- function(df, signs) {
  sign_vec <- signs[df$mediator]
  df |>
    mutate(
      mediator_id = mediator,
      mediator = mediator_label(mediator),
      alpha = alpha * sign_vec,
      beta = beta * sign_vec
    )
}

filter_by_mediation_inputs <- function(df, contrast = "all", mediator = "all") {
  out <- df
  if (!identical(contrast, "all") && "contrast" %in% names(out)) {
    out <- dplyr::filter(out, contrast == !!contrast)
  }
  if (!identical(mediator, "all") && "mediator_id" %in% names(out)) {
    out <- dplyr::filter(out, mediator_id == !!mediator)
  }
  out
}

format_number <- function(x, digits = 3) {
  out <- rep(NA_character_, length(x))
  ok <- !is.na(x)
  out[ok] <- formatC(x[ok], format = "f", digits = digits)
  out
}

interval_excludes_zero <- function(lwr, upr) {
  if (is.null(lwr) || is.null(upr)) {
    return(FALSE)
  }
  ok <- is.finite(lwr) & is.finite(upr)
  ok & ((lwr > 0 & upr > 0) | (lwr < 0 & upr < 0))
}

interval_star <- function(lwr, upr) {
  ifelse(interval_excludes_zero(lwr, upr), "*", "")
}

interval_star_note <- function() {
  "A * marks estimates whose 95% interval does not include 0."
}

format_estimate_interval <- function(est, lwr = NULL, upr = NULL, digits = 3) {
  est_txt <- paste0(format_number(est, digits), interval_star(lwr, upr))
  if (is.null(lwr) || is.null(upr)) {
    return(est_txt)
  }

  paste0(est_txt, " [", format_number(lwr, digits), ", ", format_number(upr, digits), "]")
}

present_bayesian_mediation_summary <- function(df) {
  df |>
    transmute(
      Contrast = contrast,
      `Total Effect` = format_estimate_interval(TE_est, TE_lwr, TE_upr),
      `Direct Effect` = format_estimate_interval(direct_est, direct_lwr, direct_upr),
      `Indirect Effect (Total)` = format_estimate_interval(total_IIE_est, total_IIE_lwr, total_IIE_upr),
      `Proportion Mediated` = format_estimate_interval(PM_est, PM_lwr, PM_upr)
    )
}

present_bayesian_mediation_probabilities <- function(df) {
  out <- tibble::tibble(Contrast = df$contrast)
  if ("TE_Pr_gt_0" %in% names(df)) {
    out$`Total Effect` <- format_number(df$TE_Pr_gt_0)
  }
  if ("direct_Pr_gt_0" %in% names(df)) {
    out$`Direct Effect` <- format_number(df$direct_Pr_gt_0)
  }
  if ("total_IIE_Pr_gt_0" %in% names(df)) {
    out$`Indirect Effect` <- format_number(df$total_IIE_Pr_gt_0)
  }

  out
}

present_bayesian_mediation_path <- function(df) {
  out <- df |>
    transmute(
      Contrast = contrast,
      Mediator = mediator,
      `Exposure -> Mediator` = format_estimate_interval(alpha_est, alpha_lwr, alpha_upr),
      `Mediator -> Outcome` = format_estimate_interval(beta_est, beta_lwr, beta_upr),
      `Indirect Effect` = format_estimate_interval(IIE_est, IIE_lwr, IIE_upr)
    )

  if ("IIE_Pr_gt_0" %in% names(df)) {
    out$`Pr(Indirect Effect > 0)` <- format_number(df$IIE_Pr_gt_0)
  }
  if ("IIE_Pr_lt_0" %in% names(df)) {
    out$`Pr(Indirect Effect < 0)` <- format_number(df$IIE_Pr_lt_0)
  }

  out
}

present_frequentist_mediation_summary <- function(df) {
  df |>
    transmute(
      Contrast = contrast,
      `Total Effect` = format_estimate_interval(TE),
      `Direct Effect` = format_estimate_interval(direct_effect),
      `Indirect Effect (Total)` = format_estimate_interval(total_indirect_effect),
      `Proportion Mediated` = format_estimate_interval(PM)
    )
}

present_frequentist_mediation_path <- function(df) {
  df |>
    transmute(
      Contrast = contrast,
      Mediator = mediator,
      `Exposure -> Mediator` = format_estimate_interval(alpha, digits = 4),
      `Mediator -> Outcome` = format_estimate_interval(beta, digits = 4),
      `Indirect Effect` = format_estimate_interval(indirect_component, digits = 4)
    )
}

present_bootstrap_summary_ci <- function(df, contrast_labels = NULL) {
  label <- if ("contrast" %in% names(df)) {
    df$contrast
  } else if (is.null(contrast_labels)) {
    contrast_label(df$exposure)
  } else {
    out <- unname(contrast_labels[df$exposure])
    ifelse(is.na(out), df$exposure, out)
  }

  df |>
    mutate(contrast = label) |>
    transmute(
      Contrast = contrast,
      `Total Effect` = format_estimate_interval(TE_mean, TE_lwr, TE_upr),
      `Direct Effect` = format_estimate_interval(NDE_mean, NDE_lwr, NDE_upr),
      `Indirect Effect (Total)` = format_estimate_interval(NIE_mean, NIE_lwr, NIE_upr),
      `Proportion Mediated` = format_estimate_interval(PM_mean, PM_lwr, PM_upr)
    )
}

present_bootstrap_path_ci <- function(df, contrast_labels = NULL) {
  label <- if ("contrast" %in% names(df)) {
    df$contrast
  } else if (is.null(contrast_labels)) {
    contrast_label(df$exposure)
  } else {
    out <- unname(contrast_labels[df$exposure])
    ifelse(is.na(out), df$exposure, out)
  }

  df |>
    mutate(contrast = label) |>
    transmute(
      Contrast = contrast,
      Mediator = mediator_label(mediator),
      `Exposure -> Mediator` = format_estimate_interval(alpha_mean, alpha_lwr, alpha_upr),
      `Mediator -> Outcome` = format_estimate_interval(beta_mean, beta_lwr, beta_upr),
      `Indirect Effect` = format_estimate_interval(indirect_mean, indirect_lwr, indirect_upr)
    )
}

present_bayesian_sensitivity_summary <- function(df) {
  df |>
    transmute(
      Contrast = contrast,
      `Direct Effect` = format_estimate_interval(direct_est, direct_lwr, direct_upr),
      `Indirect Effect (Total)` = format_estimate_interval(total_IIE_est, total_IIE_lwr, total_IIE_upr),
      `Total Effect` = format_estimate_interval(TE_est, TE_lwr, TE_upr),
      `Proportion Mediated` = format_estimate_interval(PM_est, PM_lwr, PM_upr)
    )
}

present_bayesian_sensitivity_path <- function(df) {
  df |>
    transmute(
      Contrast = contrast,
      `ρ` = format_number(rho, digits = 2),
      `Mediator Varied` = mediator_label(mediator_varied),
      Mediator = mediator_label(mediator),
      `Indirect Effect` = format_estimate_interval(IIE_est, IIE_lwr, IIE_upr),
      `Pr(Indirect Effect > 0)` = format_number(IIE_Pr_gt_0),
      `Pr(Indirect Effect < 0)` = format_number(IIE_Pr_lt_0)
    )
}

present_frequentist_sensitivity_selected <- function(df) {
  for (column in c(
    "adjusted_direct_lwr", "adjusted_direct_upr",
    "adjusted_indirect_lwr", "adjusted_indirect_upr",
    "adjusted_total_lwr", "adjusted_total_upr"
  )) {
    if (!column %in% names(df)) {
      df[[column]] <- NA_real_
    }
  }
  if ("varied_mediator" %in% names(df)) {
    df |>
      transmute(
        Contrast = contrast,
        `Direct Effect` = format_estimate_interval(adjusted_direct_effect, adjusted_direct_lwr, adjusted_direct_upr),
        `Indirect Effect (Total)` = format_estimate_interval(adjusted_total_indirect, adjusted_indirect_lwr, adjusted_indirect_upr),
        `Total Effect` = format_estimate_interval(adjusted_total_effect, adjusted_total_lwr, adjusted_total_upr),
        `Proportion Mediated` = format_estimate_interval(PM_adjusted)
      )
  } else {
    df |>
      transmute(
        Contrast = contrast,
        `Direct Effect` = format_estimate_interval(adjusted_direct_effect, adjusted_direct_lwr, adjusted_direct_upr),
        `Indirect Effect (Total)` = format_estimate_interval(adjusted_total_indirect, adjusted_indirect_lwr, adjusted_indirect_upr),
        `Total Effect` = format_estimate_interval(adjusted_total_effect, adjusted_total_lwr, adjusted_total_upr),
        `Proportion Mediated` = format_estimate_interval(PM_adjusted)
      )
  }
}

present_feasibility_status <- function(summary_df) {
  fraction <- min(summary_df$valid_draw_fraction, na.rm = TRUE)
  percent <- round(100 * fraction)

  if (!is.finite(fraction) || fraction <= 0) {
    list(status = "Infeasible", text = "Infeasible for the selected ρ. No posterior draws support this setting.")
  } else if (fraction < 1) {
    list(
      status = "Partially feasible",
      text = paste0("Feasible for ", percent, "% of posterior draws. Interpret this ρ setting with caution.")
    )
  } else {
    list(status = "Feasible", text = "Feasible for 100% of posterior draws.")
  }
}

pm_warning_needed <- function(summary_df) {
  if (nrow(summary_df) == 0) {
    return(FALSE)
  }

  te_near_zero <- abs(summary_df$TE_est) < 0.01
  te_crosses_zero <- summary_df$TE_lwr <= 0 & summary_df$TE_upr >= 0
  opposing_signs <- summary_df$direct_est * summary_df$total_IIE_est < 0

  any(te_near_zero | te_crosses_zero | opposing_signs, na.rm = TRUE)
}

pm_warning_needed_frequentist <- function(summary_df) {
  if (nrow(summary_df) == 0) {
    return(FALSE)
  }
  te_near_zero <- abs(summary_df$adjusted_total_effect) < 0.01
  opposing_signs <- summary_df$adjusted_direct_effect * summary_df$adjusted_total_indirect < 0
  any(te_near_zero | opposing_signs, na.rm = TRUE)
}

spatial_plot_limits <- function(data, padding_fraction = 0.06) {
  x_range <- range(data$x_plot, na.rm = TRUE)
  y_range <- range(data$y_plot, na.rm = TRUE)
  x_pad <- diff(x_range) * padding_fraction
  y_pad <- diff(y_range) * padding_fraction
  if (!is.finite(x_pad) || x_pad == 0) x_pad <- 1
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 1

  list(
    x = c(x_range[[1]] - x_pad, x_range[[2]] + x_pad),
    y = c(y_range[[1]] - y_pad, y_range[[2]] + y_pad)
  )
}

compute_current_analysis <- function() {
  example_analysis
}

compute_sensitivity <- function(dat_full_allpc,
                                exposures = exposures,
                                mediators = mediators,
                                covariates = covariates) {
  common <- residual_corr_sensitivity(
    data = dat_full_allpc,
    exposures = exposures,
    outcome = "Y",
    mediators = mediators,
    covariates = covariates,
    rho_grid = rho_grid,
    mode = "common"
  )

  one_at_a_time <- residual_corr_sensitivity(
    data = dat_full_allpc,
    exposures = exposures,
    outcome = "Y",
    mediators = mediators,
    covariates = covariates,
    rho_grid = rho_grid,
    mode = "one_at_a_time"
  )

  list(
    common = common,
    one_at_a_time = one_at_a_time,
    rho0_common = validate_rho0_matches_observed(common),
    rho0_one_at_a_time = validate_rho0_matches_observed(one_at_a_time),
    tipping_common = summarize_rho_tipping_points(common),
    tipping_one_at_a_time = summarize_rho_tipping_points(one_at_a_time),
    plots = build_rho_sensitivity_plots(common, one_at_a_time)
  )
}

analysis_cache <- local({
  value <- example_analysis
  function() {
    value
  }
})

sensitivity_cache <- local({
  value <- NULL
  function(dat_full_allpc) {
    if (is.null(value)) {
      value <<- compute_sensitivity(dat_full_allpc)
    }
    value
  }
})

compute_bayesian_sensitivity_curve <- function(artifact, mode = c("common", "one_at_a_time"), mediator = NULL) {
  mode <- match.arg(mode)
  valid_mediators <- bayes_A_bundle_mediators(artifact$sensitivity_A$residual_draws)
  raw <- compute_bayesian_sensitivity_A_from_draws(
    residual_draws = artifact$sensitivity_A$residual_draws,
    rho_grid = bayesian_sensitivity_rho_grid,
    mode = mode,
    mediators = valid_mediators,
    mediators_varied = if (identical(mode, "one_at_a_time") && !is.null(mediator)) mediator else valid_mediators
  )

  summary <- canonicalize_bayesian_A_summary(raw)
  path <- canonicalize_bayesian_A_path(raw)
  feasibility <- raw$feasibility

  if (identical(mode, "one_at_a_time")) {
    mediator <- mediator %||% valid_mediators[[1]]
    summary <- dplyr::filter(summary, mediator_varied == mediator)
    path <- dplyr::filter(path, mediator_varied == mediator)
    feasibility <- dplyr::filter(feasibility, mediator_varied == mediator)
  }

  list(summary = summary, path = path, feasibility = feasibility)
}

bayesian_sensitivity_curve_cache <- local({
  values <- list()
  function(artifact, artifact_key, mode = c("common", "one_at_a_time"), mediator = NULL) {
    mode <- match.arg(mode)
    key <- paste(artifact_key, mode, mediator %||% "all", sep = ":")
    if (is.null(values[[key]])) {
      values[[key]] <<- compute_bayesian_sensitivity_curve(artifact = artifact, mode = mode, mediator = mediator)
    }
    values[[key]]
  }
})

bayesian_sensitivity_all_mediator_curve_cache <- local({
  values <- list()
  function(artifact, artifact_key) {
    key <- paste(artifact_key, "one_at_a_time", "all_mediators", sep = ":")
    if (is.null(values[[key]])) {
      raw <- compute_bayesian_sensitivity_A_from_draws(
        residual_draws = artifact$sensitivity_A$residual_draws,
        rho_grid = bayesian_sensitivity_rho_grid,
        mode = "one_at_a_time",
        mediators = bayes_A_bundle_mediators(artifact$sensitivity_A$residual_draws),
        mediators_varied = bayes_A_bundle_mediators(artifact$sensitivity_A$residual_draws)
      )
      values[[key]] <<- list(
        summary = canonicalize_bayesian_A_summary(raw),
        path = canonicalize_bayesian_A_path(raw),
        feasibility = raw$feasibility
      )
    }
    values[[key]]
  }
})

bayesian_sensitivity_selected_cache <- local({
  values <- list()
  function(artifact, artifact_key, rho, mode = c("common", "one_at_a_time"), mediator = NULL) {
    mode <- match.arg(mode)
    key <- paste(artifact_key, mode, mediator %||% "all", sprintf("%.12g", rho), sep = ":")
    if (is.null(values[[key]])) {
      values[[key]] <<- compute_bayesian_app_sensitivity_A(
        artifact = artifact,
        rho = rho,
        mode = mode,
        mediator = mediator
      )
    }
    values[[key]]
  }
})

top_loading_table <- function(loading_df, pc, direction = c("positive", "negative"), n = 30) {
  direction <- match.arg(direction)
  if (direction == "positive") {
    loading_df |>
      arrange(desc(.data[[pc]])) |>
      slice_head(n = n) |>
      select(gene, loading = all_of(pc))
  } else {
    loading_df |>
      arrange(.data[[pc]]) |>
      slice_head(n = n) |>
      select(gene, loading = all_of(pc))
  }
}

app_panel_card <- function(title, ..., class = "") {
  tags$section(
    class = paste("app-card", class),
    h4(title),
    ...
  )
}

app_help_tooltip <- function(text) {
  tags$span(
    class = "help-tooltip",
    `data-tooltip` = text,
    `aria-label` = text,
    tabindex = "0",
    "?"
  )
}

app_label_help <- function(label, tooltip, class = NULL) {
  tags$span(
    class = paste(c("concept-help-target", class), collapse = " "),
    `data-tooltip` = tooltip,
    `aria-label` = tooltip,
    tabindex = "0",
    label
  )
}

app_icon_help_target <- function(label, icon, tooltip) {
  tags$span(
    class = "help-tooltip setup-icon-help-target",
    `data-tooltip` = tooltip,
    `aria-label` = tooltip,
    tabindex = "0",
    tags$span(class = paste("glyphicon", icon), `aria-hidden` = "true"),
    tags$span(class = "setup-icon-help-label", label)
  )
}

app_icon_only_tooltip <- function(icon, tooltip) {
  tags$span(
    class = "help-tooltip setup-icon-only-tooltip",
    `data-tooltip` = tooltip,
    `aria-label` = tooltip,
    tabindex = "0",
    tags$span(class = paste("glyphicon", icon), `aria-hidden` = "true")
  )
}

setup_heading_icon <- function(label, icon, tooltip = NULL) {
  tagList(
    if (is.null(tooltip)) {
      tags$span(class = paste("glyphicon", icon, "setup-heading-icon"), `aria-hidden` = "true")
    } else {
      app_icon_only_tooltip(icon, tooltip)
    },
    tags$span(label)
  )
}

app_tutorial_section <- function(title, ..., open = FALSE) {
  args <- c(
    list(class = "tutorial-section"),
    if (isTRUE(open)) list(open = NA) else list(),
    list(tags$summary(title), tags$div(class = "tutorial-body", ...))
  )
  do.call(tags$details, args)
}

app_tutorial_step <- function(number, title, text) {
  tags$div(
    class = "tutorial-step",
    tags$div(class = "tutorial-step-number", number),
    tags$div(
      class = "tutorial-step-copy",
      tags$h4(title),
      tags$p(text)
    )
  )
}

app_tutorial_tile <- function(title, text) {
  tags$div(
    class = "tutorial-tile",
    tags$h4(title),
    tags$p(text)
  )
}

app_setup_section <- function(title, ...) {
  tags$section(
    class = "setup-sidebar-section",
    h4(title),
    ...
  )
}

app_kv_item <- function(label, value, class = "") {
  tags$div(
    class = paste("app-kv-item", class),
    tags$span(class = "app-kv-label", label),
    tags$span(class = "app-kv-value", value)
  )
}

setup_preview_table <- function(data, n_rows = 5, n_cols = 8) {
  preview <- table_preview(data, n_rows = n_rows, n_cols = n_cols)
  tags$div(
    class = "setup-static-preview",
    tags$table(
      tags$thead(
        tags$tr(lapply(names(preview), tags$th))
      ),
      tags$tbody(
        lapply(seq_len(nrow(preview)), function(i) {
          tags$tr(lapply(preview[i, , drop = FALSE], function(x) {
            value <- x[[1]]
            tags$td(if (is.na(value)) "" else as.character(value))
          }))
        })
      )
    )
  )
}

setup_column_letters <- function(n) {
  vapply(seq_len(n), function(i) {
    out <- character(0)
    while (i > 0) {
      i <- i - 1
      out <- c(LETTERS[(i %% 26) + 1], out)
      i <- i %/% 26
    }
    paste(out, collapse = "")
  }, character(1))
}

setup_spreadsheet_table <- function(data) {
  column_letters <- setup_column_letters(ncol(data))
  header_cells <- c(
    list(tags$th(class = "sheet-corner", "")),
    lapply(column_letters, function(x) tags$th(class = "sheet-col-letter", x))
  )
  name_row <- c(
    list(tags$th(class = "sheet-row-number", "1")),
    lapply(names(data), function(x) tags$td(class = "sheet-variable-name", x))
  )
  data_rows <- lapply(seq_len(nrow(data)), function(i) {
    c(
      list(tags$th(class = "sheet-row-number", as.character(i + 1))),
      lapply(data[i, , drop = FALSE], function(x) {
        value <- x[[1]]
        tags$td(if (is.na(value)) "" else as.character(value))
      })
    )
  })

  tags$div(
    class = "setup-spreadsheet-scroll",
    tags$table(
      class = "setup-spreadsheet-table",
      tags$thead(tags$tr(header_cells)),
      tags$tbody(c(list(tags$tr(name_row)), lapply(data_rows, tags$tr)))
    )
  )
}

setup_requirement_line <- function(label, text, tooltip = NULL) {
  tags$p(
    class = "setup-requirement-line",
    tags$strong(label),
    text,
    if (!is.null(tooltip)) app_help_tooltip(tooltip)
  )
}

setup_file_example <- function(title, description, columns) {
  tags$div(
    class = "setup-file-example",
    tags$div(class = "setup-file-icon", substr(title, 1, 1)),
    tags$div(
      tags$h5(title),
      tags$p(description),
      tags$div(
        class = "setup-column-chips",
        lapply(columns, function(column) tags$span(column))
      )
    )
  )
}

setup_spreadsheet_example <- function(title, data) {
  tags$div(
    class = "setup-spreadsheet-example",
    tags$h5(title),
    setup_preview_table(data, n_rows = nrow(data), n_cols = ncol(data))
  )
}

setup_match_cue <- function() {
  tags$div(
    class = "setup-match-cue",
    tags$div(
      tags$strong("Main dataset"),
      tags$span("spot_001"),
      tags$span("spot_002")
    ),
    tags$div(
      class = "setup-match-lines",
      tags$span("matches"),
      tags$span("matches")
    ),
    tags$div(
      tags$strong("Mediator feature matrix"),
      tags$span("spot_001"),
      tags$span("spot_002")
    )
  )
}

setup_metric_row <- function(label, value) {
  tags$div(
    class = "setup-metric-row",
    tags$span(label),
    tags$strong(value)
  )
}

pca_spec_label <- function(spec) {
  if (is.null(spec) || length(spec$selected_pcs) == 0) {
    return("Not confirmed")
  }
  paste(spec$selected_pcs, collapse = ", ")
}

ui <- navbarPage(
  title = "ASPEN",
  windowTitle = "ASPEN",
  id = "main_nav",
  footer = tags$footer(
    class = "app-brand-footer",
    tags$strong("ASPEN"),
    tags$span(`aria-hidden` = "true", HTML("&middot;")),
    tags$span("Analysis of Spatial Pathways and Effect Mediation")
  ),
  header = tagList(
    tags$head(
      tags$script(HTML(
        "
        Shiny.addCustomMessageHandler('setDisabled', function(x) {
          $('#' + x.id).prop('disabled', x.disabled);
        });
        (function() {
          var bayesTimer = null;
          var bayesStartTime = null;

          function bayesProgressState(elapsedSeconds) {
            if (elapsedSeconds < 12) {
              return {
                pct: 2 + elapsedSeconds * (16 / 12),
                label: 'Preparing Bayesian model inputs'
              };
            }
            if (elapsedSeconds < 90) {
              return {
                pct: 18 + (elapsedSeconds - 12) * (27 / 78),
                label: 'Fitting mediator model (stage 1 of 3)'
              };
            }
            if (elapsedSeconds < 170) {
              return {
                pct: 45 + (elapsedSeconds - 90) * (25 / 80),
                label: 'Fitting outcome model (stage 2 of 3)'
              };
            }
            if (elapsedSeconds < 260) {
              return {
                pct: 70 + (elapsedSeconds - 170) * (18 / 90),
                label: 'Fitting total-effect model (stage 3 of 3)'
              };
            }
            return {
              pct: Math.min(96, 88 + Math.log1p(elapsedSeconds - 260) * 2),
              label: 'Finalizing posterior summaries'
            };
          }

          function renderBayesProgress(pct, label) {
            var panel = $('#bayesian-client-progress');
            if (panel.length === 0) return;
            panel.removeClass('hidden');
            panel.find('.bayes-client-progress-fill').css('width', Math.max(0, Math.min(100, pct)) + '%');
            panel.find('.bayes-client-progress-percent').text(Math.round(pct) + '%');
            panel.find('.bayes-client-progress-label').text(label);
          }

          function startBayesProgress() {
            if (bayesTimer) {
              window.clearInterval(bayesTimer);
            }
            bayesStartTime = Date.now();
            renderBayesProgress(2, 'Starting Bayesian model fitting');
            bayesTimer = window.setInterval(function() {
              var elapsed = (Date.now() - bayesStartTime) / 1000;
              var state = bayesProgressState(elapsed);
              renderBayesProgress(state.pct, state.label);
            }, 1000);
          }

          function finishBayesProgress(status, message) {
            if (bayesTimer) {
              window.clearInterval(bayesTimer);
              bayesTimer = null;
            }
            if (status === 'error') {
              renderBayesProgress(100, message || 'Bayesian fit failed');
              $('#bayesian-client-progress').addClass('bayes-client-progress-error');
            } else {
              renderBayesProgress(100, message || 'Bayesian fit complete');
              $('#bayesian-client-progress').addClass('bayes-client-progress-complete');
            }
          }

          $(document).on('click', '#run_bayesian_mediation', function() {
            $('#bayesian-client-progress')
              .removeClass('bayes-client-progress-error bayes-client-progress-complete');
            startBayesProgress();
          });

          Shiny.addCustomMessageHandler('bayesianClientProgress', function(x) {
            if (x.action === 'start') {
              startBayesProgress();
            } else if (x.action === 'complete') {
              finishBayesProgress('complete', x.message);
            } else if (x.action === 'error') {
              finishBayesProgress('error', x.message);
            }
          });
        })();
        (function() {
          function openSpatialReference(trigger) {
            var wrapper = $(trigger).closest('.spatial-reference-popover');
            $('.spatial-reference-popover').not(wrapper).removeClass('is-open');
            wrapper.addClass('is-open');
          }

          $(document).on('mouseenter focus click', '.spatial-reference-trigger', function(event) {
            if (event.type === 'click') {
              event.preventDefault();
            }
            openSpatialReference(this);
          });

          $(document).on('click', '.spatial-reference-close', function(event) {
            event.preventDefault();
            event.stopPropagation();
            $(this).closest('.spatial-reference-popover').removeClass('is-open');
          });

          $(document).on('keydown', function(event) {
            if (event.key === 'Escape') {
              $('.spatial-reference-popover').removeClass('is-open');
            }
          });
        })();
        Shiny.addCustomMessageHandler('setNavWorkflow', function(x) {
          var navItems = $('a[data-value]');
          navItems
            .removeClass('nav-complete nav-available nav-unavailable nav-utility')
            .removeAttr('aria-disabled')
            .removeAttr('tabindex');
          navItems.find('.nav-check').remove();

          x.steps.forEach(function(step) {
            var anchor = $('a[data-value=\"' + step.value + '\"]').first();
            if (anchor.length === 0) return;
            if (step.hidden) {
              anchor.parent().hide();
              return;
            }
            anchor.parent().show();
            anchor.addClass(step.workflow ? 'nav-available' : 'nav-utility');
            if (step.completed) {
              anchor.addClass('nav-complete');
              anchor.prepend('<span class=\"nav-check\">&#10003;</span>');
            }
            if (!step.available) {
              anchor.addClass('nav-unavailable');
              anchor.attr('aria-disabled', 'true');
              anchor.attr('tabindex', '-1');
            }
          });
        });
        "
      )),
      tags$style(HTML("
	      :root {
	        --app-ink: #243447;
		        --app-accent: #274690;
		        --app-accent-2: #2f80ed;
		        --app-accent-soft: #eef4ff;
		        --app-blue: #2f80ed;
		        --app-blue-soft: #eaf2ff;
		        --app-violet: #7c5cc4;
		        --app-violet-soft: #f1edff;
		        --app-red: #c2413a;
		        --app-red-soft: #fdeceb;
		        --app-slate-accent: #64748b;
		        --app-slate-soft: #f1f5f9;
		        --app-border: #dde8ef;
		        --app-muted: #5f6f7f;
		        --app-soft: #f7fafc;
		      }
	      body {
	        background: #f7fafc;
	        color: var(--app-ink);
            min-height: 100vh;
            display: flex;
            flex-direction: column;
	      }
          body > .container-fluid {
            width: 100%;
            flex: 1 0 auto;
            display: flex;
            flex-direction: column;
          }
          body > .container-fluid > .tab-content {
            flex: 1 0 auto;
          }
	      .navbar {
	        margin-bottom: 0;
	        border-bottom: 1px solid var(--app-border);
		        box-shadow: 0 2px 8px rgba(39, 70, 144, 0.08);
	      }
	      .navbar-default {
	        background: #ffffff;
	        border-color: transparent transparent var(--app-border) transparent;
	      }
	      .navbar-default .navbar-brand {
	        color: var(--app-accent);
	        font-weight: 800;
	      }
	      .navbar-default .navbar-brand:hover,
	      .navbar-default .navbar-brand:focus {
		        color: #1d3470;
	      }
	      .navbar-default .navbar-nav > li > a {
	        color: var(--app-ink);
	        font-weight: 600;
	        border-bottom: 3px solid transparent;
	        transition: background-color 120ms ease, color 120ms ease, border-color 120ms ease;
	      }
	      .navbar-default .navbar-nav > li.active > a,
	      .navbar-default .navbar-nav > li.active > a:focus,
	      .navbar-default .navbar-nav > li.active > a:hover {
	        background: var(--app-accent-soft);
	        color: var(--app-accent);
	        border-bottom-color: var(--app-accent-2);
	      }
      .navbar-default .navbar-nav > li > a.nav-unavailable {
        color: #9aa7b3;
        cursor: default;
        pointer-events: none;
      }
      .navbar-default .navbar-nav > li > a.nav-unavailable:hover,
      .navbar-default .navbar-nav > li > a.nav-unavailable:focus {
        background: transparent;
        color: #9aa7b3;
      }
	      .nav-check {
		        color: var(--app-blue);
	        font-size: 11px;
	        margin-right: 5px;
        position: relative;
        top: -1px;
      }
      .app-brand-footer {
        display: flex;
        flex-wrap: wrap;
        align-items: baseline;
        gap: 4px 8px;
        margin-top: 18px;
        padding: 12px 15px;
        border-top: 1px solid var(--app-border);
        color: var(--app-muted);
        font-size: 13px;
        line-height: 1.5;
        overflow-wrap: anywhere;
      }
      .app-brand-footer strong {
        color: var(--app-accent);
      }
      .setup-shell {
        display: grid;
        grid-template-columns: minmax(280px, 340px) minmax(0, 1fr);
        gap: 18px;
        align-items: start;
        margin-top: 12px;
      }
	      .setup-sidebar {
	        position: sticky;
	        top: 12px;
          z-index: 3500;
	        background: #ffffff;
	        border: 1px solid var(--app-border);
	        border-top: 3px solid var(--app-accent);
	        border-radius: 8px;
	        padding: 14px;
		        box-shadow: 0 4px 12px rgba(39, 70, 144, 0.07);
	      }
      .setup-main {
        min-width: 0;
        position: relative;
        z-index: 1;
      }
      .setup-sidebar-section {
        border-bottom: 1px solid #e8eef5;
        padding-bottom: 12px;
        margin-bottom: 14px;
      }
      .setup-sidebar-section:last-child {
        border-bottom: 0;
        margin-bottom: 0;
        padding-bottom: 0;
      }
      .setup-sidebar .result-export-actions {
        margin-bottom: 24px;
        padding-bottom: 18px;
        border-bottom: 1px solid #e8eef5;
      }
	      .setup-sidebar-section h4,
	      .app-card h4 {
	        margin-top: 0;
	        margin-bottom: 10px;
	        font-size: 15px;
	        font-weight: 700;
	        color: var(--app-ink);
	      }
      .setup-sidebar .form-group {
        margin-bottom: 10px;
      }
	      .setup-sidebar .btn-primary {
	        width: 100%;
	      }
      .setup-sidebar .action-button,
      .setup-sidebar .shiny-download-link {
        width: 100%;
      }
      .setup-sidebar-actions .workflow-action:first-child {
        margin-top: 0;
      }
      .setup-sidebar .workflow-action > .shiny-html-output,
      .setup-sidebar .setup-prepare-control {
        display: block;
        width: 100%;
      }
		      .btn-primary,
		      .btn.btn-primary,
		      .action-button.btn-primary,
		      input.btn-primary {
		        background-color: var(--app-accent) !important;
		        border-color: var(--app-accent) !important;
		        color: #ffffff !important;
		        font-weight: 700;
		      }
		      .btn-primary:visited,
		      .btn-primary:active,
		      .btn.btn-primary:visited,
		      .btn.btn-primary:active,
		      .action-button.btn-primary:visited,
		      .action-button.btn-primary:active {
		        color: #ffffff !important;
		      }
		      .btn-primary:hover,
		      .btn-primary:focus,
		      .btn.btn-primary:hover,
		      .btn.btn-primary:focus,
		      .action-button.btn-primary:hover,
		      .action-button.btn-primary:focus,
		      input.btn-primary:hover,
		      input.btn-primary:focus {
		        background-color: #1f3978 !important;
		        border-color: #1f3978 !important;
		        color: #ffffff !important;
		      }
		      .btn-primary[disabled],
		      .btn-primary.disabled,
		      .btn.btn-primary[disabled],
		      .btn.btn-primary.disabled,
		      .action-button.btn-primary[disabled],
		      .action-button.btn-primary.disabled,
		      input.btn-primary[disabled],
		      input.btn-primary.disabled {
		        background-color: #6f7fb3 !important;
		        border-color: #6f7fb3 !important;
		        color: #ffffff !important;
		        opacity: 0.88;
		      }
	      .btn-default {
	        background: #ffffff;
		        border-color: #c6d6f5;
	        color: var(--app-accent);
	        font-weight: 700;
	      }
	      .btn-default:hover,
	      .btn-default:focus {
	        background: var(--app-accent-soft);
		        border-color: var(--app-blue);
	        color: var(--app-accent);
	      }
		      .app-card {
		        background: #ffffff;
		        border: 1px solid var(--app-border);
		        border-radius: 8px;
	        padding: 14px;
	        margin-bottom: 14px;
		        box-shadow: 0 2px 8px rgba(22, 34, 51, 0.04);
		        min-width: 0;
		      }
		      .tutorial-shell {
		        max-width: 1080px;
		        margin: 18px auto 0 auto;
		      }
			      .tutorial-hero {
				        background: #f6f8ff;
				        border: 1px solid #d9e3ff;
				        border-left: 5px solid var(--app-blue);
			        border-radius: 8px;
			        padding: 22px 24px;
			        margin-bottom: 14px;
				        box-shadow: 0 5px 18px rgba(39, 70, 144, 0.10);
			      }
			      .tutorial-eyebrow {
				        color: var(--app-accent);
			        font-size: 12px;
		        font-weight: 700;
		        text-transform: uppercase;
		        margin-bottom: 6px;
		      }
			      .tutorial-hero h2 {
			        margin: 0 0 8px 0;
			        color: var(--app-ink);
		        font-size: 28px;
		        font-weight: 700;
		        line-height: 1.15;
		      }
		      .tutorial-hero p {
		        color: var(--app-muted);
		        font-size: 15px;
		        line-height: 1.5;
		        margin: 0;
		        max-width: 780px;
		      }
	      .tutorial-workflow {
	        display: grid;
	        grid-template-columns: repeat(5, minmax(0, 1fr));
	        gap: 10px;
	        align-items: stretch;
	        margin-bottom: 14px;
	      }
	      .tutorial-step {
	        background: #ffffff;
	        border: 1px solid var(--app-border);
	        border-top: 3px solid var(--app-blue);
	        border-radius: 8px;
	        padding: 12px;
	        box-shadow: 0 2px 8px rgba(22, 34, 51, 0.04);
	        height: 100%;
	      }
			      .tutorial-step:nth-child(2),
			      .tutorial-step:nth-child(5) {
				        border-top-color: var(--app-violet);
			      }
			      .tutorial-step:nth-child(3) {
				        border-top-color: var(--app-accent);
			      }
			      .tutorial-step-number {
		        width: 26px;
		        height: 26px;
		        border-radius: 999px;
			        background: var(--app-accent-soft);
			        color: var(--app-accent);
		        display: flex;
		        align-items: center;
		        justify-content: center;
		        font-weight: 700;
		        font-size: 12px;
			        margin-bottom: 9px;
			      }
			      .tutorial-step:nth-child(2) .tutorial-step-number,
			      .tutorial-step:nth-child(5) .tutorial-step-number {
				        background: var(--app-violet-soft);
				        color: var(--app-violet);
			      }
			      .tutorial-step:nth-child(3) .tutorial-step-number {
				        background: var(--app-accent-soft);
				        color: var(--app-accent);
			      }
			      .tutorial-step h4,
			      .tutorial-tile h4 {
			        margin: 0 0 6px 0;
			        color: var(--app-ink);
		        font-size: 14px;
		        font-weight: 700;
		      }
		      .tutorial-step p,
		      .tutorial-tile p {
		        margin: 0;
		        color: var(--app-muted);
		        font-size: 12px;
		        line-height: 1.4;
		      }
	      .tutorial-tile-grid {
	        display: grid;
	        grid-template-columns: repeat(3, minmax(0, 1fr));
	        gap: 10px;
	        align-items: stretch;
	        margin-bottom: 14px;
	      }
	      .tutorial-tile {
	        background: var(--app-soft);
	        border: 1px solid #e5edf5;
	        border-top: 3px solid var(--app-blue);
	        border-radius: 8px;
	        padding: 12px;
	        height: 100%;
	      }
			      .tutorial-tile:nth-child(2) {
				        border-top-color: var(--app-accent);
				        background: var(--app-accent-soft);
			      }
			      .tutorial-tile:nth-child(3) {
				        border-top-color: var(--app-violet);
				        background: var(--app-violet-soft);
			      }
			      .tutorial-section {
			        background: #ffffff;
	        border: 1px solid var(--app-border);
	        border-radius: 8px;
	        margin-bottom: 10px;
	        box-shadow: 0 1px 2px rgba(22, 34, 51, 0.04);
	        overflow: hidden;
	      }
		      .tutorial-section summary {
		        cursor: pointer;
		        padding: 12px 14px;
		        color: var(--app-accent);
		        font-weight: 700;
		      }
	      .tutorial-body {
	        border-top: 1px solid #eef3f8;
	        padding: 12px 14px 14px 14px;
	        color: var(--app-muted);
	        font-size: 13px;
	        line-height: 1.5;
	      }
	      .tutorial-body p {
	        margin: 0 0 10px 0;
	      }
		      .tutorial-body p:last-child,
		      .tutorial-body ol,
		      .tutorial-body ul {
		        margin-bottom: 0;
		      }
	      .setup-preview-table .dataTables_wrapper {
	        width: 100%;
	      }
      .setup-preview-table .dataTables_scroll {
        overflow: visible;
      }
      .setup-preview-table .dataTables_scrollBody {
        border-bottom: 1px solid #d9e2ec;
      }
	      .setup-main-grid {
	        display: grid;
	        grid-template-columns: repeat(2, minmax(0, 1fr));
	        gap: 14px;
	        align-items: stretch;
	      }
	      .setup-main-grid > .app-card,
	      .overview-summary-grid > .app-card,
	      .pca-detail-grid > .app-card {
	        height: 100%;
	        margin-bottom: 0;
	      }
      .setup-empty {
        min-height: 220px;
        display: flex;
        flex-direction: column;
        justify-content: center;
        background: var(--app-soft);
      }
      .setup-empty p,
      .setup-hint,
      .setup-meta {
        color: var(--app-muted);
      }
      .setup-validation-list {
        list-style: none;
        padding-left: 0;
        margin-bottom: 0;
      }
      .setup-validation-list li {
        padding: 5px 0;
        border-bottom: 1px solid #eef3f8;
      }
      .setup-validation-list li:last-child {
        border-bottom: 0;
      }
      .setup-ok {
        color: #1b6e3c;
        font-weight: 700;
        margin-right: 6px;
      }
      .setup-warning {
        color: #8a5a00;
        font-weight: 700;
        margin-right: 6px;
      }
      .setup-status {
        margin-top: 8px;
        font-size: 13px;
      }
      .setup-stage-card {
        max-width: 980px;
      }
      .setup-stage-header {
        display: flex;
        justify-content: space-between;
        align-items: center;
        gap: 14px;
        margin-bottom: 12px;
      }
      .setup-stage-title {
        min-width: 220px;
      }
      .setup-stage-header h3 {
        margin: 0 0 4px 0;
        color: var(--app-ink);
        font-size: 20px;
        font-weight: 800;
      }
      .setup-stage-header h4 {
        margin: 0;
        color: var(--app-ink);
        font-size: 16px;
        font-weight: 800;
      }
      .setup-stage-header p {
        margin: 0;
        color: var(--app-muted);
        font-size: 13px;
      }
      .setup-stage-actions {
        display: flex;
        flex-wrap: wrap;
        gap: 8px;
        justify-content: flex-end;
        align-items: center;
        flex: 0 0 auto;
      }
      .setup-stage-header > .setup-stage-actions:only-child {
        margin-left: auto;
      }
      .setup-stage-header > h4 + .setup-stage-actions {
        margin-left: auto;
      }
      .setup-stage-actions > .shiny-html-output,
      .setup-prepare-control {
        display: inline-flex;
        align-items: center;
      }
      .setup-stage-actions .btn {
        margin-bottom: 0;
      }
      .setup-compact-file-grid {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 10px;
      }
      .setup-compact-file {
        border: 1px solid #e2eaf4;
        border-radius: 8px;
        background: #ffffff;
        padding: 11px 12px;
      }
      .setup-compact-file strong {
        display: block;
        color: var(--app-ink);
        font-size: 13px;
        margin-bottom: 4px;
      }
      .setup-compact-file span {
        color: var(--app-muted);
        font-size: 12px;
      }
      .setup-link-actions {
        display: flex;
        flex-wrap: wrap;
        gap: 8px;
      }
      .setup-role-grid {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 10px 14px;
      }
      .setup-role-grid .form-group {
        margin-bottom: 0;
      }
      .setup-preview-panel {
        min-width: 0;
      }
      .setup-preview-panel h5 {
        margin: 0 0 8px 0;
        color: var(--app-ink);
        font-size: 14px;
        font-weight: 800;
      }
      .setup-preview-panel .setup-metric-grid {
        margin-bottom: 10px;
      }
      .setup-ready-row {
        color: var(--app-muted);
        font-size: 13px;
        margin-top: 7px;
      }
      .setup-collapsed-grid {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 12px;
      }
      .setup-collapsed-grid details,
      .setup-wide-details {
        background: #ffffff;
        border: 1px solid var(--app-border);
        border-radius: 8px;
        padding: 0;
        overflow: hidden;
      }
      .setup-collapsed-grid summary,
      .setup-wide-details summary {
        cursor: pointer;
        color: var(--app-accent);
        font-weight: 800;
        padding: 11px 12px;
      }
      .setup-detail-body {
        border-top: 1px solid #eef3f8;
        padding: 12px;
      }
      .setup-intro-card {
        background: linear-gradient(135deg, #f8fbff 0%, #ffffff 55%, #f6f2ff 100%);
        border-left: 4px solid var(--app-blue);
      }
      .setup-intro-card p {
        color: var(--app-muted);
        margin-bottom: 12px;
        max-width: 820px;
      }
      .setup-file-grid {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 12px;
        margin-top: 12px;
      }
      .setup-file-example {
        display: grid;
        grid-template-columns: 42px minmax(0, 1fr);
        gap: 11px;
        align-items: start;
        background: #ffffff;
        border: 1px solid #e2eaf4;
        border-radius: 8px;
        padding: 12px;
      }
      .setup-file-icon {
        width: 42px;
        height: 42px;
        border-radius: 8px;
        background: var(--app-blue-soft);
        color: var(--app-accent);
        display: flex;
        align-items: center;
        justify-content: center;
        font-weight: 900;
      }
      .setup-file-example:nth-child(2) .setup-file-icon {
        background: var(--app-violet-soft);
        color: var(--app-violet);
      }
      .setup-file-example h5 {
        margin: 0 0 4px 0;
        font-size: 14px;
        font-weight: 800;
        color: var(--app-ink);
      }
      .setup-file-example p {
        margin: 0 0 8px 0;
        font-size: 12px;
      }
      .setup-column-chips {
        display: flex;
        flex-wrap: wrap;
        gap: 6px;
      }
      .setup-column-chips span {
        background: #f3f6fb;
        border: 1px solid #dfe8f3;
        border-radius: 999px;
        color: #526477;
        font-size: 11px;
        font-weight: 700;
        padding: 4px 8px;
      }
      .setup-why-card {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 10px;
        margin-top: 12px;
      }
      .setup-why-item {
        background: #ffffff;
        border: 1px solid #e2eaf4;
        border-radius: 8px;
        padding: 11px 12px;
      }
      .setup-why-item strong {
        display: block;
        color: var(--app-ink);
        font-size: 13px;
        margin-bottom: 4px;
      }
      .setup-why-item span {
        color: var(--app-muted);
        font-size: 12px;
        line-height: 1.4;
      }
      .setup-inline-help-row {
        display: flex;
        flex-wrap: wrap;
        align-items: center;
        gap: 12px;
        color: var(--app-muted);
        font-size: 12px;
        margin: 8px 0 0 0;
      }
      .setup-section-kicker {
        margin: 0 0 8px 0;
        color: var(--app-ink);
        font-size: 13px;
        font-weight: 900;
      }
      .setup-question-strip {
        border: 1px solid #e1eaf5;
        border-radius: 8px;
        background: #fbfdff;
        padding: 10px 12px;
        margin: 0 0 12px 0;
      }
      .setup-question-flow {
        display: grid;
        grid-template-columns: minmax(0, 1fr) 34px minmax(0, 1.15fr) 34px minmax(0, 1fr);
        gap: 8px;
        align-items: center;
      }
      .setup-question-node {
        min-width: 0;
        border: 1px solid #dfe8f3;
        border-radius: 8px;
        background: #ffffff;
        padding: 9px 10px;
        text-align: center;
      }
      .setup-question-node strong {
        display: block;
        color: var(--app-ink);
        font-size: 13px;
        margin-bottom: 3px;
      }
      .setup-question-node span {
        display: block;
        color: var(--app-muted);
        font-size: 12px;
        line-height: 1.25;
      }
      .setup-question-node .concept-help-target {
        display: inline-flex;
        width: auto;
      }
      .setup-question-arrow {
        color: var(--app-accent);
        font-size: 18px;
        font-weight: 900;
        text-align: center;
      }
      .setup-question-note {
        margin: 8px 0 0 0;
        color: var(--app-muted);
        font-size: 12px;
        line-height: 1.35;
      }
      .spatial-reference-popover {
        display: inline-flex;
        margin-left: auto;
        position: relative;
      }
      .spatial-reference-trigger {
        align-items: center;
        background: #ffffff;
        border: 1px solid #d7e4f3;
        border-radius: 999px;
        color: var(--app-accent);
        cursor: help;
        display: inline-flex;
        font-size: 12px;
        font-weight: 800;
        gap: 6px;
        line-height: 1.2;
        padding: 6px 10px;
      }
      .spatial-reference-trigger:hover,
      .spatial-reference-trigger:focus {
        background: var(--app-blue-soft);
        border-color: #bcd2ee;
        outline: none;
      }
      .spatial-reference-panel {
        background: #fbfdff;
        border: 1px solid #dfe8f3;
        border-radius: 8px;
        box-shadow: 0 14px 34px rgba(36, 52, 71, 0.16);
        color: var(--app-muted);
        display: none;
        max-height: min(70vh, 520px);
        overflow: auto;
        padding: 10px 12px 12px 12px;
        position: absolute;
        right: 0;
        top: calc(100% + 8px);
        width: min(760px, calc(100vw - 48px));
        z-index: 3000;
      }
      .spatial-reference-popover.is-open .spatial-reference-panel {
        display: block;
      }
      .spatial-reference-panel::before {
        background: #fbfdff;
        border-left: 1px solid #dfe8f3;
        border-top: 1px solid #dfe8f3;
        content: '';
        height: 10px;
        position: absolute;
        right: 22px;
        top: -6px;
        transform: rotate(45deg);
        width: 10px;
        z-index: 3001;
      }
      .spatial-reference-header {
        align-items: center;
        display: flex;
        gap: 10px;
        justify-content: space-between;
        margin-bottom: 4px;
      }
      .spatial-reference-header h5 {
        margin: 0;
        color: var(--app-ink);
        font-size: 13px;
        font-weight: 900;
      }
      .spatial-reference-close {
        align-items: center;
        background: transparent;
        border: 0;
        border-radius: 999px;
        color: #7c8ca1;
        display: inline-flex;
        font-size: 18px;
        font-weight: 800;
        height: 24px;
        justify-content: center;
        line-height: 1;
        padding: 0;
        width: 24px;
      }
      .spatial-reference-close:hover,
      .spatial-reference-close:focus {
        background: #eef4ff;
        color: var(--app-accent);
        outline: none;
      }
      .spatial-reference-panel p {
        margin: 0 0 8px 0;
        color: var(--app-muted);
        font-size: 12px;
        line-height: 1.4;
      }
      .setup-sidebar .spatial-reference-popover {
        display: flex;
        margin: 10px 0 0 0;
        width: 100%;
      }
      .setup-sidebar .spatial-reference-trigger {
        justify-content: center;
        width: 100%;
      }
      .setup-sidebar .spatial-reference-panel {
        left: 0;
        right: auto;
      }
      .setup-sidebar .spatial-reference-panel::before {
        left: 22px;
        right: auto;
      }
      .spatial-reference-table {
        width: 100%;
        border-collapse: collapse;
        background: #ffffff;
        border: 1px solid #e2eaf4;
        border-radius: 8px;
        overflow: hidden;
        font-size: 12px;
      }
      .spatial-reference-table th {
        background: #f3f6fb;
        color: var(--app-ink);
        font-size: 11px;
        font-weight: 900;
        padding: 7px 8px;
        text-align: left;
        text-transform: uppercase;
      }
      .spatial-reference-table td {
        border-top: 1px solid #eef3f8;
        color: var(--app-muted);
        line-height: 1.35;
        padding: 7px 8px;
        vertical-align: top;
      }
      .spatial-reference-table td:first-child {
        color: #5c6d82;
        font-weight: 700;
      }
      .spatial-reference-boundary {
        margin-top: 8px;
        color: var(--app-muted);
        font-size: 12px;
        line-height: 1.4;
      }
      .spatial-reference-boundary span {
        display: block;
        margin-top: 3px;
      }
      .setup-example-grid {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 12px;
        margin-top: 12px;
        align-items: stretch;
      }
      .setup-spreadsheet-example {
        display: flex;
        flex-direction: column;
        min-width: 0;
      }
      .setup-spreadsheet-example h5 {
        margin: 0 0 7px 0;
        color: var(--app-ink);
        font-size: 13px;
        font-weight: 800;
      }
      .setup-spreadsheet-example > h5 + .setup-card-note {
        min-height: 52px;
      }
      .setup-match-cue {
        display: grid;
        grid-template-columns: minmax(0, 1fr) 110px minmax(0, 1fr);
        gap: 10px;
        align-items: center;
        margin-top: 12px;
        padding: 11px 12px;
        border: 1px dashed #cbd9ea;
        border-radius: 8px;
        background: #fbfdff;
      }
      .setup-match-cue strong,
      .setup-match-cue span {
        display: block;
      }
      .setup-match-cue strong {
        color: var(--app-ink);
        font-size: 12px;
        margin-bottom: 6px;
      }
      .setup-match-cue span {
        color: var(--app-muted);
        font-size: 12px;
        font-variant-numeric: tabular-nums;
      }
      .setup-match-lines {
        text-align: center;
        color: var(--app-blue);
        font-size: 11px;
        font-weight: 800;
        text-transform: uppercase;
      }
      .setup-match-lines span::before,
      .setup-match-lines span::after {
        content: '';
        display: inline-block;
        width: 18px;
        height: 1px;
        background: #9db6dc;
        vertical-align: middle;
        margin: 0 5px;
      }
      .setup-template-actions {
        display: flex;
        flex-wrap: wrap;
        gap: 8px;
        margin-top: 12px;
      }
      .setup-template-actions .btn {
        font-weight: 700;
      }
      .setup-example-footer {
        margin-top: auto;
        padding-top: 10px;
      }
      .setup-example-footer .btn {
        min-height: 34px;
        display: inline-flex;
        align-items: center;
      }
      .help-tooltip.setup-icon-help-target {
        display: inline-flex;
        align-items: center;
        gap: 5px;
        width: auto;
        height: auto;
        margin-left: 0;
        padding: 0;
        border: 0;
        border-radius: 0;
        background: transparent;
        color: var(--app-muted);
        cursor: help;
        font-size: 12px;
        font-weight: 700;
        line-height: 1.35;
      }
      .help-tooltip.setup-icon-only-tooltip {
        display: inline-flex;
        align-items: center;
        justify-content: center;
        color: #7c8ca1;
        cursor: help;
        line-height: 1;
      }
      .setup-spreadsheet-example h5 .setup-icon-only-tooltip {
        margin-left: 0;
        background: transparent;
        border: 0;
      }
      .setup-inline-help-row .setup-icon-only-tooltip {
        background: transparent;
        border: 0;
      }
      .setup-icon-help-target:hover,
      .setup-icon-help-target:focus,
      .setup-icon-only-tooltip:hover,
      .setup-icon-only-tooltip:focus {
        color: var(--app-accent);
        outline: none;
      }
      .setup-icon-help-target .glyphicon,
      .setup-icon-only-tooltip .glyphicon,
      .setup-heading-icon {
        color: #7c8ca1;
        font-size: 12px;
        flex: 0 0 auto;
      }
      .setup-icon-only-tooltip:hover .glyphicon,
      .setup-icon-only-tooltip:focus .glyphicon {
        color: var(--app-accent);
      }
      .setup-heading-icon {
        margin-right: 6px;
      }
      .setup-spreadsheet-example h5 {
        display: flex;
        align-items: center;
        gap: 5px;
        margin: 0 0 9px 0;
        color: var(--app-ink);
        font-size: 14px;
        font-weight: 800;
      }
      .setup-spreadsheet-example > .setup-card-note {
        margin: -2px 0 8px 0;
        min-height: 52px;
      }
      .setup-spreadsheet-scroll {
        width: 100%;
        max-width: 100%;
        overflow-x: auto;
        border: 1px solid #dfe7f0;
        border-radius: 8px;
        background: #ffffff;
      }
      .setup-spreadsheet-table {
        width: max-content;
        min-width: 100%;
        border-collapse: collapse;
        color: var(--app-ink);
        font-size: 12px;
      }
      .setup-spreadsheet-table th,
      .setup-spreadsheet-table td {
        border-right: 1px solid #e6edf5;
        border-bottom: 1px solid #e6edf5;
        padding: 7px 9px;
        text-align: left;
        white-space: nowrap;
      }
      .setup-spreadsheet-table tr:last-child th,
      .setup-spreadsheet-table tr:last-child td {
        border-bottom: 0;
      }
      .setup-spreadsheet-table th:last-child,
      .setup-spreadsheet-table td:last-child {
        border-right: 0;
      }
      .sheet-corner,
      .sheet-col-letter,
      .sheet-row-number {
        background: #f3f6fa;
        color: #6c7b8e;
        font-weight: 800;
      }
      .sheet-corner,
      .sheet-row-number {
        width: 34px;
        min-width: 34px;
        text-align: center !important;
      }
      .sheet-col-letter {
        text-align: center !important;
      }
      .sheet-variable-name {
        background: #f9fbfd;
        color: var(--app-ink);
        font-weight: 800;
      }
      .setup-requirement-stack {
        display: grid;
        gap: 5px;
        margin-top: 9px;
      }
      .setup-requirement-line {
        display: flex;
        align-items: center;
        gap: 5px;
        margin: 0;
        color: var(--app-muted);
        font-size: 12px;
        line-height: 1.35;
      }
      .setup-requirement-line strong {
        color: var(--app-ink);
      }
      .setup-workflow-list {
        display: grid;
        gap: 9px;
      }
      .setup-workflow-step {
        border: 1px solid #e5edf5;
        border-radius: 8px;
        background: #ffffff;
        padding: 10px 11px;
      }
      .setup-workflow-step strong {
        display: block;
        color: var(--app-ink);
        font-size: 13px;
        margin-bottom: 3px;
      }
      .setup-workflow-step span {
        color: var(--app-muted);
        font-size: 12px;
        line-height: 1.35;
      }
      .setup-static-preview {
        width: 100%;
        overflow-x: auto;
        border: 1px solid #e2eaf4;
        border-radius: 8px;
      }
      .setup-static-preview table {
        width: max-content;
        min-width: 100%;
        border-collapse: collapse;
        font-size: 12px;
      }
      .setup-static-preview th,
      .setup-static-preview td {
        border-bottom: 1px solid #eef3f8;
        padding: 7px 9px;
        text-align: left;
        white-space: nowrap;
      }
      .setup-static-preview th {
        background: #f6f8fb;
        color: var(--app-muted);
        font-weight: 800;
      }
      .setup-static-preview tr:last-child td {
        border-bottom: 0;
      }
      .setup-card-note {
        color: var(--app-muted);
        font-size: 12px;
        margin: 8px 0 0 0;
      }
      .setup-requirement-stack {
        display: grid;
        gap: 5px;
        margin-top: 9px;
      }
      .setup-requirement-stack .setup-card-note {
        margin: 0;
      }
      .setup-compact-guidance {
        display: flex;
        flex-wrap: wrap;
        align-items: center;
        justify-content: space-between;
        gap: 10px 14px;
        margin-top: 12px;
        padding: 10px 12px;
        background: #ffffff;
        border: 1px solid #e2eaf4;
        border-radius: 8px;
      }
      .setup-compact-guidance h5 {
        margin: 0 0 3px 0;
        color: var(--app-ink);
        font-size: 13px;
        font-weight: 800;
      }
      .setup-compact-guidance p {
        margin: 0;
        color: var(--app-muted);
        font-size: 12px;
      }
      .setup-compact-help {
        display: flex;
        flex-wrap: wrap;
        gap: 10px;
        align-items: center;
        font-size: 12px;
        color: var(--app-muted);
        font-weight: 700;
      }
      .setup-metric-grid {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 8px;
      }
      .setup-metric-row {
        background: var(--app-soft);
        border: 1px solid #e5edf5;
        border-radius: 8px;
        padding: 9px 10px;
      }
      .setup-metric-row span {
        display: block;
        color: var(--app-muted);
        font-size: 11px;
        font-weight: 700;
        text-transform: uppercase;
      }
      .setup-metric-row strong {
        display: block;
        color: var(--app-ink);
        font-size: 16px;
        margin-top: 3px;
      }
      .setup-success-banner,
      .setup-error-banner,
      .setup-warning-banner {
        border-radius: 8px;
        padding: 11px 12px;
        margin-bottom: 10px;
      }
      .setup-success-banner {
        background: #f3fbf8;
        border: 1px solid #cdece0;
        color: #1b6e3c;
      }
      .setup-error-banner {
        background: #fff5f5;
        border: 1px solid #f3cccc;
        color: #9b1c1c;
      }
      .setup-warning-banner {
        background: #fff9ed;
        border: 1px solid #f2ddb7;
        color: #8a5a00;
      }
      .setup-detail-toggle {
        color: var(--app-muted);
        font-size: 12px;
      }
      .setup-detail-toggle summary {
        cursor: pointer;
        font-weight: 800;
        color: var(--app-accent);
        margin-bottom: 8px;
      }
      .help-tooltip {
        display: inline-flex;
        align-items: center;
        justify-content: center;
        position: relative;
        width: 18px;
        height: 18px;
        margin-left: 6px;
        border-radius: 999px;
        background: var(--app-blue-soft);
        border: 1px solid #c6d6f5;
        color: var(--app-accent);
        cursor: help;
        font-size: 12px;
        font-weight: 800;
        line-height: 1;
      }
      .help-tooltip:hover::after,
      .help-tooltip:focus::after {
        content: attr(data-tooltip);
        position: absolute;
        left: 50%;
        top: calc(100% + 9px);
        transform: translateX(-50%);
        width: max-content;
        max-width: 280px;
        background: #243447;
        color: #ffffff;
        border-radius: 8px;
        box-shadow: 0 8px 22px rgba(22, 34, 51, 0.18);
        font-size: 12px;
        font-weight: 600;
        line-height: 1.4;
        padding: 8px 10px;
        text-align: left;
        white-space: normal;
        z-index: 2000;
      }
      .help-tooltip:hover::before,
      .help-tooltip:focus::before {
        content: '';
        position: absolute;
        left: 50%;
        top: calc(100% + 3px);
        transform: translateX(-50%);
        border: 6px solid transparent;
        border-bottom-color: #243447;
        z-index: 2001;
      }
      .concept-help-target {
        position: relative;
        display: inline-flex;
        align-items: baseline;
        gap: 6px;
        width: max-content;
        max-width: 100%;
        color: inherit;
        border-bottom: 1px dotted rgba(105, 122, 145, 0.65);
        cursor: help;
        outline: none;
        text-decoration: none;
      }
      .concept-help-target:hover,
      .concept-help-target:focus {
        border-bottom-color: var(--app-accent);
      }
      .concept-help-target:hover::after,
      .concept-help-target:focus::after {
        content: attr(data-tooltip);
        position: absolute;
        z-index: 2200;
        left: 50%;
        bottom: calc(100% + 9px);
        transform: translateX(-50%);
        width: max-content;
        max-width: min(320px, 80vw);
        white-space: normal;
        border: 1px solid #ccd8e6;
        border-radius: 8px;
        background: #ffffff;
        color: var(--app-ink);
        box-shadow: 0 12px 28px rgba(31, 42, 68, 0.16);
        padding: 8px 10px;
        font-size: 12px;
        font-weight: 700;
        line-height: 1.35;
        text-align: left;
        pointer-events: none;
      }
      .concept-help-target:hover::before,
      .concept-help-target:focus::before {
        content: '';
        position: absolute;
        z-index: 2201;
        left: 50%;
        bottom: calc(100% + 3px);
        transform: translateX(-50%);
        border-left: 6px solid transparent;
        border-right: 6px solid transparent;
        border-top: 6px solid #ffffff;
        pointer-events: none;
      }
      .setup-sidebar .concept-help-target:hover::after,
      .setup-sidebar .concept-help-target:focus::after {
        left: 0;
        transform: none;
      }
      .setup-sidebar .concept-help-target:hover::before,
      .setup-sidebar .concept-help-target:focus::before {
        left: 9px;
        transform: none;
      }
      .setup-sidebar,
      .setup-sidebar-section,
      .app-card {
        overflow: visible;
      }
      .setup-sidebar .help-tooltip:hover::after,
      .setup-sidebar .help-tooltip:focus::after {
        left: 0;
        right: auto;
        transform: none;
      }
      .setup-sidebar .help-tooltip:hover::before,
      .setup-sidebar .help-tooltip:focus::before {
        left: 9px;
        transform: none;
      }
      .setup-requirement-list {
        margin: 8px 0 0 0;
        padding-left: 20px;
      }
      .setup-requirement-list > li {
        margin-bottom: 10px;
      }
      .setup-requirement-list ul {
        color: var(--app-muted);
        font-size: 13px;
        margin-top: 5px;
        padding-left: 18px;
      }
	      .workflow-context-card {
	        background: var(--app-soft);
	        border: 1px solid #e5edf5;
	        border-left: 3px solid var(--app-blue);
	        border-radius: 8px;
        padding: 11px 12px;
        margin-bottom: 10px;
      }
      .workflow-context-card h5 {
        margin: 0 0 8px 0;
        font-size: 13px;
        font-weight: 700;
        color: #243447;
      }
      .workflow-context-card .app-kv-item {
        padding: 5px 0;
        font-size: 12px;
      }
      .workflow-action {
        margin-top: 10px;
      }
      .workflow-action .btn,
      .workflow-action .btn-primary {
        width: 100%;
      }
      .workflow-note {
        color: var(--app-muted);
        font-size: 12px;
        margin: 6px 0 0 0;
      }
	      .workflow-note.warning {
	        color: #8a5a00;
	      }
	      .workflow-note.ok {
	        color: var(--app-accent);
	      }
      .overview-shell {
        display: grid;
        grid-template-columns: minmax(280px, 340px) minmax(0, 1fr);
        gap: 18px;
        align-items: start;
        margin-top: 12px;
      }
      .overview-main {
        min-width: 0;
      }
      .overview-spatial-card {
        padding: 14px 16px 10px 16px;
      }
      .overview-plot-grid {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 12px;
        align-items: start;
      }
      .overview-plot-panel {
        min-width: 0;
      }
      .overview-plot-panel .shiny-plot-output {
        width: 100% !important;
      }
      .overview-metric-row {
        display: grid;
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 10px;
        margin-bottom: 12px;
      }
	      .overview-metric {
	        background: var(--app-blue-soft);
	        border: 1px solid #d6e6ff;
	        border-radius: 8px;
	        padding: 11px 12px;
	      }
	      .overview-metric-value {
	        display: block;
	        color: var(--app-accent);
        font-size: 24px;
        font-weight: 700;
        line-height: 1.1;
      }
      .overview-metric-label {
        display: block;
        color: var(--app-muted);
        font-size: 12px;
        margin-top: 4px;
      }
      .app-kv-item {
        display: flex;
        justify-content: space-between;
        gap: 12px;
        padding: 8px 0;
        border-bottom: 1px solid #eef3f8;
      }
      .app-kv-item:last-child {
        border-bottom: 0;
      }
      .app-kv-label {
        color: var(--app-muted);
      }
	      .app-kv-value {
	        color: var(--app-ink);
	        font-weight: 600;
	        text-align: right;
	      }
      .pca-shell {
        display: grid;
        grid-template-columns: minmax(280px, 340px) minmax(0, 1fr);
        gap: 18px;
        align-items: start;
        margin-top: 12px;
      }
      .pca-main {
        min-width: 0;
      }
      .pca-main-actions {
        display: flex;
        justify-content: flex-end;
        margin-bottom: 12px;
      }
      .pca-section {
        margin-bottom: 14px;
      }
      .pca-section > h4 {
        margin-top: 0;
        margin-bottom: 10px;
        font-size: 15px;
        font-weight: 700;
        color: #243447;
      }
      .pca-view-tabs .nav-tabs {
        border-bottom: 1px solid var(--app-border);
        margin-bottom: 14px;
      }
	      .pca-view-tabs .nav-tabs > li > a {
	        color: var(--app-muted);
	        font-weight: 700;
	        border-radius: 8px 8px 0 0;
	        border: 1px solid transparent;
	        padding: 9px 14px;
	      }
	      .pca-view-tabs .nav-tabs > li > a:hover,
	      .pca-view-tabs .nav-tabs > li > a:focus {
	        color: var(--app-accent);
	        background: var(--app-blue-soft);
	      }
	      .pca-view-tabs .nav-tabs > li.active > a,
	      .pca-view-tabs .nav-tabs > li.active > a:focus,
	      .pca-view-tabs .nav-tabs > li.active > a:hover {
	        color: var(--app-accent);
	        background: var(--app-accent-soft);
	        border-color: #d9e3ff #d9e3ff var(--app-accent-soft) #d9e3ff;
	      }
      .pca-view-tabs .tab-content {
        min-width: 0;
      }
      .pca-loading-controls {
        display: grid;
        grid-template-columns: repeat(2, minmax(170px, 1fr));
        gap: 12px;
        align-items: start;
        margin-bottom: 12px;
      }
      .pca-loading-controls .form-group {
        margin-bottom: 0;
      }
	      .pca-detail-grid {
	        display: grid;
	        grid-template-columns: minmax(360px, 0.95fr) minmax(0, 1.05fr);
	        gap: 14px;
	        align-items: stretch;
	      }
      .pca-detail-map-card .shiny-plot-output {
        width: 100% !important;
      }
      .pca-pc-row {
        padding: 9px 0;
        border-bottom: 1px solid #e8eef5;
      }
      .pca-pc-row:last-child {
        border-bottom: 0;
      }
      .pca-pc-row .checkbox {
        margin-top: 0;
        margin-bottom: 4px;
        font-weight: 700;
      }
      .pca-pc-row .form-group {
        margin-bottom: 0;
      }
	      .pca-compact-summary {
	        padding: 5px 0;
	        border-bottom: 1px solid #eef3f8;
	      }
	      .pca-compact-summary:last-child {
	        border-bottom: 0;
	      }
	      .pca-summary-label {
	        color: var(--app-muted);
	        font-size: 12px;
	        margin-bottom: 5px;
	      }
	      .pca-selected-pills {
	        display: flex;
	        flex-wrap: nowrap;
	        gap: 4px;
	        align-items: center;
	        overflow: hidden;
	      }
		      .pca-pill {
	        display: inline-flex;
	        align-items: center;
	        justify-content: center;
		        border: 1px solid #c6d6f5;
		        background: var(--app-accent-soft);
		        border-radius: 8px;
		        color: var(--app-accent);
	        font-size: 12px;
	        font-weight: 700;
	        line-height: 1;
	        padding: 5px 7px;
	        white-space: nowrap;
	      }
	      .pca-summary-value {
	        color: #243447;
	        font-size: 12px;
	        font-weight: 600;
	        overflow-wrap: anywhere;
	      }
	      .pca-map-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
        gap: 14px;
      }
      .pca-map-grid.pc-count-1 {
        max-width: 560px;
        margin: 0 auto;
      }
      .pca-map-card {
        padding: 12px;
      }
      .pca-map-card .shiny-plot-output {
        width: 100% !important;
      }
      .pca-card-meta {
        color: var(--app-muted);
        font-size: 12px;
        margin-top: -5px;
        margin-bottom: 8px;
      }
      .pca-secondary-grid {
        display: grid;
        grid-template-columns: minmax(280px, 0.85fr) minmax(0, 1.15fr);
        gap: 14px;
      }
	      .loading-result-grid {
	        display: grid;
	        grid-template-columns: repeat(2, minmax(0, 1fr));
	        gap: 14px;
	        align-items: stretch;
	      }
      .loading-result-grid.single {
        grid-template-columns: minmax(280px, 680px);
        justify-content: center;
      }
	      .loading-panel {
	        background: #ffffff;
	        border: 1px solid #e5edf5;
	        border-top: 3px solid var(--app-blue);
	        border-radius: 8px;
	        padding: 12px;
	        min-width: 0;
	        height: 100%;
	      }
	      .loading-panel.negative {
	        border-top-color: var(--app-red);
	      }
	      .loading-panel h5 {
        margin-top: 0;
        margin-bottom: 10px;
	        color: var(--app-ink);
        font-weight: 700;
      }
      .loading-bar-row {
        padding: 7px 0;
        border-bottom: 1px solid #eef3f8;
      }
      .loading-bar-row:last-child {
        border-bottom: 0;
      }
      .loading-bar-head {
        display: flex;
        justify-content: space-between;
        gap: 12px;
        margin-bottom: 4px;
        font-size: 13px;
      }
	      .loading-gene {
	        color: var(--app-ink);
        font-weight: 600;
        overflow-wrap: anywhere;
      }
      .loading-value {
        color: var(--app-muted);
        font-variant-numeric: tabular-nums;
      }
	      .loading-bar-track {
	        height: 7px;
		        background: #edf2fb;
        border-radius: 999px;
        overflow: hidden;
      }
	      .loading-bar-fill {
	        height: 100%;
	        border-radius: 999px;
		        background: var(--app-blue);
	      }
	      .loading-bar-fill.negative {
	        background: var(--app-red);
	      }
      .mediation-shell {
        display: grid;
        grid-template-columns: minmax(280px, 340px) minmax(0, 1fr);
        gap: 18px;
        align-items: start;
        margin-top: 12px;
      }
      .mediation-main {
        min-width: 0;
      }
      .mediation-header {
        display: flex;
        justify-content: space-between;
        align-items: flex-start;
        gap: 14px;
        margin-bottom: 12px;
      }
	      .mediation-header h3 {
        margin: 0 0 4px 0;
        font-size: 20px;
        font-weight: 700;
	        color: var(--app-ink);
      }
      .mediation-header p {
        margin: 0;
      }
      .mediation-context {
        color: var(--app-muted);
        font-size: 13px;
        text-align: right;
      }
      .mediation-contrast-panel {
        margin-bottom: 14px;
      }
	      .mediation-contrast-panel h4 {
        margin-top: 0;
        margin-bottom: 10px;
        font-size: 15px;
        font-weight: 700;
	        color: var(--app-ink);
	      }
	      .effect-card-grid {
	        display: grid;
	        grid-template-columns: repeat(4, minmax(0, 1fr));
	        gap: 12px;
	        align-items: stretch;
	      }
	      .effect-card {
	        background: #ffffff;
	        border: 1px solid #e5edf5;
	        border-top: 3px solid var(--app-blue);
	        border-radius: 8px;
	        padding: 12px;
	        min-width: 0;
	        height: 100%;
	      }
	      .effect-card:nth-child(2) {
		        border-top-color: var(--app-violet);
	      }
	      .effect-card:nth-child(3) {
		        border-top-color: var(--app-accent);
	      }
	      .effect-card:nth-child(4) {
		        border-top-color: var(--app-slate-accent);
	      }
      .effect-label {
        color: var(--app-muted);
        font-size: 12px;
        margin-bottom: 6px;
      }
	      .effect-value {
	        color: var(--app-ink);
        font-size: 23px;
        font-weight: 700;
        line-height: 1.1;
        font-variant-numeric: tabular-nums;
      }
      .effect-interval,
      .effect-evidence {
        color: var(--app-muted);
        font-size: 12px;
        margin-top: 5px;
        font-variant-numeric: tabular-nums;
      }
      .mediation-result-grid {
        display: grid;
        grid-template-columns: minmax(0, 1fr);
        gap: 14px;
      }
      .mediation-dag-card {
        margin-bottom: 14px;
        overflow: visible;
      }
      .mediation-dag-layout {
        display: grid;
        grid-template-columns: minmax(0, 1fr);
        gap: 12px;
        align-items: stretch;
        overflow: visible;
      }
      .mediation-dag-canvas {
        position: relative;
        min-height: 340px;
        border: 1px solid #e2eaf4;
        border-radius: 10px;
        background: linear-gradient(180deg, #ffffff 0%, #f8fbff 100%);
        overflow: visible;
      }
      .mediation-node {
        position: absolute;
        z-index: 2;
        width: 138px;
        min-height: 70px;
        border: 1px solid #d7e3f1;
        border-radius: 10px;
        background: #ffffff;
        box-shadow: 0 10px 26px rgba(31, 42, 68, 0.08);
        padding: 12px 10px;
        text-align: center;
      }
      .mediation-node-x {
        left: 5%;
        bottom: 18px;
      }
      .mediation-node-m {
        left: 50%;
        top: 14px;
        transform: translateX(-50%);
        width: min(260px, 48%);
        min-height: 72px;
        border-top: 4px solid var(--app-accent);
      }
      .mediation-node-y {
        right: 5%;
        bottom: 18px;
      }
      .mediation-node-label {
        color: var(--app-muted);
        font-size: 12px;
        font-weight: 700;
        text-transform: uppercase;
      }
      .mediation-node-value {
        margin-top: 5px;
        color: var(--app-ink);
        font-size: 18px;
        font-weight: 800;
      }
      .mediation-node-detail {
        margin-top: 3px;
        color: var(--app-muted);
        font-size: 12px;
      }
      .mediation-arrow-svg {
        position: absolute;
        inset: 0;
        width: 100%;
        height: 100%;
        z-index: 1;
        overflow: visible;
        pointer-events: none;
      }
      .mediation-dag-line {
        fill: none;
        stroke-width: 4;
        stroke-linecap: round;
        vector-effect: non-scaling-stroke;
      }
      .mediation-dag-line-alpha {
        stroke: var(--app-accent);
      }
      .mediation-dag-line-beta {
        stroke: var(--app-violet);
      }
      .mediation-dag-line-direct {
        stroke: var(--app-slate-accent);
      }
      .mediation-arrow-label {
        position: absolute;
        z-index: 3;
        width: 176px;
        max-width: 220px;
        min-height: 70px;
        border: 1px solid #dfe7f0;
        border-radius: 10px;
        background: rgba(255, 255, 255, 0.94);
        color: var(--app-ink);
        font-size: 12px;
        font-weight: 800;
        line-height: 1.2;
        padding: 8px 10px;
        box-shadow: 0 8px 20px rgba(31, 42, 68, 0.08);
        display: flex;
        align-items: center;
        justify-content: center;
        text-align: center;
      }
      .mediation-arrow-label small {
        display: none;
      }
      .mediation-arrow-label-xm {
        left: 9%;
        top: 36%;
      }
      .mediation-arrow-label-my {
        right: 9%;
        top: 36%;
      }
      .mediation-arrow-label-xy {
        left: 50%;
        bottom: 13px;
        transform: translateX(-50%);
      }
      .mediation-arrow-label-xy .interval-tooltip {
        font-size: 18px;
        font-weight: 900;
      }
      .mediation-dag-meta {
        display: flex;
        flex-wrap: wrap;
        gap: 8px;
      }
      .mediation-dag-chip {
        border: 1px solid #dce6f2;
        border-radius: 999px;
        background: #ffffff;
        color: var(--app-muted);
        font-size: 12px;
        font-weight: 700;
        padding: 6px 10px;
      }
      .mediation-dag-formula {
        border: 1px solid #dce6f2;
        border-radius: 10px;
        background: #f7fbff;
        padding: 12px;
      }
      .mediation-dag-formula-head {
        display: flex;
        flex-wrap: wrap;
        align-items: center;
        gap: 8px;
      }
      .mediation-dag-summary-strip {
        display: grid;
        grid-template-columns: minmax(0, 1fr);
        gap: 12px;
        align-items: start;
      }
      .mediation-dag-formula strong {
        color: var(--app-ink);
      }
      .mediation-dag-formula code {
        color: var(--app-accent);
        background: transparent;
        padding: 0;
        font-weight: 800;
      }
      .mediation-path-symbol {
        color: var(--app-ink);
        font-size: 18px;
        font-weight: 900;
      }
      .mediation-path-symbol-wrap {
        display: inline-flex;
        align-items: center;
        justify-content: center;
        gap: 6px;
      }
      .mediation-path-lines {
        display: grid;
        gap: 3px;
        margin-top: 5px;
      }
      .mediation-path-line {
        display: flex;
        justify-content: space-between;
        gap: 10px;
        color: var(--app-muted);
        font-size: 13px;
        font-weight: 700;
        font-variant-numeric: tabular-nums;
      }
      .mediation-path-est {
        color: var(--app-ink);
        font-size: 14px;
      }
      .interval-tooltip {
        position: relative;
        cursor: help;
        border-bottom: 1px dotted #91a4bd;
      }
      .interval-tooltip:hover::after,
      .interval-tooltip:focus::after {
        content: attr(data-tooltip);
        position: absolute;
        z-index: 2000;
        left: 50%;
        bottom: calc(100% + 9px);
        transform: translateX(-50%);
        width: max-content;
        max-width: 260px;
        white-space: pre-line;
        border: 1px solid #ccd8e6;
        border-radius: 8px;
        background: #ffffff;
        color: var(--app-ink);
        box-shadow: 0 12px 28px rgba(31, 42, 68, 0.16);
        padding: 8px 10px;
        font-size: 12px;
        font-weight: 700;
        line-height: 1.35;
        text-align: left;
        pointer-events: none;
      }
      .interval-tooltip:hover::before,
      .interval-tooltip:focus::before {
        content: '';
        position: absolute;
        z-index: 2001;
        left: 50%;
        bottom: calc(100% + 3px);
        transform: translateX(-50%);
        border-left: 6px solid transparent;
        border-right: 6px solid transparent;
        border-top: 6px solid #ffffff;
        pointer-events: none;
      }
      .mediation-indirect-list {
        display: flex;
        flex-wrap: wrap;
        gap: 7px;
        margin-top: 10px;
      }
      .mediation-indirect-row {
        display: flex;
        justify-content: space-between;
        gap: 12px;
        min-width: 150px;
        border: 1px solid #e2eaf4;
        border-radius: 8px;
        background: #ffffff;
        color: var(--app-muted);
        font-size: 12px;
        font-weight: 700;
        font-variant-numeric: tabular-nums;
        padding: 7px 9px;
      }
      .mediation-indirect-row.total {
        border-color: #cfdcf2;
        background: #f5f8ff;
        color: var(--app-ink);
      }
      .mediation-indirect-value {
        color: var(--app-accent);
        font-weight: 900;
      }
      .mediation-effect-grid {
        display: grid;
        grid-template-columns: repeat(3, minmax(0, 1fr));
        gap: 10px;
      }
      .mediation-effect-card {
        border: 1px solid #dce6f2;
        border-radius: 10px;
        background: #ffffff;
        padding: 12px;
        min-width: 0;
      }
      .mediation-effect-card h5 {
        display: flex;
        align-items: baseline;
        gap: 8px;
        margin: 0 0 9px 0;
        color: var(--app-ink);
        font-size: 14px;
        font-weight: 900;
      }
      .mediation-effect-card code {
        color: var(--app-accent);
        background: transparent;
        padding: 0;
        font-size: 12px;
        font-weight: 900;
      }
      .mediation-effect-main {
        display: flex;
        justify-content: space-between;
        gap: 10px;
        align-items: center;
        border-radius: 8px;
        background: #f5f8ff;
        color: var(--app-muted);
        font-size: 12px;
        font-weight: 800;
        padding: 8px 9px;
        font-variant-numeric: tabular-nums;
      }
      .mediation-effect-value {
        color: var(--app-accent);
        font-size: 15px;
        font-weight: 900;
      }
      .mediation-interval-guide {
        display: inline-flex;
        align-items: center;
        gap: 6px;
        margin-top: 10px;
        color: var(--app-muted);
        font-size: 12px;
        font-weight: 700;
      }
      .sensitivity-dag-canvas {
        position: relative;
        min-height: 280px;
        border: 1px solid #e2eaf4;
        border-radius: 10px;
        background: linear-gradient(180deg, #ffffff 0%, #f8fbff 100%);
        overflow: hidden;
        margin-bottom: 14px;
      }
      .sensitivity-dag-canvas .mediation-node {
        width: 132px;
        min-height: 70px;
        padding: 12px 14px;
      }
      .sensitivity-dag-canvas .mediation-node-x,
      .sensitivity-dag-canvas .mediation-node-y {
        bottom: 18px;
      }
      .sensitivity-dag-canvas .mediation-node-m {
        top: 14px;
        width: min(230px, 46%);
        min-height: 70px;
      }
      .sensitivity-dag-canvas .mediation-arrow-label {
        width: 170px;
        min-height: 64px;
      }
      .sensitivity-dag-label {
        position: absolute;
        z-index: 3;
        min-width: 86px;
        border: 1px solid #dfe7f0;
        border-radius: 10px;
        background: rgba(255, 255, 255, 0.94);
        color: var(--app-ink);
        font-size: 14px;
        font-weight: 900;
        line-height: 1.15;
        padding: 8px 10px;
        box-shadow: 0 8px 20px rgba(31, 42, 68, 0.08);
        text-align: center;
      }
      .sensitivity-dag-label-a {
        left: 18%;
        top: 37%;
      }
      .sensitivity-dag-label-b {
        right: 18%;
        top: 37%;
      }
      .sensitivity-dag-label-c {
        left: 50%;
        bottom: 18px;
        transform: translateX(-50%);
      }
      .sensitivity-rho-box {
        position: absolute;
        z-index: 3;
        right: 31%;
        top: 30%;
        border: 1px solid rgba(47, 85, 151, 0.25);
        border-radius: 10px;
        background: #f3f7ff;
        color: var(--app-accent);
        font-size: 13px;
        font-weight: 900;
        padding: 8px 12px;
        box-shadow: 0 8px 18px rgba(47, 85, 151, 0.08);
      }
      .sensitivity-dag-note {
        margin: -2px 0 12px;
        color: var(--app-muted);
        font-size: 13px;
      }
      .mediation-path-mini-list {
        display: grid;
        gap: 8px;
      }
      .mediation-coeff-list {
        display: grid;
        gap: 6px;
        margin-top: 10px;
      }
      .mediation-coeff-row {
        display: grid;
        grid-template-columns: minmax(42px, 0.65fr) repeat(3, minmax(56px, 1fr));
        gap: 5px;
        align-items: center;
        border: 1px solid #e2eaf4;
        border-radius: 8px;
        background: #f9fbff;
        color: var(--app-ink);
        font-size: 12px;
        font-variant-numeric: tabular-nums;
        padding: 7px 8px;
      }
      .mediation-coeff-row:hover {
        border-color: #c7d8f2;
        background: #f3f7ff;
      }
      .mediation-coeff-name {
        font-weight: 800;
      }
      .mediation-coeff-value {
        white-space: nowrap;
      }
      .mediation-coeff-indirect {
        color: var(--app-accent);
        font-weight: 800;
      }
      .mediation-path-mini {
        border: 1px solid #e3ebf4;
        border-radius: 10px;
        background: #ffffff;
        padding: 10px;
      }
      .mediation-path-mini-head {
        display: flex;
        justify-content: space-between;
        gap: 8px;
        color: var(--app-ink);
        font-weight: 800;
      }
      .mediation-path-mini-effect {
        color: var(--app-accent);
        font-variant-numeric: tabular-nums;
      }
      .mediation-path-mini-detail {
        margin-top: 5px;
        color: var(--app-muted);
        font-size: 12px;
        font-variant-numeric: tabular-nums;
      }
      .mediation-view-tabs .nav-tabs {
        border-bottom: 1px solid var(--app-border);
        margin-bottom: 14px;
      }
	      .mediation-view-tabs .nav-tabs > li > a {
	        color: var(--app-muted);
	        font-weight: 700;
	        border-radius: 8px 8px 0 0;
	        border: 1px solid transparent;
	        padding: 9px 14px;
	      }
	      .mediation-view-tabs .nav-tabs > li > a:hover,
	      .mediation-view-tabs .nav-tabs > li > a:focus {
	        color: var(--app-accent);
	        background: var(--app-blue-soft);
	      }
	      .mediation-view-tabs .nav-tabs > li.active > a,
	      .mediation-view-tabs .nav-tabs > li.active > a:focus,
	      .mediation-view-tabs .nav-tabs > li.active > a:hover {
	        color: var(--app-accent);
	        background: var(--app-accent-soft);
	        border-color: #d9e3ff #d9e3ff var(--app-accent-soft) #d9e3ff;
	      }
      .mediation-view-tabs .tab-content {
        min-width: 0;
      }
      .mediation-path-table {
        width: 100%;
        border-collapse: collapse;
        font-size: 13px;
      }
      .mediation-path-table th {
        color: var(--app-muted);
        font-weight: 700;
        border-bottom: 1px solid #d9e2ec;
        padding: 8px 8px;
        text-align: left;
        white-space: nowrap;
      }
      .mediation-path-table td {
        border-bottom: 1px solid #eef3f8;
        padding: 8px 8px;
        vertical-align: top;
      }
      .mediation-path-table tr:last-child td {
        border-bottom: 0;
      }
      .mediation-header-actions {
        display: flex;
        justify-content: flex-end;
        align-items: center;
        gap: 8px;
        flex-shrink: 0;
      }
      .mediation-empty-action {
        margin-bottom: 14px;
      }
      .mediation-run-strip {
        display: flex;
        flex-direction: column;
        align-items: flex-start;
        gap: 10px;
        margin-bottom: 14px;
      }
      .mediation-run-action {
        margin-bottom: 20px;
      }
      .mediation-run-action .workflow-action {
        margin-top: 0;
      }
      .bayes-client-progress {
        background: var(--app-soft);
        border: 1px solid #d6e6ff;
        border-radius: 8px;
        margin-top: 12px;
        padding: 10px;
      }
      .bayes-client-progress-head {
        display: flex;
        justify-content: space-between;
        gap: 12px;
        color: var(--app-ink);
        font-size: 12px;
        font-weight: 700;
        margin-bottom: 8px;
      }
      .bayes-client-progress-track {
        background: #e4ecf7;
        border-radius: 999px;
        height: 8px;
        overflow: hidden;
      }
      .bayes-client-progress-fill {
        background: linear-gradient(90deg, var(--app-accent), var(--app-blue));
        border-radius: 999px;
        height: 100%;
        transition: width 900ms ease;
        width: 0%;
      }
      .bayes-client-progress-complete .bayes-client-progress-fill {
        background: linear-gradient(90deg, #1b6e3c, #2f80ed);
      }
      .bayes-client-progress-error .bayes-client-progress-fill {
        background: var(--app-red);
      }
      .sensitivity-shell {
        display: grid;
        grid-template-columns: minmax(280px, 340px) minmax(0, 1fr);
        gap: 18px;
        align-items: start;
        margin-top: 12px;
      }
      .sensitivity-main {
        min-width: 0;
      }
      .sensitivity-header {
        display: flex;
        justify-content: space-between;
        align-items: flex-start;
        gap: 14px;
        margin-bottom: 12px;
      }
	      .sensitivity-header h3 {
        margin: 0 0 4px 0;
        font-size: 20px;
        font-weight: 700;
	        color: var(--app-ink);
      }
      .sensitivity-header p {
        margin: 0;
      }
      .sensitivity-context {
        color: var(--app-muted);
        font-size: 13px;
        text-align: right;
      }
      .sensitivity-header-actions {
        display: flex;
        flex-direction: column;
        align-items: flex-end;
        gap: 8px;
      }
      .sensitivity-result-grid {
        display: grid;
        grid-template-columns: minmax(0, 1fr);
        gap: 14px;
      }
      .mediation-primary-tabs > .tabbable > .nav-tabs {
        display: flex;
        flex-wrap: wrap;
        gap: 4px;
        border-bottom: 2px solid var(--app-border);
        margin-bottom: 16px;
      }
      .mediation-primary-tabs > .tabbable > .nav-tabs > li > a {
        color: var(--app-muted);
        font-size: 15px;
        font-weight: 800;
        padding: 10px 16px;
        border: 1px solid transparent;
        border-bottom: 3px solid transparent;
      }
      .mediation-primary-tabs > .tabbable > .nav-tabs > li.active > a,
      .mediation-primary-tabs > .tabbable > .nav-tabs > li.active > a:hover,
      .mediation-primary-tabs > .tabbable > .nav-tabs > li.active > a:focus {
        color: var(--app-accent);
        background: #ffffff;
        border-color: transparent transparent var(--app-accent);
      }
	      .sensitivity-view-tabs .nav-tabs {
        border-bottom: 1px solid var(--app-border);
        margin-bottom: 14px;
        display: flex;
        flex-wrap: wrap;
      }
	      .sensitivity-view-tabs .nav-tabs > li > a {
	        color: var(--app-muted);
	        font-size: 13px;
	        font-weight: 700;
	        border-radius: 8px 8px 0 0;
	        border: 1px solid transparent;
	        padding: 9px 14px;
	      }
	      .sensitivity-view-tabs .nav-tabs > li > a:hover,
	      .sensitivity-view-tabs .nav-tabs > li > a:focus {
	        color: var(--app-violet);
	        background: var(--app-violet-soft);
	      }
	      .sensitivity-view-tabs .nav-tabs > li.active > a,
	      .sensitivity-view-tabs .nav-tabs > li.active > a:focus,
	      .sensitivity-view-tabs .nav-tabs > li.active > a:hover {
	        color: var(--app-violet);
	        background: var(--app-violet-soft);
	        border-color: #e3daf8 #e3daf8 var(--app-violet-soft) #e3daf8;
	      }
      .sensitivity-view-tabs .tab-content {
        min-width: 0;
      }
      .sensitivity-tipping-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(210px, 1fr));
        gap: 10px;
        margin-bottom: 12px;
      }
	      .sensitivity-tipping-card {
	        background: var(--app-slate-soft);
	        border: 1px solid #dce5f0;
	        border-top: 3px solid var(--app-slate-accent);
	        border-radius: 8px;
	        padding: 11px 12px;
	      }
      .sensitivity-tipping-label {
        color: var(--app-muted);
        font-size: 12px;
        margin-bottom: 5px;
      }
	      .sensitivity-tipping-value {
	        color: var(--app-accent);
	        font-size: 20px;
        font-weight: 700;
        font-variant-numeric: tabular-nums;
      }
      .sensitivity-decomp-panel {
        margin-bottom: 14px;
      }
	      .sensitivity-decomp-panel h4 {
        margin-top: 0;
        margin-bottom: 10px;
        font-size: 15px;
        font-weight: 700;
	        color: var(--app-ink);
      }
	      .sensitivity-effect-grid {
	        display: grid;
	        grid-template-columns: repeat(4, minmax(0, 1fr));
	        gap: 12px;
	        align-items: stretch;
	        margin-bottom: 12px;
	      }
      .sensitivity-mini-table {
        width: 100%;
        border-collapse: collapse;
        font-size: 13px;
      }
      .sensitivity-mini-table th {
        color: var(--app-muted);
        border-bottom: 1px solid #d9e2ec;
        padding: 7px 8px;
        text-align: left;
      }
      .sensitivity-mini-table td {
        border-bottom: 1px solid #eef3f8;
        padding: 7px 8px;
      }
      .sensitivity-mini-table tr:last-child td {
        border-bottom: 0;
      }
	      .sensitivity-warning {
	        color: #8a5a00;
	      }
      @media (max-width: 900px) {
        .setup-shell {
          grid-template-columns: 1fr;
        }
        .setup-sidebar {
          position: static;
        }
        .setup-main-grid {
          grid-template-columns: 1fr;
        }
        .overview-shell {
          grid-template-columns: 1fr;
        }
        .tutorial-workflow,
        .tutorial-tile-grid,
        .setup-file-grid,
        .setup-metric-grid,
        .setup-why-card,
        .setup-example-grid,
        .setup-question-flow,
        .setup-match-cue,
        .setup-compact-file-grid,
        .setup-role-grid,
        .setup-collapsed-grid {
          grid-template-columns: 1fr;
        }
        .setup-question-arrow {
          transform: rotate(90deg);
        }
        .spatial-reference-popover {
          margin-left: 0;
          width: 100%;
        }
        .spatial-reference-trigger {
          justify-content: center;
          width: 100%;
        }
        .spatial-reference-panel {
          left: 0;
          right: auto;
          width: min(100%, calc(100vw - 32px));
        }
        .spatial-reference-panel::before {
          left: 28px;
          right: auto;
        }
        .spatial-reference-table,
        .spatial-reference-table thead,
        .spatial-reference-table tbody,
        .spatial-reference-table tr,
        .spatial-reference-table th,
        .spatial-reference-table td {
          display: block;
        }
        .spatial-reference-table thead {
          display: none;
        }
        .spatial-reference-table td {
          border-top: 1px solid #eef3f8;
        }
        .spatial-reference-table td::before {
          content: attr(data-label);
          display: block;
          color: var(--app-ink);
          font-size: 11px;
          font-weight: 900;
          margin-bottom: 3px;
          text-transform: uppercase;
        }
        .setup-stage-header {
          display: block;
        }
        .setup-stage-actions {
          justify-content: flex-start;
          margin-top: 10px;
        }
        .setup-match-lines {
          text-align: left;
        }
	        .tutorial-hero {
	          padding: 18px;
	        }
	        .tutorial-hero h2 {
	          font-size: 23px;
	        }
	        .overview-plot-grid,
	        .overview-metric-row {
          grid-template-columns: 1fr;
        }
        .pca-shell,
        .pca-secondary-grid,
        .pca-loading-controls,
        .pca-detail-grid,
        .loading-result-grid,
        .loading-result-grid.single {
          grid-template-columns: 1fr;
        }
        .pca-map-grid.pc-count-1 {
          max-width: none;
        }
        .mediation-shell,
        .mediation-dag-layout,
        .effect-card-grid {
          grid-template-columns: 1fr;
        }
        .mediation-dag-canvas {
          min-height: 380px;
        }
        .mediation-node {
          width: 112px;
        }
        .mediation-node-m {
          width: min(230px, 58%);
        }
        .mediation-node-value {
          font-size: 16px;
        }
        .mediation-dag-summary-strip,
        .mediation-effect-grid,
        .mediation-coeff-row {
          grid-template-columns: 1fr;
        }
        .mediation-arrow-label {
          width: 128px;
          max-width: 128px;
          font-size: 11px;
          padding: 5px 7px;
        }
        .mediation-header {
          display: block;
        }
        .mediation-header-actions {
          justify-content: flex-start;
          margin-top: 8px;
        }
        .mediation-context {
          text-align: left;
          margin-top: 6px;
        }
        .sensitivity-shell,
        .sensitivity-effect-grid {
          grid-template-columns: 1fr;
        }
        .sensitivity-header {
          display: block;
        }
        .sensitivity-context {
          text-align: left;
          margin-top: 6px;
        }
        .sensitivity-header-actions {
          align-items: flex-start;
        }
      }
      "))
	    )
	  ),
	  tabPanel(
	    title = "Step 1 Data Preparation",
	    value = "Data Setup",
    div(
      class = "setup-shell",
      tags$aside(
        class = "setup-sidebar",
        uiOutput("data_setup_sidebar")
	      ),
	      tags$main(
	        class = "setup-main",
	        uiOutput("data_setup_main_content")
	      )
	    )
  ),
	  tabPanel(
	    title = "Step 2 Data Visualization",
	    value = "Overview",
	    div(
	      class = "overview-shell",
      tags$aside(
        class = "setup-sidebar",
        tags$div(
          class = "setup-sidebar-section setup-sidebar-actions",
          tags$div(class = "workflow-action",
            actionButton("go_to_pca_from_overview", "Continue to Feature Summarization", class = "btn-primary")
          )
        ),
        app_setup_section("Model Setup", uiOutput("overview_variables"))
      ),
	      tags$main(
        class = "overview-main",
        app_panel_card(
          "Spatial Overview",
          class = "overview-spatial-card",
          div(
            class = "overview-plot-grid",
            div(class = "overview-plot-panel", plotOutput("overview_exposure_spatial", height = 410)),
            div(class = "overview-plot-panel", plotOutput("overview_outcome_spatial", height = 410))
          )
        )
      )
    )
  ),
  tabPanel(
    title = "Step 3 Feature Summarization",
    value = "PCA",
    div(
      class = "pca-shell",
      tags$aside(
        class = "setup-sidebar",
        tags$div(
          class = "setup-sidebar-section setup-sidebar-actions",
          uiOutput("pca_confirmation_controls")
        ),
        app_setup_section(
          tagList(
            "PCs for Mediation",
            app_help_tooltip("Principal component analysis (PCA) turns the mediator feature matrix into summary axes called PCs. The selected PCs are used as mediators to study indirect pathways between exposure and outcome.")
          ),
          uiOutput("pca_pc_controls"),
          uiOutput("pca_selection_notice")
        )
      ),
	      tags$main(
	        class = "pca-main",
	        div(
          class = "pca-view-tabs",
          tabsetPanel(
            id = "pca_view",
            type = "tabs",
            tabPanel(
              "Summary Axis Details",
              div(
                class = "pca-loading-controls",
                uiOutput("loading_pc_control"),
                selectInput("loading_n", "Number of loading features", choices = c(10, 30, 50), selected = 10)
              ),
              div(
                class = "pca-detail-grid",
                app_panel_card(
                  "PC Score Map",
                  class = "pca-detail-map-card",
                  uiOutput("pca_detail_meta"),
                  plotOutput("pca_detail_spatial", height = 430)
                ),
                app_panel_card(
                  tagList(
                    "Loadings",
                    app_help_tooltip("Loadings show which input variables contribute most strongly to the selected PC. Larger absolute values indicate stronger contribution to that PC pattern.")
                  ),
                  uiOutput("loading_results")
                )
              )
            ),
            tabPanel(
              "Explained Variance",
              app_panel_card(
                "Explained Variance",
                tags$p(
                  class = "setup-meta",
                  "This shows how much of the overall gene or feature variation is summarized by each PC. Taller bars mean that PC captures more variation; the line shows how much variation is captured as PCs are added together."
                ),
                plotOutput("scree_plot", height = 430)
              )
            )
          )
        )
      )
    )
	  ),
	  tabPanel(
	    title = "Step 4 Mediation",
	    value = "Mediation",
	    div(
	      class = "mediation-shell",
	      tags$aside(
	        class = "setup-sidebar",
	        uiOutput("result_export_controls"),
	        uiOutput("mediation_run_controls"),
	        uiOutput("mediation_input_context"),
		        app_setup_section(
		          app_label_help(
		            "Framework",
		            "Bayesian analysis summarizes uncertainty with posterior distributions and credible intervals. Frequentist analysis gives point estimates and bootstrap confidence intervals and is usually faster, so it can be useful for quick comparison. Bayesian is the default."
		          ),
		          radioButtons(
		            "analysis_framework",
		            NULL,
		            choices = c("Bayesian" = "bayesian", "Frequentist" = "frequentist"),
		            selected = "bayesian"
		          )
		        ),
		        conditionalPanel(
		          condition = "input.mediation_primary_view == 'mediation'",
		          uiOutput("mediation_view_controls")
		        ),
	        conditionalPanel(
	          condition = "input.analysis_framework == 'frequentist'",
	          app_setup_section(
		            "Bootstrap CI",
		            numericInput("bootstrap_B", "Bootstrap replicates", value = 1000, min = 100, max = 5000, step = 100),
		            uiOutput("bootstrap_status")
	          )
	        ),
	        conditionalPanel(
	          condition = "input.mediation_primary_view == 'sensitivity'",
	          app_setup_section(
	            tagList(
	              "Mediator-outcome confounding",
	              app_help_tooltip("Sensitivity analysis varies structural mediator-outcome residual correlation (ρ), preserving the fitted conditional outcome variance and total effect. Common correlation varies all mediator-outcome correlations together; One mediator varies one correlation at a time. Correlated mediators can transmit either adjustment to multiple indirect effects.")
	            ),
	            radioButtons(
	              "sensitivity_scenario",
	              NULL,
	              choices = c("Common correlation" = "common", "One mediator" = "one_at_a_time"),
	              selected = "common"
	            ),
	            uiOutput("sensitivity_contrast_control"),
	            uiOutput("sensitivity_mediator_control")
	          )
	        )
	      ),
	      tags$main(
	        class = "mediation-main",
	        div(
	          class = "mediation-primary-tabs",
	          tabsetPanel(
	            id = "mediation_primary_view",
	            type = "tabs",
	            tabPanel("Mediation Results", value = "mediation", uiOutput("mediation_main_content")),
	            tabPanel(
	              "Sensitivity Analysis",
	              value = "sensitivity",
	              div(class = "sensitivity-view-tabs", uiOutput("sensitivity_workflow_tabs"))
	            )
	          )
	        )
	      )
	    )
	  )
)

server <- function(input, output, session) {
  analysis_state <- reactiveVal(NULL)
  setup_status <- reactiveVal(list(
    ok = FALSE,
    error = FALSE,
    message = ""
  ))
  data_setup_page <- reactiveVal("upload")
  bayesian_fit_state <- reactiveVal(NULL)
  bayesian_fit_cache <- reactiveVal(list())
  bayesian_fit_running <- reactiveVal(FALSE)
	  frequentist_fit_state <- reactiveVal(NULL)
	  frequentist_fit_cache <- reactiveVal(list())
	  frequentist_bootstrap_cache <- reactiveVal(list())
	  frequentist_sensitivity_cache <- reactiveVal(list())
  frequentist_fit_running <- reactiveVal(FALSE)
  sensitivity_warmup_status <- reactiveVal(NULL)
  confirmed_pca_spec <- reactiveVal(NULL)

  clear_workflow_results <- function() {
    confirmed_pca_spec(NULL)
    bayesian_fit_state(NULL)
    bayesian_fit_cache(list())
    frequentist_fit_state(NULL)
    frequentist_fit_cache(list())
    frequentist_bootstrap_cache(list())
    frequentist_sensitivity_cache(list())
    sensitivity_warmup_status(NULL)
  }

  set_analysis_unprepared <- function(message = "Upload both datasets, confirm variable roles, then prepare the analysis.", notify = FALSE) {
    clear_workflow_results()
    analysis_state(NULL)
    data_setup_page("upload")
    setup_status(list(ok = FALSE, error = FALSE, message = message))
    updateTabsetPanel(session, "main_nav", selected = "Data Setup")
    if (isTRUE(notify)) {
      showNotification(message, type = "message")
    }
  }

  safe_read_table <- function(file, sheet) {
    if (is.null(file)) {
      return(list(data = NULL, error = "No file uploaded."))
    }
    tryCatch(
      list(data = read_uploaded_table(file$datapath, sheet = sheet), error = NULL),
      error = function(e) list(data = NULL, error = conditionMessage(e))
    )
  }

  uploaded_main_cache <- reactiveVal(list(data = NULL, error = "No file uploaded.", defaults = list()))
  uploaded_feature_cache <- reactiveVal(list(data = NULL, error = "No file uploaded.", defaults = list()))

  output$main_sheet_ui <- renderUI({
    req(input$main_file)
    ext <- tolower(tools::file_ext(input$main_file$name))
    if (!ext %in% c("xlsx", "xls")) {
      return(NULL)
    }
    sheets <- tryCatch(available_excel_sheets(input$main_file$datapath), error = function(e) character(0))
    if (length(sheets) == 0) {
      return(helpText("Excel sheet list is unavailable."))
    }
    selectInput("main_sheet", "Main dataset sheet", choices = sheets, selected = sheets[[1]])
  })

  output$feature_sheet_ui <- renderUI({
    req(input$feature_file)
    ext <- tolower(tools::file_ext(input$feature_file$name))
    if (!ext %in% c("xlsx", "xls")) {
      return(NULL)
    }
    sheets <- tryCatch(available_excel_sheets(input$feature_file$datapath), error = function(e) character(0))
    if (length(sheets) == 0) {
      return(helpText("Excel sheet list is unavailable."))
    }
    selectInput("feature_sheet", "Feature matrix sheet", choices = sheets, selected = sheets[[1]])
  })

  observeEvent(list(input$main_file, input$main_sheet), {
    if (identical(input$data_source %||% "example", "upload") && !is.null(input$main_file)) {
      x <- safe_read_table(input$main_file, input$main_sheet)
      x$defaults <- list()
      uploaded_main_cache(x)
    }
  }, ignoreInit = TRUE)

  observeEvent(list(input$feature_file, input$feature_sheet), {
    if (identical(input$data_source %||% "example", "upload") && !is.null(input$feature_file)) {
      x <- safe_read_table(input$feature_file, input$feature_sheet)
      x$defaults <- list()
      uploaded_feature_cache(x)
    }
  }, ignoreInit = TRUE)

  template_main_data <- function() {
    data.frame(
      spot_id = sprintf("spot_%02d", 1:8),
      exposure = rep(c("control", "adjacent", "far", "control"), length.out = 8),
      outcome = c(1.12, 1.35, 1.72, 1.08, 1.44, 1.81, 1.20, 1.62),
      x_coord = c(10, 18, 27, 36, 12, 24, 33, 42),
      y_coord = c(8, 9, 11, 13, 20, 21, 23, 25),
      age = c(61, 58, 64, 55, 67, 60, 63, 59),
      region = rep(c("tumor", "margin"), length.out = 8),
      check.names = FALSE
    )
  }

  template_feature_data <- function() {
    data.frame(
      spot_id = sprintf("spot_%02d", 1:8),
      Gene1 = c(5.1, 5.3, 6.2, 5.0, 6.1, 6.7, 5.4, 6.4),
      Gene2 = c(2.4, 2.9, 3.2, 2.2, 3.0, 3.5, 2.7, 3.4),
      Gene3 = c(8.2, 8.0, 7.4, 8.5, 7.6, 7.1, 8.1, 7.3),
      Gene4 = c(1.1, 1.4, 1.8, 1.0, 1.7, 2.0, 1.3, 1.9),
      check.names = FALSE
    )
  }

  illustrative_main_data <- function() {
    data.frame(
      observation_id = c("spot_001", "spot_002", "spot_003"),
      exposure = c("inside", "adjacent", "far"),
      outcome = c(0.42, 0.31, 0.58),
      x_coord = c(5754, 4982, 6221),
      y_coord = c(9979, 7819, 8540),
      covariate = c("...", "...", "..."),
      check.names = FALSE
    )
  }

  illustrative_feature_data <- function() {
    data.frame(
      observation_id = c("spot_001", "spot_002", "spot_003"),
      Gene_A = c(1.21, 0.37, 1.04),
      Gene_B = c(-0.53, 1.04, 0.18),
      Gene_C = c(0.82, -0.12, -0.46),
      Gene_D = c(0.44, 0.91, -0.08),
      check.names = FALSE
    )
  }

  output$download_main_template <- downloadHandler(
    filename = function() "main_dataset_template.csv",
    content = function(file) {
      utils::write.csv(template_main_data(), file, row.names = FALSE)
    }
  )

  output$download_feature_template <- downloadHandler(
    filename = function() "mediator_feature_matrix_template.csv",
    content = function(file) {
      utils::write.csv(template_feature_data(), file, row.names = FALSE)
    }
  )

  upload_data_ready <- reactive({
    source <- input$data_source %||% "example"
    if (identical(source, "example")) {
      return(TRUE)
    }
    main <- uploaded_main_cache()
    feature <- uploaded_feature_cache()
    !is.null(main$data) && is.null(main$error) &&
      !is.null(feature$data) && is.null(feature$error)
  })

  observeEvent(input$confirm_data_upload, {
    if (!isTRUE(upload_data_ready())) {
      showNotification("Upload both readable files before confirming data upload.", type = "warning")
      return()
    }
    data_setup_page("review")
  }, ignoreInit = TRUE)

  observeEvent(input$back_to_upload, {
    data_setup_page("upload")
  }, ignoreInit = TRUE)

  output$confirm_data_upload_controls <- renderUI({
    ready <- isTRUE(upload_data_ready())
    button_args <- list(
      inputId = "confirm_data_upload",
      label = "Confirm Data Upload",
      class = "btn-primary"
    )
    if (!ready) {
      button_args$disabled <- "disabled"
    }

    tagList(
      do.call(actionButton, button_args),
      if (identical(input$data_source %||% "example", "upload") && !ready) {
        tags$div(class = "setup-status setup-meta", "Upload both readable files to continue.")
      }
    )
  })

  output$data_setup_sidebar <- renderUI({
    source <- input$data_source %||% "example"
    page <- data_setup_page()

    if (identical(page, "review")) {
      return(tagList(
        tags$div(
          class = "setup-sidebar-section setup-sidebar-actions",
          tags$div(class = "workflow-action", uiOutput("data_setup_prepare_top")),
          tags$div(class = "workflow-action", actionButton("back_to_upload", "Back to Upload", class = "btn-default"))
        ),
        app_setup_section(
          tagList(
            "Variable Roles",
            if (identical(source, "example")) {
              tagList(tags$span(class = "setup-ok", HTML("&#10003;")), tags$span(class = "setup-meta", " Automatically detected"))
            }
          ),
          uiOutput("observation_id_ui"),
          uiOutput("exposure_ui"),
          uiOutput("outcome_ui"),
          uiOutput("spatial_x_ui"),
          uiOutput("spatial_y_ui"),
          uiOutput("feature_id_ui"),
          uiOutput("covariates_ui"),
          uiOutput("reference_ui"),
          uiOutput("exposure_levels_ui")
        )
      ))
    }

    tagList(
      tags$div(
        class = "setup-sidebar-section setup-sidebar-actions",
        uiOutput("confirm_data_upload_controls")
      ),
      app_setup_section(
        "Data Source",
        radioButtons(
          "data_source",
          NULL,
          choiceNames = list(
            "Example Dataset",
            tagList(
              "Upload My Data",
              app_icon_only_tooltip(
                "glyphicon-folder-open",
                "CSV, TSV, TXT, XLSX, and XLS files are supported."
              )
            )
          ),
          choiceValues = c("example", "upload"),
          selected = source
        ),
        if (identical(source, "upload")) {
          tagList(
            tags$strong("1. Main dataset"),
            fileInput("main_file", NULL, accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
            uiOutput("main_sheet_ui"),
            tags$strong("2. Mediator feature matrix"),
            fileInput("feature_file", NULL, accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
            uiOutput("feature_sheet_ui")
          )
        } else {
          tags$p(class = "setup-meta", "Built-in spatial transcriptomics example with preset roles.")
        }
      ),
      app_setup_section(
        "Templates",
        tags$div(class = "workflow-action", downloadButton("download_main_template", "Download main template")),
        tags$div(class = "workflow-action", downloadButton("download_feature_template", "Download feature template")),
        spatial_platform_reference()
      )
    )
  })

  main_input <- reactive({
    if (identical(input$data_source %||% "example", "example")) {
      return(list(data = example_inputs$main_data, error = NULL, defaults = example_inputs$defaults))
    }
    uploaded_main_cache()
  })

  feature_input <- reactive({
    if (identical(input$data_source %||% "example", "example")) {
      return(list(data = example_inputs$feature_data, error = NULL, defaults = example_inputs$defaults))
    }
    uploaded_feature_cache()
  })

  observeEvent(input$data_source, {
    if (identical(input$data_source %||% "example", "example")) {
      set_analysis_unprepared("")
    } else {
      set_analysis_unprepared()
    }
  }, ignoreInit = TRUE)

  observeEvent(input$main_file, {
    if (identical(input$data_source %||% "example", "upload")) {
      set_analysis_unprepared("Main dataset selected. Confirm roles and prepare the analysis.")
    }
  }, ignoreInit = TRUE)

  observeEvent(input$feature_file, {
    if (identical(input$data_source %||% "example", "upload")) {
      set_analysis_unprepared("Feature matrix selected. Confirm matching and prepare the analysis.")
    }
  }, ignoreInit = TRUE)

  main_data_current <- reactive({
    x <- main_input()
    validate(need(is.null(x$error), x$error))
    x$data
  })

  feature_data_current <- reactive({
    x <- feature_input()
    validate(need(is.null(x$error), x$error))
    x$data
  })

  main_data_for_setup_ui <- reactive({
    x <- main_input()
    if (!is.null(x$error)) NULL else x$data
  })

  feature_data_for_setup_ui <- reactive({
    x <- feature_input()
    if (!is.null(x$error)) NULL else x$data
  })

  normalize_setup_name <- function(x) {
    tolower(gsub("[^a-z0-9]+", "", x))
  }

  setup_find_column <- function(data, aliases, predicate = NULL) {
    if (is.null(data)) {
      return("")
    }
    cols <- names(data)
    normalized <- normalize_setup_name(cols)
    alias_keys <- normalize_setup_name(aliases)
    matches <- cols[normalized %in% alias_keys]
    if (!is.null(predicate)) {
      matches <- matches[vapply(matches, function(col) isTRUE(predicate(data[[col]])), logical(1))]
    }
    if (length(matches) == 1) {
      matches[[1]]
    } else {
      ""
    }
  }

  setup_role_suggestions <- reactive({
    if (identical(input$data_source %||% "example", "example")) {
      return(example_inputs$defaults)
    }

    main_data <- main_data_for_setup_ui()
    feature_data <- feature_data_for_setup_ui()
    if (is.null(main_data)) {
      return(list())
    }

    is_numeric_col <- function(x) is.numeric(x)
    is_categorical_col <- function(x) {
      vals <- unique(stats::na.omit(as.character(x)))
      !is.numeric(x) && length(vals) >= 2 && length(vals) <= 20
    }

    id <- setup_find_column(main_data, c("spot_id", "barcode", "observation_id", "obs_id", "cell_id", "sample_id", "id"))
    exposure <- setup_find_column(main_data, c("treatment", "exposure", "condition", "group", "x_cat", "x"), is_categorical_col)
    if (!nzchar(exposure)) {
      candidates <- names(main_data)[vapply(main_data, is_categorical_col, logical(1))]
      candidates <- setdiff(candidates, id)
      exposure <- if (length(candidates) == 1) candidates[[1]] else ""
    }

    spatial_x <- setup_find_column(main_data, c("x_coord", "xcoord", "spatial_x", "imagecol", "x_coordinate", "x"), is_numeric_col)
    spatial_y <- setup_find_column(main_data, c("y_coord", "ycoord", "spatial_y", "imagerow", "y_coordinate", "y"), is_numeric_col)
    outcome <- setup_find_column(main_data, c("outcome", "response", "y_value", "y"), is_numeric_col)
    if (!nzchar(outcome)) {
      candidates <- names(main_data)[vapply(main_data, is_numeric_col, logical(1))]
      candidates <- setdiff(candidates, c(id, spatial_x, spatial_y))
      outcome <- if (length(candidates) == 1) candidates[[1]] else ""
    }

    feature_id <- ""
    if (!is.null(feature_data)) {
      feature_id <- setup_find_column(feature_data, c("spot_id", "barcode", "observation_id", "obs_id", "cell_id", "sample_id", "id"))
    }

    list(
      observation_id = id,
      exposure = exposure,
      outcome = outcome,
      spatial_x = spatial_x,
      spatial_y = spatial_y,
      covariates = character(0),
      reference = "",
      feature_observation_id = feature_id
    )
  })

  current_setup_spec <- reactive({
    main_data <- main_data_for_setup_ui()
    feature_data <- feature_data_for_setup_ui()
    if (is.null(main_data) || is.null(feature_data)) {
      return(NULL)
    }
    roles <- selected_setup_roles()
    if (is.null(roles) || !nzchar(roles$exposure %||% "") || !nzchar(roles$outcome %||% "") ||
        !nzchar(roles$spatial_x %||% "") || !nzchar(roles$spatial_y %||% "") ||
        !nzchar(roles$reference %||% "")) {
      return(NULL)
    }

    tryCatch(
      build_analysis_spec(
        data_source = input$data_source %||% "example",
        main_data = main_data,
        feature_data = feature_data,
        roles = roles,
        feature_id = selected_feature_id(),
        excluded_feature_columns = character(0)
      ),
      error = function(e) structure(list(message = conditionMessage(e)), class = "setup_spec_error")
    )
  })

  setup_spec_is_error <- function(x) {
    inherits(x, "setup_spec_error")
  }

  default_role <- function(name, fallback = "") {
    setup_role_suggestions()[[name]] %||% fallback
  }

  none_choices <- function(cols) {
    c("None" = "", cols)
  }

  required_choices <- function(cols) {
    if (identical(input$data_source %||% "example", "example")) {
      cols
    } else {
      c("Choose column" = "", cols)
    }
  }

  setup_input_or_default <- function(input_id, role_name, cols = NULL, fallback = "") {
    value <- input[[input_id]] %||% ""
    default <- default_role(role_name, fallback)
    if (is.null(value) || length(value) == 0 || !nzchar(value[[1]])) {
      value <- default
    }
    if (!is.null(cols) && nzchar(value[[1]]) && !value[[1]] %in% cols) {
      value <- default
    }
    value %||% fallback
  }

  selected_setup_roles <- reactive({
    main_data <- main_data_for_setup_ui()
    if (is.null(main_data)) {
      return(NULL)
    }
    main_cols <- names(main_data)
    roles <- list(
      observation_id = setup_input_or_default("observation_id", "observation_id", main_cols),
      exposure = setup_input_or_default("exposure_var", "exposure", main_cols),
      outcome = setup_input_or_default("outcome_var", "outcome", main_cols),
      spatial_x = setup_input_or_default("spatial_x_var", "spatial_x", main_cols),
      spatial_y = setup_input_or_default("spatial_y_var", "spatial_y", main_cols),
      covariates = input$covariates %||% default_role("covariates", character(0)),
      reference = input$reference_category %||% default_role("reference", "")
    )
    roles$covariates <- intersect(roles$covariates %||% character(0), main_cols)
    if (nzchar(roles$exposure %||% "") && roles$exposure %in% main_cols) {
      exposure_levels_now <- unique(as.character(main_data[[roles$exposure]]))
      default_reference <- default_role("reference", exposure_levels_now[[1]] %||% "")
      if (!nzchar(roles$reference %||% "") || !roles$reference %in% exposure_levels_now) {
        roles$reference <- if (default_reference %in% exposure_levels_now) default_reference else exposure_levels_now[[1]]
      }
    }
    roles
  })

  selected_feature_id <- reactive({
    feature_data <- feature_data_for_setup_ui()
    if (is.null(feature_data)) {
      return("")
    }
    setup_input_or_default("feature_observation_id", "feature_observation_id", names(feature_data))
  })

  output$main_data_summary <- renderTable({
    data <- main_data_current()
    numeric_cols <- numeric_column_names(data)
    data.frame(
      item = c("observations", "columns", "numeric columns"),
      value = c(nrow(data), ncol(data), length(numeric_cols)),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$feature_data_summary <- renderTable({
    data <- feature_data_current()
    numeric_cols <- numeric_column_names(data)
    data.frame(
      item = c("observations", "columns", "numeric columns"),
      value = c(nrow(data), ncol(data), length(numeric_cols)),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$main_data_preview <- renderDT({
    datatable(
      table_preview(main_data_current(), n_rows = 6, n_cols = 6),
      rownames = FALSE,
      class = "compact stripe nowrap",
      options = list(pageLength = 6, dom = "t", scrollX = TRUE, autoWidth = TRUE)
    )
  }, server = FALSE)

  output$feature_data_preview <- renderDT({
    datatable(
      table_preview(feature_data_current(), n_rows = 6, n_cols = 6),
      rownames = FALSE,
      class = "compact stripe nowrap",
      options = list(pageLength = 6, dom = "t", scrollX = TRUE, autoWidth = TRUE)
    )
  }, server = FALSE)

  output$main_column_types <- renderDT({
    datatable(column_type_table(main_data_current()), rownames = FALSE, options = list(pageLength = 8, dom = "tip"))
  }, server = FALSE)

  output$feature_numeric_status <- renderText({
    data <- feature_data_current()
    numeric_cols <- numeric_column_names(data)
    paste0(length(numeric_cols), " numeric columns detected before excluding ID/non-feature columns.")
  })

  output$observation_id_ui <- renderUI({
    data <- main_data_for_setup_ui()
    if (is.null(data)) return(NULL)
    cols <- names(data)
    selectInput(
      "observation_id",
      tagList("Observation ID", app_help_tooltip("Column that identifies the same spatial observation in both files. Leave as None only when both files are already in exactly the same row order.")),
      choices = none_choices(cols),
      selected = default_role("observation_id")
    )
  })

  output$exposure_ui <- renderUI({
    data <- main_data_for_setup_ui()
    if (is.null(data)) return(NULL)
    cols <- names(data)
    selectInput(
      "exposure_var",
      tagList("Exposure (X)", app_help_tooltip("Categorical group or condition being compared, such as exposure group, region, or tissue class.")),
      choices = required_choices(cols),
      selected = default_role("exposure", if (identical(input$data_source %||% "example", "example")) cols[[1]] else "")
    )
  })

  output$outcome_ui <- renderUI({
    data <- main_data_for_setup_ui()
    if (is.null(data)) return(NULL)
    cols <- names(data)
    selectInput(
      "outcome_var",
      tagList("Outcome (Y)", app_help_tooltip("Numeric response variable that the analysis tries to explain.")),
      choices = required_choices(cols),
      selected = default_role("outcome", if (identical(input$data_source %||% "example", "example")) cols[[1]] else "")
    )
  })

  output$spatial_x_ui <- renderUI({
    data <- main_data_for_setup_ui()
    if (is.null(data)) return(NULL)
    cols <- names(data)
    selectInput(
      "spatial_x_var",
      tagList("Spatial X Coordinate", app_help_tooltip("Numeric horizontal spatial coordinate used for maps and spatial adjustment.")),
      choices = required_choices(cols),
      selected = default_role("spatial_x", if (identical(input$data_source %||% "example", "example")) cols[[1]] else "")
    )
  })

  output$spatial_y_ui <- renderUI({
    data <- main_data_for_setup_ui()
    if (is.null(data)) return(NULL)
    cols <- names(data)
    selectInput(
      "spatial_y_var",
      tagList("Spatial Y Coordinate", app_help_tooltip("Numeric vertical spatial coordinate used for maps and spatial adjustment.")),
      choices = required_choices(cols),
      selected = default_role("spatial_y", if (identical(input$data_source %||% "example", "example")) cols[[1]] else "")
    )
  })

  output$covariates_ui <- renderUI({
    data <- main_data_for_setup_ui()
    if (is.null(data)) return(NULL)
    cols <- names(data)
    selected <- default_role("covariates", character(0))
    selectizeInput(
      "covariates",
      tagList("Additional Covariates", app_help_tooltip("Optional extra variables to adjust for. Leave blank if you do not need additional adjustment.")),
      choices = cols,
      selected = selected,
      multiple = TRUE
    )
  })

  output$feature_id_ui <- renderUI({
    data <- feature_data_for_setup_ui()
    if (is.null(data)) return(NULL)
    cols <- names(data)
    selectizeInput(
      "feature_observation_id",
      tagList("Feature Matrix Observation ID", app_help_tooltip("Column in the mediator feature matrix that matches the main dataset observation ID.")),
      choices = none_choices(cols),
      selected = default_role("feature_observation_id"),
      options = list(maxOptions = 2000)
    )
  })

  exposure_levels <- reactive({
    req(input$exposure_var)
    req(nzchar(input$exposure_var %||% ""))
    data <- main_data_for_setup_ui()
    req(data)
    unique(as.character(data[[input$exposure_var]]))
  })

  output$reference_ui <- renderUI({
    levels <- exposure_levels()
    if (length(levels) == 0) {
      return(NULL)
    }
    selected <- default_role("reference", levels[[1]])
    if (!selected %in% levels) {
      selected <- levels[[1]]
    }
    selectInput(
      "reference_category",
      tagList("Reference Category", app_help_tooltip("Baseline group. Contrasts compare other exposure groups against this group.")),
      choices = levels,
      selected = selected
    )
  })

  output$exposure_levels_ui <- renderUI({
    NULL
  })

  output$data_setup_role_controls <- renderUI({
    tagList(
      if (identical(input$data_source %||% "example", "example")) {
        tags$div(class = "setup-success-banner", tags$span(class = "setup-ok", HTML("&#10003;")), tags$strong("Automatically detected"))
      },
      tags$div(
        class = "setup-role-grid",
        uiOutput("observation_id_ui"),
        uiOutput("exposure_ui"),
        uiOutput("outcome_ui"),
        uiOutput("spatial_x_ui"),
        uiOutput("spatial_y_ui"),
        uiOutput("feature_id_ui"),
        uiOutput("covariates_ui"),
        uiOutput("reference_ui")
      ),
        uiOutput("exposure_levels_ui")
    )
  })

  output$data_setup_prepare_top <- renderUI({
    spec <- display_setup_spec()
    ready <- !is.null(spec) && !setup_spec_is_error(spec)
    button_args <- list(
      inputId = "prepare_analysis",
      label = "Confirm Setup and Continue",
      class = "btn-primary"
    )
    if (!isTRUE(ready)) {
      button_args$disabled <- "disabled"
    }
    tags$div(
      class = "setup-prepare-control",
      do.call(actionButton, button_args)
    )
  })

  observeEvent(input$prepare_analysis, {
    roles <- selected_setup_roles()

    result <- tryCatch({
      withProgress(message = "Preparing analysis data", value = 0, {
        setProgress(0.08, detail = "Checking selected variables and uploaded files")
        spec <- build_analysis_spec(
          data_source = input$data_source %||% "example",
          main_data = main_data_current(),
          feature_data = feature_data_current(),
          roles = roles,
          feature_id = selected_feature_id(),
          excluded_feature_columns = character(0)
        )
        setProgress(0.35, detail = "Matching observations across datasets")
        setProgress(0.55, detail = "Preparing PCA inputs from the feature matrix")
        analysis_obj <- compute_analysis_from_spec(spec)
        setProgress(0.95, detail = "Preparing app views")
        list(ok = TRUE, spec = spec, analysis = analysis_obj, message = "Analysis data prepared successfully.")
      })
    }, error = function(e) {
      list(ok = FALSE, message = conditionMessage(e))
    })

    if (isTRUE(result$ok)) {
      analysis_state(result$analysis)
      setup_status(list(ok = TRUE, error = FALSE, message = result$message))
	      confirmed_pca_spec(NULL)
	      bayesian_fit_state(NULL)
	      bayesian_fit_cache(list())
	      frequentist_fit_state(NULL)
	      frequentist_fit_cache(list())
	      frequentist_bootstrap_cache(list())
	      frequentist_sensitivity_cache(list())
	      sensitivity_warmup_status(NULL)
      updateTabsetPanel(session, "main_nav", selected = "Overview")
    } else {
      setup_status(list(ok = FALSE, error = TRUE, message = result$message))
    }
  }, ignoreInit = TRUE)

  output$prepare_message <- renderUI({
    if (identical(input$data_source %||% "example", "upload") &&
        (is.null(input$main_file) || is.null(input$feature_file))) {
      return(tags$div(class = "setup-status setup-meta", "Upload both datasets, then prepare the analysis."))
    }
    status <- setup_status()
    if (isTRUE(status$ok)) {
      return(NULL)
    }
    color <- if (isTRUE(status$ok)) {
      "#1b6e3c"
    } else if (isTRUE(status$error)) {
      "#9b1c1c"
    } else {
      "#5f6f7f"
    }
    if (!nzchar(status$message %||% "")) {
      return(NULL)
    }
    tags$div(class = "setup-status", style = paste0("color:", color, ";"), status$message)
  })

  output$validation_status <- renderTable({
    validation_status_table(analysis()$spec)
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  display_setup_spec <- reactive({
    spec <- current_setup_spec()
    if (is.null(spec) && identical(input$data_source %||% "example", "example")) {
      return(example_spec)
    }
    spec
  })

  spatial_reference_row <- function(examples, observation, features) {
    tags$tr(
      tags$td(`data-label` = "Examples", examples),
      tags$td(`data-label` = "One observation", observation),
      tags$td(`data-label` = "Features for mediator summarization", features)
    )
  }

  spatial_platform_reference <- function() {
    tags$div(
      class = "spatial-reference-popover",
      tags$button(
        type = "button",
        class = "spatial-reference-trigger",
        tags$span(class = "glyphicon glyphicon-info-sign", `aria-hidden` = "true"),
        tags$span("Working with spatial platform data?")
      ),
      tags$div(
        class = "spatial-reference-panel",
        role = "dialog",
        `aria-label` = "Spatial platform data reference",
        tags$div(
          class = "spatial-reference-header",
          tags$h5("Working with spatial platform data?"),
          tags$button(
            type = "button",
            class = "spatial-reference-close",
            `aria-label` = "Close spatial platform reference",
            HTML("&times;")
          )
        ),
        tags$p("Processed spatial data can take different forms; use the examples below to identify what corresponds to one analysis observation and the features used for mediator summarization."),
        tags$table(
          class = "spatial-reference-table",
          tags$thead(tags$tr(
            tags$th("Examples"),
            tags$th("One observation"),
            tags$th("Features for mediator summarization")
          )),
          tags$tbody(
            spatial_reference_row(
              "Visium / Visium HD · Curio Seeker · Stereo-seq · Slide-seq · DBiT-seq",
              "Spot / bin / bead",
              "Gene/transcript abundance measurements"
            ),
            spatial_reference_row(
              "Xenium · CosMx RNA · MERFISH / MERSCOPE · Molecular Cartography",
              "Segmented cell / entity",
              "Gene/transcript abundance measurements"
            ),
            spatial_reference_row(
              "CODEX / PhenoCycler · Hyperion IMC · MIBI",
              "Segmented cell",
              "Protein-marker measurements"
            ),
            spatial_reference_row(
              "GeoMx DSP",
              "Region / AOI / compartment",
              "Molecular measurements"
            )
          )
        ),
        tags$div(
          class = "spatial-reference-boundary",
          tags$span("Choose the features that represent the molecular profile relevant to your mediation question."),
          tags$span("Platform-specific preprocessing must be completed before upload. Native platform files are not imported directly."),
          tags$span("The app centers and scales the uploaded numeric features before summarizing features.")
        )
      )
    )
  }

  output$data_setup_main_content <- renderUI({
    source <- input$data_source %||% "example"
    page <- data_setup_page()

    preparation_guide <- app_panel_card(
      "Prepare your data",
      class = "setup-intro-card",
      tags$p(class = "setup-section-kicker", "1 · Define your mediation question"),
      tags$div(
        class = "setup-question-strip",
        tags$div(
          class = "setup-question-flow",
          tags$div(
            class = "setup-question-node",
            tags$strong("Exposure (X)"),
            tags$span("What differs?")
          ),
          tags$div(class = "setup-question-arrow", "->"),
          tags$div(
            class = "setup-question-node",
            tags$strong("Mediator features"),
            tags$span("What may mediate?"),
            tags$span(
              "Summarized into ",
              app_label_help("PCs", "PCA is a dimension-reduction method; PCs are the summary axes produced from the uploaded feature matrix.")
            )
          ),
          tags$div(class = "setup-question-arrow", "->"),
          tags$div(
            class = "setup-question-node",
            tags$strong("Outcome (Y)"),
            tags$span("What changes?")
          )
        ),
        tags$p(
          class = "setup-question-note",
          "Mediation asks how much of the relationship between exposure (X) and outcome (Y) operates through the mediator (M). Choose the exposure, outcome, and high-dimensional feature profile for your research question; the app summarizes those features into mediators."
        )
      ),
      tags$p(
        class = "setup-section-kicker",
        "2 · Prepare two files"
      ),
      tags$div(
        class = "setup-example-grid",
        tags$div(
          class = "setup-spreadsheet-example",
          tags$h5(
            setup_heading_icon(
              "1. Main dataset",
              "glyphicon-list-alt"
            )
          ),
          tags$p(
            class = "setup-card-note",
            "One row per observation, with the exposure group, outcome value, spatial coordinates, and optional covariates."
          ),
          setup_spreadsheet_table(illustrative_main_data()),
          tags$div(
            class = "setup-requirement-stack",
            setup_requirement_line(
              "Required: ",
              "exposure, outcome, x/y coordinates",
              tooltip = "Exposure must be categorical with at least two observed categories. Outcome and spatial coordinates must be numeric and contain no missing or non-finite values."
            ),
            tags$p(
              class = "setup-requirement-line",
              tags$strong("Recommended: "),
              "observation ID"
            ),
            setup_requirement_line("Optional: ", "covariates")
          ),
        ),
        tags$div(
          class = "setup-spreadsheet-example",
          tags$h5(
            setup_heading_icon(
              "2. Mediator feature matrix",
              "glyphicon-th"
            )
          ),
          tags$p(
            class = "setup-card-note",
            "High-dimensional features representing the profile of interest, measured on the same observations; these features are summarized into mediator variables."
          ),
          setup_spreadsheet_table(illustrative_feature_data()),
          tags$div(
            class = "setup-requirement-stack",
            setup_requirement_line(
              "Required: ",
              "numeric feature columns",
              tooltip = "At least three numeric feature columns are required. Feature values must be finite and must not contain missing values. Features with zero variance are not accepted. The app does not normalize or transform uploaded feature values, so features should be appropriately preprocessed before upload. Eligible numeric features are standardized before dimension reduction. If the outcome is derived from molecular features, consider whether those same features should also be included in the feature matrix used to construct the mediator PCs."
            ),
            tags$p(
              class = "setup-requirement-line",
              tags$strong("Recommended: "),
              "matching observation ID"
            )
          ),
        )
      ),
      tags$div(
        class = "setup-inline-help-row",
        tags$span(
          "Column names and order do not need to match these examples."
        )
      )
    )

    if (identical(page, "upload")) {
      return(preparation_guide)
    }

    main <- main_input()
    feature <- feature_input()

    read_errors <- c(main$error, feature$error)
    read_errors <- read_errors[!is.na(read_errors) & nzchar(read_errors)]
    if (length(read_errors) > 0) {
      return(tagList(
        app_panel_card(
          "Review & Configure",
          tags$p("The selected files could not be read. Return to Upload to replace them.")
        ),
        app_panel_card(
          "Data upload issue",
          class = "setup-empty",
          if (identical(source, "upload") && (is.null(input$main_file) || is.null(input$feature_file))) {
            tags$div(class = "setup-warning-banner", "Upload both the main dataset and mediator feature matrix to prepare custom data.")
          } else {
            tags$div(class = "setup-error-banner", paste(unique(read_errors), collapse = " "))
          }
        )
      ))
    }

    tagList(
      tags$section(
        class = "app-card",
        tags$div(
          class = "setup-stage-header",
          tags$h4("Review & Configure")
        ),
        tags$p(
          class = "setup-meta",
          "We suggested column roles automatically. Please review and change any incorrect selections before continuing."
        ),
        uiOutput("prepare_message")
      ),
      app_panel_card(
        "Data Preview",
        tags$div(
          class = "setup-main-grid",
          tags$div(
            class = "setup-preview-panel",
            tags$h5("Main Dataset"),
            tags$div(
              class = "setup-metric-grid",
              setup_metric_row("Observations", format(nrow(main$data), big.mark = ",")),
              setup_metric_row("Columns", format(ncol(main$data), big.mark = ","))
            ),
            setup_preview_table(main$data, n_rows = 4, n_cols = 6)
          ),
          tags$div(
            class = "setup-preview-panel",
            tags$h5("Mediator Feature Matrix"),
            tags$div(
              class = "setup-metric-grid",
              setup_metric_row("Observations", format(nrow(feature$data), big.mark = ",")),
              setup_metric_row("Numeric features", format(length(numeric_column_names(feature$data)), big.mark = ","))
            ),
            setup_preview_table(feature$data, n_rows = 4, n_cols = 6)
          )
        )
      ),
      app_panel_card(
        "Matching and Validation",
        uiOutput("validation_status_items")
      ),
      tags$details(
        class = "setup-wide-details",
        tags$summary("View validation and matching details"),
        tags$div(
          class = "setup-detail-body",
          uiOutput("matching_summary"),
          uiOutput("validation_detail_items")
        )
      )
    )
  })

  output$validation_status_items <- renderUI({
    spec <- display_setup_spec()
    if (is.null(spec)) {
      return(tags$p(class = "setup-meta", "Validation will appear after data files and variable roles are available."))
    }
    if (setup_spec_is_error(spec)) {
      return(tags$div(class = "setup-error-banner", spec$message))
    }

    ready_text <- if (identical(input$data_source %||% "example", "example")) {
      "Data ready for analysis."
    } else {
      "Data ready for analysis."
    }
    tagList(
      tags$div(class = "setup-success-banner", tags$span(class = "setup-ok", HTML("&#10003;")), tags$strong(ready_text)),
      tags$div(
        class = "setup-metric-grid",
        setup_metric_row("Matched observations", format(spec$matching$matched_n, big.mark = ",")),
        setup_metric_row("Mediator features", format(spec$mediator_features$n_features, big.mark = ",")),
        setup_metric_row("Exposure categories", paste(spec$exposure$levels, collapse = ", ")),
        setup_metric_row("Reference", spec$exposure$reference)
      ),
      if (!is.null(spec$matching$warning)) {
        tags$div(class = "setup-warning-banner", spec$matching$warning)
      }
    )
  })

  output$validation_detail_items <- renderUI({
    spec <- display_setup_spec()
    if (is.null(spec)) {
      return(tags$p(class = "setup-meta", "Validation details will appear after data files and variable roles are available."))
    }
    if (setup_spec_is_error(spec)) {
      return(tags$div(class = "setup-error-banner", spec$message))
    }

    rows <- validation_status_table(spec)
    tags$ul(
      class = "setup-validation-list",
      lapply(seq_len(nrow(rows)), function(i) {
        tags$li(
          tags$span(class = "setup-ok", HTML("&#10003;")),
          tags$strong(rows$Check[[i]]),
          ": ",
          rows$Status[[i]]
        )
      })
    )
  })

  output$matching_summary <- renderUI({
    spec <- display_setup_spec()
    if (is.null(spec)) {
      return(tags$p(class = "setup-meta", "Matching summary will appear after both datasets are available."))
    }
    if (setup_spec_is_error(spec)) {
      return(tags$p(class = "setup-meta", "Resolve the validation issue before matching can be summarized."))
    }

    unmatched <- (spec$matching$main_n - spec$matching$matched_n) +
      (spec$matching$feature_n - spec$matching$matched_n)
    tagList(
      tags$div(
        class = "setup-metric-grid",
        setup_metric_row("Main rows", format(spec$matching$main_n, big.mark = ",")),
        setup_metric_row("Feature rows", format(spec$matching$feature_n, big.mark = ",")),
        setup_metric_row("Matched rows", format(spec$matching$matched_n, big.mark = ",")),
        setup_metric_row("Unmatched rows", format(unmatched, big.mark = ","))
      ),
      if (!is.null(spec$matching$warning)) {
        tags$p(class = "setup-card-note", tags$span(class = "setup-warning", "!"), spec$matching$warning)
      } else {
        tags$p(class = "setup-card-note", "Observation IDs are being used for matching.")
      }
    )
  })

  output$matching_summary_table <- renderTable({
    spec <- display_setup_spec()
    validate(need(!is.null(spec) && !setup_spec_is_error(spec), ""))
    unmatched <- (spec$matching$main_n - spec$matching$matched_n) +
      (spec$matching$feature_n - spec$matching$matched_n)
    data.frame(
      item = c("Main dataset", "Feature matrix", "Matched", "Unmatched"),
      value = c(spec$matching$main_n, spec$matching$feature_n, spec$matching$matched_n, unmatched),
      check.names = FALSE
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  analysis <- reactive({
    analysis_state()
  })

  current_pca_spec <- reactive({
    selected <- selected_pc_ids(input, analysis())
    signs <- pc_signs_for_selection(selected, analysis())
    list(
      selected_pcs = selected,
      pc_signs = signs
    )
  })

  pca_specs_equal <- function(a, b) {
    if (is.null(a) || is.null(b)) {
      return(FALSE)
    }
    if (!identical(a$selected_pcs, b$selected_pcs)) {
      return(FALSE)
    }
    selected <- a$selected_pcs
    identical(
      unname(a$pc_signs[selected]),
      unname(b$pc_signs[selected])
    )
  }

  pca_config_confirmed <- reactive({
    !is.null(confirmed_pca_spec())
  })

  pca_has_unconfirmed_changes <- reactive({
    spec <- confirmed_pca_spec()
    !is.null(spec) && !pca_specs_equal(spec, current_pca_spec())
  })

  confirmed_selected_pc_ids <- reactive({
    spec <- confirmed_pca_spec()
    if (is.null(spec)) {
      return(character(0))
    }
    spec$selected_pcs
  })

  confirmed_pc_sign_vector <- reactive({
    spec <- confirmed_pca_spec()
    if (is.null(spec)) {
      return(pc_sign_vector(input, analysis()))
    }
    spec$pc_signs
  })

  mediation_sign_vector_confirmed <- reactive({
    pc_signs <- confirmed_pc_sign_vector()
    c(
      PC1_R = pc_signs[["PC1"]],
      PC2_R = pc_signs[["PC2"]],
      PC3 = pc_signs[["PC3"]]
    )
  })

  confirmed_pca_context <- reactive({
    spec <- confirmed_pca_spec()
    if (is.null(spec)) {
      return(NULL)
    }
    a <- analysis()
    list(
      exposure = a$spec$original_labels$exposure,
      reference = a$spec$exposure$reference,
      outcome = a$spec$original_labels$outcome,
      coordinates = paste(c(a$spec$original_labels$spatial_x, a$spec$original_labels$spatial_y), collapse = ", "),
      covariates = if (length(a$spec$covariates) == 0) "None" else paste(unname(a$spec$covariates), collapse = ", "),
      mediators = pca_spec_label(spec)
    )
  })

  sensitivity <- reactive({
    freq <- frequentist_current_result()
    validate(need(!is.null(freq), "Run Frequentist Mediation first for the current data and selected PCs."))
    key <- freq$key
    cache <- frequentist_sensitivity_cache()
    if (!is.null(cache[[key]])) {
      return(cache[[key]])
    }
    value <- compute_frequentist_sensitivity_with_progress(freq$result)
    cache[[key]] <- value
    frequentist_sensitivity_cache(cache)
    value
  })

  compute_frequentist_sensitivity_with_progress <- function(frequentist_result) {
    withProgress(message = "Preparing sensitivity analysis", value = 0, {
      setProgress(0.03, detail = "Starting residual-correlation sensitivity")
      result <- compute_frequentist_dynamic_sensitivity(
        frequentist_result = frequentist_result,
        rho_grid = bayesian_sensitivity_rho_grid,
        progress_callback = function(value, detail) {
          setProgress(value, detail = detail)
        }
      )
      setProgress(1, detail = "Sensitivity analysis ready")
      result
    })
  }

  set_sensitivity_warmup_status <- function(framework, state, message = NULL) {
    status <- sensitivity_warmup_status() %||% list()
    status[[framework]] <- list(state = state, message = message)
    sensitivity_warmup_status(status)
  }

  warm_bayesian_sensitivity_cache <- function(result, rho = 0) {
    if (is.null(result) || is.null(result$artifact)) {
      return(invisible(FALSE))
    }
    set_sensitivity_warmup_status("bayesian", "running")
    ok <- tryCatch({
      bayesian_sensitivity_curve_cache(
        artifact = result$artifact,
        artifact_key = result$key,
        mode = "common"
      )
      bayesian_sensitivity_selected_cache(
        artifact = result$artifact,
        artifact_key = result$key,
        rho = rho,
        mode = "common"
      )
      TRUE
    }, error = function(e) {
      set_sensitivity_warmup_status("bayesian", "error", conditionMessage(e))
      FALSE
    })
    if (isTRUE(ok)) {
      set_sensitivity_warmup_status("bayesian", "ready")
    }
    invisible(ok)
  }

  warm_frequentist_sensitivity_cache <- function(result) {
    if (is.null(result) || is.null(result$result)) {
      return(invisible(FALSE))
    }
    set_sensitivity_warmup_status("frequentist", "running")
    ok <- tryCatch({
      key <- result$key
      cache <- frequentist_sensitivity_cache()
      if (is.null(cache[[key]])) {
        cache[[key]] <- compute_frequentist_sensitivity_with_progress(result$result)
        frequentist_sensitivity_cache(cache)
      }
      TRUE
    }, error = function(e) {
      set_sensitivity_warmup_status("frequentist", "error", conditionMessage(e))
      FALSE
    })
    if (isTRUE(ok)) {
      set_sensitivity_warmup_status("frequentist", "ready")
    }
    invisible(ok)
  }

	  output$sensitivity_contrast_control <- renderUI({
	    a <- analysis()
	    contrast_choices <- stats::setNames(unname(a$contrast_labels), unname(a$contrast_labels))
    selected <- input$sensitivity_contrast_filter %||% unname(contrast_choices)[[1]]
    if (!selected %in% unname(contrast_choices)) {
      selected <- unname(contrast_choices)[[1]]
    }
    selectInput("sensitivity_contrast_filter", "Contrast", choices = contrast_choices, selected = selected)
  })

	  contrast_label_current <- function(exposure_names) {
    labels <- analysis()$contrast_labels
    out <- unname(labels[exposure_names])
    ifelse(is.na(out), exposure_names, out)
  }

  filter_sensitivity_contrast <- function(df) {
    contrast <- input$sensitivity_contrast_filter %||% "all"
    if (identical(contrast, "all") || !"contrast" %in% names(df)) {
      return(df)
    }
    dplyr::filter(df, contrast == !!contrast)
  }

  bayesian_artifact_available <- reactive({
    if (!isTRUE(pca_config_confirmed())) {
      return(FALSE)
    }
    entry <- bayesian_validated_artifact_for_selection(confirmed_selected_pc_ids())
    isTRUE(analysis()$bayesian_artifact_allowed) && !is.null(entry) &&
      isTRUE(all(confirmed_pc_sign_vector()[entry$pc_ids] == entry$pc_signs))
  })

  bayesian_current_mediators <- reactive({
    if (isTRUE(bayesian_artifact_available())) {
      entry <- bayesian_validated_artifact_for_selection(confirmed_selected_pc_ids())
      return(bayes_A_bundle_mediators(entry$artifact$sensitivity_A$residual_draws))
    }
    confirmed_selected_pc_ids()
  })

  bayesian_dynamic_settings <- reactive({
    bayesian_dynamic_default_settings()
  })

  bayesian_current_cache_key <- reactive({
    req(pca_config_confirmed())
    bayesian_dynamic_cache_key(
      analysis(),
      selected_pcs = confirmed_selected_pc_ids(),
      pc_signs = confirmed_pc_sign_vector(),
      settings = bayesian_dynamic_settings()
    )
  })

  bayesian_current_result <- reactive({
    if (!isTRUE(pca_config_confirmed())) {
      return(NULL)
    }
    if (isTRUE(bayesian_artifact_available())) {
      entry <- bayesian_validated_artifact_for_selection(confirmed_selected_pc_ids())
      return(list(
        source = "validated_artifact",
        key = paste("emory_validated", entry$key,
                    paste(entry$pc_signs, collapse = ","), structural_rho_version(), sep = "_"),
        artifact = entry$artifact
      ))
    }

    key <- bayesian_current_cache_key()
    state <- bayesian_fit_state()
    if (!is.null(state) && identical(state$key, key)) {
      return(state)
    }

    cached <- bayesian_fit_cache()[[key]]
    if (!is.null(cached)) {
      bayesian_fit_state(cached)
      return(cached)
    }

    NULL
  })

  bayesian_current_artifact <- reactive({
    result <- bayesian_current_result()
    if (is.null(result)) NULL else result$artifact
  })

  pc_signs_are_default <- reactive({
    if (!isTRUE(pca_config_confirmed())) {
      return(FALSE)
    }
    selected <- confirmed_selected_pc_ids()
    all(confirmed_pc_sign_vector()[selected] == pc_signs_for_selection(selected, analysis())[selected])
  })

  frequentist_example_default_available <- reactive({
    if (!isTRUE(pca_config_confirmed())) {
      return(FALSE)
    }
    isTRUE(analysis()$bayesian_artifact_allowed) &&
      validated_selection_selected(confirmed_selected_pc_ids()) &&
      isTRUE(pc_signs_are_default())
  })

  frequentist_current_mediators <- reactive({
    if (isTRUE(frequentist_example_default_available())) {
      return(validated_pc_ids)
    }
    confirmed_selected_pc_ids()
  })

  frequentist_current_cache_key <- reactive({
    req(pca_config_confirmed())
    if (isTRUE(frequentist_example_default_available())) {
      return("emory_frequentist_PC1_PC2_PC3_reversed12_structural_rho")
    }
    frequentist_dynamic_cache_key(
      analysis(),
      selected_pcs = confirmed_selected_pc_ids(),
      pc_signs = confirmed_pc_sign_vector()
    )
  })

  frequentist_current_result <- reactive({
    if (!isTRUE(pca_config_confirmed())) {
      return(NULL)
    }
    if (isTRUE(frequentist_example_default_available())) {
      return(list(
        source = "validated_example",
        key = frequentist_current_cache_key(),
        result = make_frequentist_example_fit(analysis(), confirmed_pc_sign_vector())
      ))
    }

    key <- frequentist_current_cache_key()
    state <- frequentist_fit_state()
    if (!is.null(state) && identical(state$key, key)) {
      return(state)
    }

    cached <- frequentist_fit_cache()[[key]]
    if (!is.null(cached)) {
      frequentist_fit_state(cached)
      return(cached)
    }

    NULL
  })

	  frequentist_unavailable_message <- function() {
	    if (!isTRUE(pca_config_confirmed())) {
	      return("Confirm selected PCs in PCA before running Frequentist Mediation.")
	    }
	    "Click Run Frequentist Mediation to fit the frequentist model for the confirmed prepared data and selected PCs."
	  }

	  mediation_result_available <- reactive({
	    framework <- input$analysis_framework %||% "bayesian"
	    if (identical(framework, "bayesian")) {
	      return(!is.null(bayesian_current_artifact()))
	    }
	    !is.null(frequentist_current_result())
	  })

	  output$mediation_result_filters <- renderUI({
	    if (!isTRUE(mediation_result_available())) {
	      return(NULL)
	    }

	    a <- analysis()
	    contrast_choices <- stats::setNames(unname(a$contrast_labels), unname(a$contrast_labels))
	    selected_contrast <- input$mediation_contrast_filter %||% unname(contrast_choices)[[1]]
	    if (!selected_contrast %in% unname(contrast_choices)) {
	      selected_contrast <- unname(contrast_choices)[[1]]
	    }

	    mediator_ids <- if (identical(input$analysis_framework %||% "bayesian", "frequentist")) {
	      frequentist_current_mediators()
	    } else {
	      bayesian_current_mediators()
	    }
	    mediator_choices <- mediator_filter_choices(mediator_ids)
	    selected_mediator <- input$mediation_mediator_filter %||% "all"
	    if (!selected_mediator %in% unname(mediator_choices)) {
	      selected_mediator <- "all"
	    }

	    tagList(
	      selectInput(
	        "mediation_contrast_filter",
	        app_label_help(
	          "Contrast",
	          "A contrast compares one exposure group with the baseline/reference group defined for this analysis."
	        ),
	        choices = contrast_choices,
	        selected = selected_contrast
	      ),
	      selectInput(
	        "mediation_mediator_filter",
	        app_label_help(
	          "Mediator",
	          "Choose All to view the combined mediator pattern, or select one mediator PC when you want to inspect a specific path in detail."
	        ),
	        choices = mediator_choices,
	        selected = selected_mediator
	      )
	    )
	  })

  output$mediation_view_controls <- renderUI({
    if (!isTRUE(mediation_result_available())) {
      return(NULL)
    }
    app_setup_section("Mediation View", uiOutput("mediation_result_filters"))
  })

	  sensitivity_mediator_ids <- reactive({
    if (identical(input$analysis_framework %||% "bayesian", "bayesian")) {
      artifact <- bayesian_current_artifact()
      if (!is.null(artifact)) {
        return(bayes_A_bundle_mediators(artifact$sensitivity_A$residual_draws))
      }
      return(bayesian_current_mediators())
    }
    frequentist_current_mediators()
  })

  output$sensitivity_mediator_control <- renderUI({
    if (!identical(input$sensitivity_scenario %||% "common", "one_at_a_time")) {
      return(NULL)
    }
    choices <- mediator_filter_choices(sensitivity_mediator_ids())[-1]
    if (length(choices) == 0) {
      return(tags$p(class = "setup-meta", "Confirm mediators and fit mediation first."))
    }
    selected <- input$sensitivity_one_mediator %||% unname(choices)[[1]]
    if (!selected %in% unname(choices)) {
      selected <- unname(choices)[[1]]
    }
    selectInput("sensitivity_one_mediator", "Mediator", choices = choices, selected = selected)
  })

  output$sensitivity_rho_control <- renderUI({
    current_rho <- isolate(input$sensitivity_active_rho)
    sliderInput(
      "sensitivity_active_rho",
      "Residual correlation (ρ)",
      min = min(bayesian_sensitivity_rho_grid),
      max = max(bayesian_sensitivity_rho_grid),
      value = current_rho %||% 0,
      step = 0.05
    )
  })

  bayesian_unavailable_message <- function() {
    if (!isTRUE(pca_config_confirmed())) {
      return("Confirm selected PCs in PCA before running Bayesian Mediation.")
    }
    "Click Run Bayesian Mediation to fit the Bayesian model for the confirmed prepared data and selected PCs."
  }

  output$bayesian_mediation_guard <- renderUI({
    if (isTRUE(bayesian_artifact_available())) {
      return(tags$div(style = "color:#1b6e3c;", "Bayesian results are available for the current selection."))
    }
    tagList(
      fluidRow(
        column(
          width = 4,
          actionButton(
            "run_bayesian_mediation",
            "Run Bayesian Mediation",
            class = "btn-primary"
          )
        ),
        column(width = 8, uiOutput("bayesian_fit_status"))
      ),
      tags$div(style = "color:#8a5a00;", bayesian_unavailable_message())
    )
  })

  output$bayesian_sensitivity_guard <- renderUI({
    if (!is.null(bayesian_current_artifact())) {
      return(NULL)
    }
    tags$div(style = "color:#9b1c1c;", tags$strong("Bayesian sensitivity unavailable: "), "Run Bayesian Mediation first to enable sensitivity analysis.")
  })

  output$frequentist_sensitivity_guard <- renderUI({
    if (!is.null(frequentist_current_result())) {
      return(NULL)
    }
    tags$div(style = "color:#9b1c1c;", tags$strong("Frequentist sensitivity unavailable: "), "Run Frequentist Mediation first for the current data and selected PCs.")
  })

  output$mediation_input_context <- renderUI({
    context <- confirmed_pca_context()
    if (is.null(context)) {
      return(tags$div(
        class = "workflow-context-card",
        h5("Analysis Input"),
        tags$p(class = "workflow-note", "Confirm mediators in PCA before running mediation.")
      ))
    }

    tagList(
      tags$div(
        class = "workflow-context-card",
        h5("Analysis Input"),
        app_kv_item("Exposure", context$exposure),
        app_kv_item("Reference", context$reference),
        app_kv_item("Mediators", context$mediators),
        app_kv_item("Covariates", context$covariates)
      ),
      if (isTRUE(pca_has_unconfirmed_changes())) {
        tags$p(class = "workflow-note warning", "PCA edits are unconfirmed. Mediation still uses the confirmed mediator configuration.")
      }
    )
  })

  output$sensitivity_baseline_context <- renderUI({
    framework <- input$analysis_framework %||% "bayesian"
    context <- confirmed_pca_context()
    contrast <- input$sensitivity_contrast_filter %||% "all"
    contrast_label_ui <- if (identical(contrast, "all")) "All contrasts" else contrast

    if (is.null(context)) {
      return(tags$div(
        class = "workflow-context-card",
        h5("Analysis Input"),
        tags$p(class = "workflow-note", "Confirm mediators in PCA before running sensitivity analysis.")
      ))
    }

    tags$div(
      class = "workflow-context-card",
      h5("Analysis Input"),
      app_kv_item("Framework", if (identical(framework, "bayesian")) "Bayesian" else "Frequentist"),
      app_kv_item("Contrast", contrast_label_ui),
      app_kv_item("Mediators", context$mediators)
    )
  })

  observeEvent(input$go_to_pca_from_overview, {
    updateTabsetPanel(session, "main_nav", selected = "PCA")
  }, ignoreInit = TRUE)

  observe({
	    data_done <- isTRUE(setup_status()$ok) && !is.null(analysis())
	    pca_done <- isTRUE(pca_config_confirmed()) && !isTRUE(pca_has_unconfirmed_changes())
	    mediation_done <- pca_done && (
	      !is.null(bayesian_current_artifact()) ||
	        !is.null(frequentist_current_result())
	    )
    steps <- list(
      list(value = "Data Setup", workflow = TRUE, completed = data_done, available = TRUE),
      list(value = "Overview", workflow = TRUE, completed = data_done, available = data_done),
      list(value = "PCA", workflow = TRUE, completed = pca_done, available = data_done),
      list(value = "Mediation", workflow = TRUE, completed = mediation_done, available = pca_done)
    )
    session$sendCustomMessage("setNavWorkflow", list(steps = steps))
  })

  output$bayesian_fit_status <- renderUI({
    if (isTRUE(bayesian_fit_running())) {
      return(tags$div(
        style = "color:#8a5a00;",
        "Bayesian fitting is running across 3 Stan models. The progress bar shows approximate overall progress; detailed iteration logs are written to the server console."
      ))
    }
    if (!is.null(bayesian_current_artifact())) {
      result <- bayesian_current_result()
      diagnostics <- bayesian_app_mcmc_diagnostics(result$artifact)
      status <- if (!is.null(diagnostics) && all(diagnostics$overall_pass, na.rm = TRUE)) {
        "Bayesian fit is complete and diagnostics passed."
      } else if (!is.null(diagnostics)) {
        "Bayesian fit is complete with diagnostic warnings."
      } else {
        "Bayesian results are ready."
      }
      return(tags$div(style = "color:#1b6e3c;", status))
    }
    tags$div(style = "color:#666;", "No Bayesian fit has been run for the current setup.")
  })

  output$bayesian_diagnostics_notice <- renderUI({
    artifact <- bayesian_current_artifact()
    diagnostics <- bayesian_app_mcmc_diagnostics(artifact)
    if (is.null(artifact) || is.null(diagnostics)) {
      return(NULL)
    }
    color <- if (all(diagnostics$overall_pass, na.rm = TRUE)) "#1b6e3c" else "#8a5a00"
    tags$details(
      tags$summary("Bayesian fit diagnostics"),
      tags$div(
        style = paste0("color:", color, "; margin: 6px 0;"),
        if (all(diagnostics$overall_pass, na.rm = TRUE)) {
          "All monitored MCMC diagnostics passed."
        } else {
          "One or more monitored MCMC diagnostics did not pass. Interpret posterior summaries with caution."
        }
      ),
      tableOutput("bayesian_diagnostics_table")
    )
  })

  output$bayesian_diagnostics_table <- renderTable({
    req(bayesian_current_artifact())
    diagnostics <- bayesian_app_mcmc_diagnostics(bayesian_current_artifact())
    req(diagnostics)
    diagnostics |>
      dplyr::transmute(
        model,
        max_Rhat = round(max_Rhat, 3),
        min_bulk_ESS = round(min_bulk_ESS, 0),
        min_tail_ESS = round(min_tail_ESS, 0),
        divergent_transitions,
        max_treedepth_hits,
        overall_pass
      )
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  observeEvent(input$run_bayesian_mediation, {
    if (!isTRUE(pca_config_confirmed())) {
      showNotification("Confirm mediator configuration in PCA before running Bayesian mediation.", type = "warning")
      return(NULL)
    }
    req(!isTRUE(bayesian_artifact_available()))
    if (isTRUE(bayesian_fit_running())) {
      showNotification("Bayesian mediation is already running for this session.", type = "warning")
      return(NULL)
    }

    key <- bayesian_current_cache_key()
    cached <- bayesian_fit_cache()[[key]]
    if (!is.null(cached)) {
      bayesian_fit_state(cached)
      showNotification("Reused Bayesian mediation results for the current setup.", type = "message")
      session$sendCustomMessage("bayesianClientProgress", list(action = "complete", message = "Reused cached Bayesian results"))
      warm_bayesian_sensitivity_cache(cached)
      return(NULL)
    }

    bayesian_fit_running(TRUE)
    on.exit({
      bayesian_fit_running(FALSE)
    }, add = TRUE)

    settings <- bayesian_dynamic_settings()
    result <- tryCatch({
      artifact <- make_bayesian_dynamic_artifact(
        analysis_obj = analysis(),
        selected_pcs = confirmed_selected_pc_ids(),
        pc_signs = confirmed_pc_sign_vector(),
        settings = settings,
        progress_callback = function(step, detail = NULL) {
          invisible(NULL)
        }
      )
      list(ok = TRUE, artifact = artifact)
    }, error = function(e) {
      list(ok = FALSE, message = conditionMessage(e))
    })

    if (isTRUE(result$ok)) {
      state <- list(source = "dynamic_fit", key = key, artifact = result$artifact)
      cache <- bayesian_fit_cache()
      cache[[key]] <- state
      bayesian_fit_cache(cache)
      bayesian_fit_state(state)
      showNotification("Bayesian mediation fit completed.", type = "message")
      session$sendCustomMessage("bayesianClientProgress", list(action = "complete", message = "Bayesian fit complete"))
      warm_bayesian_sensitivity_cache(state)
    } else {
      showNotification(paste("Bayesian mediation failed:", result$message), type = "error", duration = 12)
      session$sendCustomMessage("bayesianClientProgress", list(action = "error", message = "Bayesian fit failed"))
    }
  }, ignoreInit = TRUE)

  output$frequentist_mediation_guard <- renderUI({
    if (isTRUE(frequentist_example_default_available())) {
      return(tags$div(style = "color:#1b6e3c;", "Frequentist results are available for the current selection."))
    }
    tagList(
      fluidRow(
        column(
          width = 4,
          actionButton(
            "run_frequentist_mediation",
            "Run Frequentist Mediation",
            class = "btn-primary"
          )
        ),
        column(width = 8, uiOutput("frequentist_fit_status"))
      ),
      tags$div(style = "color:#8a5a00;", frequentist_unavailable_message())
    )
  })

  output$frequentist_fit_status <- renderUI({
    if (isTRUE(frequentist_fit_running())) {
      return(tags$div(style = "color:#8a5a00;", "Frequentist mediation is running."))
    }
    if (!is.null(frequentist_current_result())) {
      return(tags$div(style = "color:#1b6e3c;", "Frequentist mediation fit is ready."))
    }
    tags$div(style = "color:#666;", "No frequentist fit has been run for the current setup.")
  })

	  observeEvent(input$run_frequentist_mediation, {
    if (!isTRUE(pca_config_confirmed())) {
      showNotification("Confirm mediator configuration in PCA before running Frequentist mediation.", type = "warning")
      return(NULL)
    }
    if (isTRUE(frequentist_fit_running())) {
      showNotification("Frequentist mediation is already running for this session.", type = "warning")
      return(NULL)
    }

    key <- frequentist_current_cache_key()
    cached <- frequentist_fit_cache()[[key]]
    if (!is.null(cached)) {
      frequentist_fit_state(cached)
      showNotification("Reused frequentist mediation results for the current setup.", type = "message")
      warm_frequentist_sensitivity_cache(cached)
      ensure_frequentist_bootstrap(cached)
      return(NULL)
    }

    frequentist_fit_running(TRUE)
    on.exit({
      frequentist_fit_running(FALSE)
    }, add = TRUE)

    result <- tryCatch({
      fit <- make_frequentist_dynamic_fit(
        analysis_obj = analysis(),
        selected_pcs = confirmed_selected_pc_ids(),
        pc_signs = confirmed_pc_sign_vector()
      )
      list(ok = TRUE, fit = fit)
    }, error = function(e) {
      list(ok = FALSE, message = conditionMessage(e))
    })

    if (isTRUE(result$ok)) {
      state <- list(source = "dynamic_fit", key = key, result = result$fit)
      cache <- frequentist_fit_cache()
      cache[[key]] <- state
      frequentist_fit_cache(cache)
      frequentist_fit_state(state)
      showNotification("Frequentist mediation fit completed.", type = "message")
      warm_frequentist_sensitivity_cache(state)
      ensure_frequentist_bootstrap(state)
    } else {
      showNotification(paste("Frequentist mediation failed:", result$message), type = "error", duration = 12)
    }
  }, ignoreInit = TRUE)

  overview_spatial_df <- reactive({
    analysis()$spatial_inputs$plot_df
  })

  overview_spatial_limits <- reactive({
    spatial_plot_limits(overview_spatial_df())
  })

  overview_exposure_plot <- reactive({
    spec <- analysis()$spec
    limits <- overview_spatial_limits()
    plot <- overview_spatial_df() |>
      ggplot(aes(x = x_plot, y = y_plot, color = X_cat)) +
      geom_point(size = 1.35, alpha = 0.85) +
      coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE) +
      labs(
        title = paste("Exposure Spatial Distribution:", spec$original_labels$exposure),
        x = spec$original_labels$spatial_x,
        y = spec$original_labels$spatial_y,
        color = spec$original_labels$exposure
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        plot.margin = margin(4, 4, 4, 4),
        legend.margin = margin(0, 0, 0, 0)
      )
    if (isTRUE(analysis()$spatial_inputs$flip_y)) {
      plot <- plot + scale_y_continuous(labels = function(z) -z)
    }
    plot
  })

  output$overview_exposure_spatial <- renderPlot(overview_exposure_plot())

  overview_outcome_plot <- reactive({
    spec <- analysis()$spec
    limits <- overview_spatial_limits()
    plot <- overview_spatial_df() |>
      ggplot(aes(x = x_plot, y = y_plot, color = Y)) +
      geom_point(size = 1.35, alpha = 0.9) +
      coord_equal(xlim = limits$x, ylim = limits$y, expand = FALSE) +
      scale_color_gradient(low = "#2C7BB6", high = "#D7191C") +
      labs(
        title = paste("Outcome Spatial Distribution:", spec$original_labels$outcome),
        x = spec$original_labels$spatial_x,
        y = spec$original_labels$spatial_y,
        color = spec$original_labels$outcome
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        plot.margin = margin(4, 4, 4, 4),
        legend.margin = margin(0, 0, 0, 0)
      )
    if (isTRUE(analysis()$spatial_inputs$flip_y)) {
      plot <- plot + scale_y_continuous(labels = function(z) -z)
    }
    plot
  })

  output$overview_outcome_spatial <- renderPlot(overview_outcome_plot())

  output$overview_dataset <- renderUI({
    a <- analysis()
    tags$div(
      class = "overview-metric-row",
      tags$div(
        class = "overview-metric",
        tags$span(class = "overview-metric-value", format(nrow(a$dat_full_allpc), big.mark = ",")),
        tags$span(class = "overview-metric-label", "observations")
      ),
      tags$div(
        class = "overview-metric",
        tags$span(class = "overview-metric-value", format(a$spec$mediator_features$n_features, big.mark = ",")),
        tags$span(class = "overview-metric-label", "mediator features")
      )
    )
  })

  output$overview_exposure <- renderUI({
    counts <- analysis()$dat_full_allpc |>
      count(X_cat, name = "spots") |>
      rename(category = X_cat)

    tags$div(
      lapply(seq_len(nrow(counts)), function(i) {
        app_kv_item(
          as.character(counts$category[[i]]),
          format(counts$spots[[i]], big.mark = ",")
        )
      })
    )
  })

  output$overview_variables <- renderUI({
    a <- analysis()
    tags$div(
      app_kv_item("Exposure", a$spec$original_labels$exposure),
      app_kv_item("Reference", a$spec$exposure$reference),
      app_kv_item("Outcome", a$spec$original_labels$outcome),
      app_kv_item("Coordinates", paste(c(a$spec$original_labels$spatial_x, a$spec$original_labels$spatial_y), collapse = ", ")),
      app_kv_item(
        "Additional Covariates",
        if (length(a$spec$covariates) == 0) "none" else paste(unname(a$spec$covariates), collapse = ", ")
      )
    )
  })

  explained_variance_plot <- reactive({
    analysis()$scree_df |>
      slice_head(n = 30) |>
      ggplot(aes(x = PC, y = var_explained)) +
      geom_col(fill = "#4C78A8") +
      geom_line(aes(y = cum_var_explained), color = "#F58518", linewidth = 0.8) +
      geom_point(aes(y = cum_var_explained), color = "#F58518", size = 1.6) +
      labs(x = "Principal component", y = "Proportion variance explained") +
      theme_bw(base_size = 12)
  })

  output$scree_plot <- renderPlot(explained_variance_plot())

  observeEvent(input$confirm_pca_mediators, {
    spec <- current_pca_spec()
    if (length(spec$selected_pcs) == 0) {
      showNotification("Select at least one PC before confirming mediators.", type = "warning")
      return(NULL)
    }

    previous <- confirmed_pca_spec()
    changed <- !pca_specs_equal(previous, spec)
    confirmed_pca_spec(spec)
    if (isTRUE(changed)) {
	      bayesian_fit_state(NULL)
	      frequentist_fit_state(NULL)
	      frequentist_bootstrap_cache(list())
	      sensitivity_warmup_status(NULL)
    }
    bayes_ready <- bayesian_current_result()
    if (!is.null(bayes_ready)) {
      warm_bayesian_sensitivity_cache(bayes_ready)
    }
    freq_ready <- frequentist_current_result()
    if (!is.null(freq_ready)) {
      warm_frequentist_sensitivity_cache(freq_ready)
    }
    updateTabsetPanel(session, "main_nav", selected = "Mediation")
  }, ignoreInit = TRUE)

  output$pca_confirmation_controls <- renderUI({
    tags$div(
      class = "workflow-action",
      actionButton("confirm_pca_mediators", "Confirm PCs and Continue", class = "btn-primary")
    )
  })

  output$pca_selection_notice <- renderUI({
    selected <- selected_pc_ids(input, analysis())
    if (length(selected) == 0) {
      return(tags$p(class = "setup-meta", tags$span(class = "setup-warning", "!"), "Select at least one PC to display maps and loading genes."))
    }
    NULL
  })

  output$pca_pc_controls <- renderUI({
    candidates <- pc_candidate_ids(analysis(), max_pcs = 6)

    tagList(
      lapply(seq_along(candidates), function(i) {
        pc <- candidates[[i]]
        selected_id <- paste0("pca_pc", i, "_selected")
        is_selected <- input[[selected_id]] %||% (i <= 2)

        tags$div(
          class = "pca-pc-row",
          checkboxInput(selected_id, pc, value = isTRUE(is_selected))
        )
      })
    )
  })

  oriented_spatial_df <- reactive({
    orient_spatial_scores(
      analysis()$spatial_inputs$plot_df,
      pc_sign_vector(input, analysis()),
      pca_fit = analysis()$pca$pca_fit
    )
  })

  oriented_loading_df <- reactive({
    orient_loading_data(
      analysis()$loading_df,
      pc_sign_vector(input, analysis()),
      pca_fit = analysis()$pca$pca_fit
    )
  })

  for (i in seq_len(6)) {
    local({
      pc_index <- i
      pc_id <- paste0("PC", pc_index)
      output[[paste0("pc", pc_index, "_spatial")]] <- renderPlot({
        req(pc_index <= length(analysis()$pca$var_explained))
        req(pc_id %in% names(oriented_spatial_df()))
        plot_pc_spatial(
          oriented_spatial_df(),
          pc_id,
          analysis()$pca$var_explained[[pc_index]],
          pc_id,
          analysis()$spatial_inputs$category_polygon_sf,
          flip_y = analysis()$spatial_inputs$flip_y
        ) +
          labs(title = NULL) +
          theme(
            plot.title = element_blank(),
            plot.margin = margin(4, 4, 4, 4)
          )
      })
    })
  }

  output$pc_spatial_outputs <- renderUI({
    pc_ids <- selected_pc_ids(input, analysis())
    if (length(pc_ids) == 0) {
      return(app_panel_card(
        "No PCs Selected",
        class = "setup-empty",
        tags$p("Select one or more PCs in the sidebar to display spatial score maps.")
      ))
    }
    n_plots <- length(pc_ids)
    plot_ids <- paste0(tolower(pc_ids), "_spatial")

    tags$div(
      class = paste("pca-map-grid", paste0("pc-count-", n_plots)),
      lapply(seq_len(n_plots), function(i) {
        pc_index <- as.integer(sub("^PC", "", pc_ids[[i]]))
        app_panel_card(
          pc_ids[[i]],
          class = "pca-map-card",
          tags$p(
            class = "pca-card-meta",
            paste0(sprintf("%.2f", 100 * analysis()$pca$var_explained[[pc_index]]), "% variance explained")
          ),
          plotOutput(plot_ids[[i]], height = 330)
        )
      })
    )
  })

  output$pca_detail_meta <- renderUI({
    pc <- input$loading_pc %||% "PC1"
    pc_index <- as.integer(sub("^PC", "", pc))
    req(!is.na(pc_index), pc_index <= length(analysis()$pca$var_explained))
    tags$p(
      class = "pca-card-meta",
      paste0(sprintf("%.2f", 100 * analysis()$pca$var_explained[[pc_index]]), "% variance explained")
    )
  })

  output$pca_detail_spatial <- renderPlot({
    pc <- input$loading_pc %||% "PC1"
    pc_index <- as.integer(sub("^PC", "", pc))
    req(!is.na(pc_index), pc_index <= length(analysis()$pca$var_explained))
    req(pc %in% names(oriented_spatial_df()))

    plot_pc_spatial(
      oriented_spatial_df(),
      pc,
      analysis()$pca$var_explained[[pc_index]],
      pc,
      analysis()$spatial_inputs$category_polygon_sf,
      flip_y = analysis()$spatial_inputs$flip_y
    ) +
      labs(title = NULL) +
      theme(
        plot.title = element_blank(),
        plot.margin = margin(4, 4, 4, 4)
      )
  })

  output$loading_pc_control <- renderUI({
    candidates <- pc_candidate_ids(analysis(), max_pcs = 6)
    if (length(candidates) == 0) {
      return(tags$p(class = "setup-meta", "No PCs are available to inspect."))
    }
    current <- input$loading_pc %||% candidates[[1]]
    if (!current %in% candidates) {
      current <- candidates[[1]]
    }
    selectInput("loading_pc", "PC to inspect", choices = candidates, selected = current)
  })

  loading_bar_panel <- function(title, loading_df, direction) {
    values <- abs(loading_df$loading)
    max_value <- max(values, na.rm = TRUE)
    if (!is.finite(max_value) || max_value == 0) {
      max_value <- 1
    }

	    tags$section(
	      class = paste("loading-panel", direction),
      h5(title),
      lapply(seq_len(nrow(loading_df)), function(i) {
        width <- 100 * abs(loading_df$loading[[i]]) / max_value
        tags$div(
          class = "loading-bar-row",
          tags$div(
            class = "loading-bar-head",
            tags$span(class = "loading-gene", loading_df$gene[[i]]),
            tags$span(class = "loading-value", formatC(loading_df$loading[[i]], format = "f", digits = 5))
          ),
          tags$div(
            class = "loading-bar-track",
            tags$div(
              class = paste("loading-bar-fill", direction),
              style = paste0("width:", sprintf("%.1f", width), "%;")
            )
          )
        )
      })
    )
  }

  output$loading_results <- renderUI({
    req(input$loading_pc, input$loading_n)
    n <- as.integer(input$loading_n)
    pos <- top_loading_table(oriented_loading_df(), input$loading_pc, "positive", n = n)
    neg <- top_loading_table(oriented_loading_df(), input$loading_pc, "negative", n = n)

    tags$div(
      class = "loading-result-grid",
      loading_bar_panel("Top Positive Loadings", pos, "positive"),
      loading_bar_panel("Top Negative Loadings", neg, "negative")
    )
  })

  mediation_summary <- reactive({
    result <- frequentist_current_result()
    validate(need(!is.null(result), frequentist_unavailable_message()))
    result$result$fit$summary |>
      transmute(
        contrast,
        TE,
        direct_effect = NDE,
        total_indirect_effect = NIE,
        PM
      ) |>
      filter_by_mediation_inputs(
        contrast = input$mediation_contrast_filter %||% "all",
        mediator = "all"
      )
	  })

	  mediation_path <- reactive({
	    result <- frequentist_current_result()
	    validate(need(!is.null(result), frequentist_unavailable_message()))
	    result$result$fit$path |>
	      transmute(
        contrast,
        mediator,
        alpha = alpha_X_to_M,
        beta = beta_M_to_Y,
        indirect_component
	      ) |>
	      mutate(
	        mediator_id = mediator,
	        mediator = mediator_label(mediator)
	      ) |>
	      filter_by_mediation_inputs(
	        contrast = input$mediation_contrast_filter %||% "all",
	        mediator = input$mediation_mediator_filter %||% "all"
	      )
	  })

	  bayesian_mediation_summary <- reactive({
	    artifact <- bayesian_current_artifact()
	    validate(need(!is.null(artifact), bayesian_unavailable_message()))
	    get_bayesian_mediation_summary(artifact) |>
	      filter_by_mediation_inputs(
	        contrast = input$mediation_contrast_filter %||% "all",
	        mediator = "all"
	      )
	  })

	  bayesian_mediation_path <- reactive({
	    artifact <- bayesian_current_artifact()
	    validate(need(!is.null(artifact), bayesian_unavailable_message()))
	    path <- get_bayesian_mediation_path(artifact)
	    # Both cached and runtime artifacts already use the confirmed axes.
	    path <- path |>
	        mutate(
	          mediator_id = mediator,
	          mediator = mediator_label(mediator)
	        )
	    path |>
	      filter_by_mediation_inputs(
	        contrast = input$mediation_contrast_filter %||% "all",
	        mediator = input$mediation_mediator_filter %||% "all"
	      )
	  })

	  bootstrap_cache_key <- function(freq, B) {
	    paste(freq$key, paste0("B=", as.integer(B)), sep = ":")
	  }

	  compute_bootstrap_result <- function(freq, B) {
	    boot <- withProgress(message = "Running bootstrap confidence intervals", value = 0, {
	      bootstrap_frequentist_dynamic(
	        frequentist_result = freq$result,
	        B = B,
	        seed = 123,
	        rho_grid = bayesian_sensitivity_rho_grid,
	        progress_callback = function(done, total) {
	          incProgress(
	            amount = 1 / total,
	            detail = paste0(done, " / ", total, " resamples with sensitivity")
	          )
	        }
	      )
	    })
	    list(
	      key = freq$key,
	      B = B,
	      ci_method = "empirical percentile quantiles at 2.5% and 97.5%",
	      ci = boot$ci,
	      bootstrap = boot$bootstrap,
	      sensitivity_bundle = boot$sensitivity_bundle,
	      sensitivity_ci = summarize_frequentist_sensitivity_bootstrap_ci(
	        boot$sensitivity_bundle,
	        rho_grid = bayesian_sensitivity_rho_grid
	      )
	    )
	  }

	  ensure_frequentist_bootstrap <- function(freq = frequentist_current_result()) {
	    if (is.null(freq)) {
	      return(NULL)
	    }
	    B <- as.integer(input$bootstrap_B %||% 1000)
	    validate(need(B > 0, "Bootstrap replicates must be positive."))
	    key <- bootstrap_cache_key(freq, B)
	    cache <- frequentist_bootstrap_cache()
	    if (is.null(cache[[key]])) {
	      cache[[key]] <- compute_bootstrap_result(freq, B)
	      frequentist_bootstrap_cache(cache)
	    }
	    invisible(cache[[key]])
	  }

	  bootstrap_result <- reactive({
	    freq <- frequentist_current_result()
	    if (is.null(freq)) {
	      return(NULL)
	    }
	    B <- as.integer(input$bootstrap_B %||% 1000)
	    validate(need(B > 0, "Bootstrap replicates must be positive."))
	    key <- bootstrap_cache_key(freq, B)
	    cache <- frequentist_bootstrap_cache()
	    cache[[key]]
	  })

	  observeEvent(input$analysis_framework, {
	    if (identical(input$analysis_framework, "frequentist")) {
	      ensure_frequentist_bootstrap()
	    }
	  }, ignoreInit = TRUE)

	  observeEvent(input$bootstrap_B, {
	    if (identical(input$analysis_framework, "frequentist")) {
	      ensure_frequentist_bootstrap()
	    }
	  }, ignoreInit = TRUE)

	  bootstrap_result_or_null <- reactive({
	    tryCatch(
	      bootstrap_result(),
	      shiny.silent.error = function(e) NULL,
      validation = function(e) NULL,
      error = function(e) NULL
    )
  })

  current_bootstrap_summary_ci <- reactive({
    boot <- bootstrap_result_or_null()
    if (is.null(boot) || !identical(boot$key, frequentist_current_cache_key())) {
      return(NULL)
    }
    boot$ci$summary_ci
  })

  current_bootstrap_path_ci <- reactive({
    boot <- bootstrap_result_or_null()
    if (is.null(boot) || !identical(boot$key, frequentist_current_cache_key())) {
      return(NULL)
    }
    boot$ci$path_ci
  })

  current_bootstrap_sensitivity_ci <- reactive({
    boot <- bootstrap_result_or_null()
    if (is.null(boot) || !identical(boot$key, frequentist_current_cache_key())) {
      return(NULL)
    }
    boot$sensitivity_ci
  })

  mediation_path_with_bootstrap_ci <- reactive({
    rows <- mediation_path()
    ci <- current_bootstrap_path_ci()
    if (is.null(ci)) {
      return(rows)
    }
    rows |>
      left_join(
        ci |>
          select(
            contrast,
            mediator,
            alpha_lwr,
            alpha_upr,
            beta_lwr,
            beta_upr,
            indirect_lwr,
            indirect_upr
          ),
        by = c("contrast", "mediator")
      )
  })

  mediation_interval_text <- function(lwr, upr, percent = FALSE) {
    if (is.null(lwr) || is.null(upr) || is.na(lwr) || is.na(upr)) {
      return(NULL)
    }
    if (isTRUE(percent)) {
      return(paste0("[", format_number(100 * lwr, 1), "%, ", format_number(100 * upr, 1), "%]"))
    }
    paste0("[", format_number(lwr, 3), ", ", format_number(upr, 3), "]")
  }

  mediation_value_text <- function(x, lwr = NA_real_, upr = NA_real_, percent = FALSE) {
    if (isTRUE(percent)) {
      return(paste0(format_number(100 * x, 1), interval_star(lwr, upr), "%"))
    }
    paste0(format_number(x, 3), interval_star(lwr, upr))
  }

  mediation_probability_text <- function(p) {
    if (is.null(p) || is.na(p)) {
      return(NULL)
    }
    paste0("Pr(effect > 0): ", format_number(100 * p, 1), "%")
  }

  mediation_effect_card <- function(label, value, lwr = NA_real_, upr = NA_real_, probability = NULL, percent = FALSE) {
    interval <- mediation_interval_text(lwr, upr, percent = percent)
    tags$div(
      class = "effect-card",
      tags$div(class = "effect-label", label),
      tags$div(class = "effect-value", mediation_value_text(value, lwr, upr, percent = percent)),
      if (!is.null(interval)) tags$div(class = "effect-interval", interval),
      if (!is.null(probability)) tags$div(class = "effect-evidence", mediation_probability_text(probability))
    )
  }

  mediation_summary_panel <- function(contrast, cards) {
    app_panel_card(
      contrast,
      class = "mediation-contrast-panel",
      tags$div(class = "effect-card-grid", cards)
    )
  }

  mediation_dag_number <- function(x, lwr = NA_real_, upr = NA_real_, digits = 3) {
    if (is.null(x) || length(x) == 0 || is.na(x[[1]])) {
      return("--")
    }
    paste0(format_number(x[[1]], digits = digits), interval_star(lwr, upr))
  }

  mediation_dag_interval <- function(lwr, upr) {
    interval <- mediation_interval_text(lwr, upr)
    if (is.null(interval)) {
      return(NULL)
    }
    tags$span(interval)
  }

  mediation_dag_interval_plain <- function(lwr, upr) {
    interval <- mediation_interval_text(lwr, upr)
    if (is.null(interval)) {
      return("interval unavailable")
    }
    interval
  }

  mediation_dag_cell <- function(est, lwr = NA_real_, upr = NA_real_) {
    tagList(
      mediation_dag_number(est),
      if (!is.null(mediation_dag_interval(lwr, upr))) {
        tagList(" ", mediation_dag_interval(lwr, upr))
      }
    )
  }

  mediation_dag_path_values <- function(row, framework) {
    if (identical(framework, "bayesian")) {
      return(list(
        alpha = row$alpha_est,
        alpha_lwr = row$alpha_lwr,
        alpha_upr = row$alpha_upr,
        beta = row$beta_est,
        beta_lwr = row$beta_lwr,
        beta_upr = row$beta_upr,
        indirect = row$IIE_est,
        indirect_lwr = row$IIE_lwr,
        indirect_upr = row$IIE_upr
      ))
    }
    value_or_na <- function(name) {
      if (!name %in% names(row)) {
        return(NA_real_)
      }
      row[[name]]
    }
    list(
      alpha = row$alpha,
      alpha_lwr = value_or_na("alpha_lwr"),
      alpha_upr = value_or_na("alpha_upr"),
      beta = row$beta,
      beta_lwr = value_or_na("beta_lwr"),
      beta_upr = value_or_na("beta_upr"),
      indirect = row$indirect_component,
      indirect_lwr = value_or_na("indirect_lwr"),
      indirect_upr = value_or_na("indirect_upr")
    )
  }

  mediation_dag_tooltip <- function(label, est, lwr = NA_real_, upr = NA_real_) {
    paste0(label, ": ", mediation_dag_number(est, lwr, upr), " ", mediation_dag_interval_plain(lwr, upr))
  }

  mediation_dag_mediator_rows <- function(path_rows, framework) {
    if (nrow(path_rows) == 0) {
      return(tags$div(class = "mediation-node-detail", "Selected PCs"))
    }
    rows <- lapply(seq_len(nrow(path_rows)), function(i) {
      row <- path_rows[i, ]
      values <- mediation_dag_path_values(row, framework)
      tooltip <- paste(
        mediation_dag_tooltip("X -> M", values$alpha, values$alpha_lwr, values$alpha_upr),
        mediation_dag_tooltip("M -> Y", values$beta, values$beta_lwr, values$beta_upr),
        mediation_dag_tooltip("Indirect effect", values$indirect, values$indirect_lwr, values$indirect_upr),
        sep = "\n"
      )
      tags$div(
        class = "mediation-coeff-row",
        title = tooltip,
        tags$span(class = "mediation-coeff-name", row$mediator),
        tags$span(class = "mediation-coeff-value", paste0("a=", mediation_dag_number(values$alpha, values$alpha_lwr, values$alpha_upr))),
        tags$span(class = "mediation-coeff-value", paste0("b=", mediation_dag_number(values$beta, values$beta_lwr, values$beta_upr))),
        tags$span(class = "mediation-coeff-value mediation-coeff-indirect", paste0("a*b=", mediation_dag_number(values$indirect, values$indirect_lwr, values$indirect_upr)))
      )
    })
    tags$div(class = "mediation-coeff-list", rows)
  }

  mediation_path_help_text <- function(path) {
    if (identical(path, "alpha")) {
      return("a is the exposure-to-mediator path. It shows how the selected exposure contrast is associated with each mediator PC.")
    }
    "b is the mediator-to-outcome path. It shows how each mediator PC is associated with the outcome after accounting for exposure and covariates."
  }

  mediation_dag_path_label <- function(symbol, path_rows = NULL, framework = NULL, path = c("alpha", "beta"), help_style = c("icon", "label")) {
    path <- match.arg(path)
    help_style <- match.arg(help_style)
    path_help <- mediation_path_help_text(path)
    symbol_help <- if (identical(help_style, "label")) {
      app_label_help(tags$span(class = "mediation-path-symbol", symbol), path_help, class = "mediation-path-concept")
    } else {
      tags$span(
        class = "mediation-path-symbol-wrap",
        tags$span(class = "mediation-path-symbol", symbol),
        app_help_tooltip(path_help)
      )
    }
    if (is.null(path_rows) || nrow(path_rows) == 0) {
      return(tags$div(
        symbol_help
      ))
    }
    rows <- lapply(seq_len(nrow(path_rows)), function(i) {
      row <- path_rows[i, ]
      values <- mediation_dag_path_values(row, framework)
      if (identical(path, "alpha")) {
        est <- values$alpha
        lwr <- values$alpha_lwr
        upr <- values$alpha_upr
        tooltip <- mediation_dag_tooltip("X -> M", est, lwr, upr)
      } else {
        est <- values$beta
        lwr <- values$beta_lwr
        upr <- values$beta_upr
        tooltip <- mediation_dag_tooltip("M -> Y", est, lwr, upr)
      }
      tags$div(
        class = "mediation-path-line",
        tags$span(row$mediator),
        tags$span(
          class = "mediation-path-est interval-tooltip",
          `data-tooltip` = tooltip,
          tabindex = "0",
          mediation_dag_number(est, lwr, upr)
        )
      )
    })
    tags$div(
      symbol_help,
      tags$div(class = "mediation-path-lines", rows)
    )
  }

  mediation_dag_indirect_rows <- function(path_rows, framework, total_indirect, total_indirect_lwr = NA_real_, total_indirect_upr = NA_real_) {
    rows <- lapply(seq_len(nrow(path_rows)), function(i) {
      row <- path_rows[i, ]
      values <- mediation_dag_path_values(row, framework)
      tags$div(
        class = "mediation-indirect-row",
        tags$span(paste0(row$mediator, " indirect")),
        tags$span(
          class = "mediation-indirect-value interval-tooltip",
          `data-tooltip` = mediation_dag_tooltip("Indirect effect", values$indirect, values$indirect_lwr, values$indirect_upr),
          tabindex = "0",
          mediation_dag_number(values$indirect, values$indirect_lwr, values$indirect_upr)
        )
      )
    })
    tagList(
      tags$div(
        class = "mediation-indirect-row total",
        tags$span("Total indirect"),
        tags$span(
          class = "mediation-indirect-value interval-tooltip",
          `data-tooltip` = mediation_dag_tooltip("Total indirect effect", total_indirect, total_indirect_lwr, total_indirect_upr),
          tabindex = "0",
          mediation_dag_number(total_indirect, total_indirect_lwr, total_indirect_upr)
        )
      ),
      rows
    )
  }

  mediation_dag_effect_card <- function(title, formula, label, est, lwr = NA_real_, upr = NA_real_, extra = NULL, help_text = NULL, help_style = c("icon", "label")) {
    help_style <- match.arg(help_style)
    heading <- tagList(tags$span(title), tags$code(formula))
    heading <- if (!is.null(help_text) && identical(help_style, "label")) {
      app_label_help(heading, help_text, class = "mediation-effect-concept")
    } else {
      tagList(heading, if (!is.null(help_text)) app_help_tooltip(help_text))
    }
    tags$div(
      class = "mediation-effect-card",
      tags$h5(heading),
      tags$div(
        class = "mediation-effect-main",
        tags$span(label),
        tags$span(
          class = "mediation-effect-value interval-tooltip",
          `data-tooltip` = mediation_dag_tooltip(title, est, lwr, upr),
          tabindex = "0",
          mediation_dag_number(est, lwr, upr)
        )
      ),
      extra
    )
  }

  mediation_dag_effect_grid <- function(path_rows, framework,
                                        indirect_est, indirect_lwr, indirect_upr,
                                        direct_est, direct_lwr, direct_upr,
                                        total_est, total_lwr, total_upr,
                                        help_style = c("icon", "label")) {
    help_style <- match.arg(help_style)
    tags$div(
      class = "mediation-effect-grid",
      mediation_dag_effect_card(
        "Indirect effect",
        "a x b",
        "Total indirect",
        indirect_est,
        indirect_lwr,
        indirect_upr,
        help_style = help_style,
        help_text = "The indirect effect is the part of the exposure-outcome association that goes through the mediator path. For each mediator, it is calculated as a x b; total indirect combines the selected mediators.",
        extra = tags$div(
          class = "mediation-indirect-list",
          lapply(seq_len(nrow(path_rows)), function(i) {
            row <- path_rows[i, ]
            values <- mediation_dag_path_values(row, framework)
            tags$div(
              class = "mediation-indirect-row",
              tags$span(row$mediator),
              tags$span(
                class = "mediation-indirect-value interval-tooltip",
                `data-tooltip` = mediation_dag_tooltip("Indirect effect", values$indirect, values$indirect_lwr, values$indirect_upr),
                tabindex = "0",
                mediation_dag_number(values$indirect, values$indirect_lwr, values$indirect_upr)
              )
            )
          })
        )
      ),
      mediation_dag_effect_card(
        "Direct effect",
        "c'",
        "Direct path",
        direct_est,
        direct_lwr,
        direct_upr
      ),
      mediation_dag_effect_card(
        "Total effect",
        "a x b + c'",
        "Total",
        total_est,
        total_lwr,
        total_upr,
        help_style = help_style,
        help_text = "The total effect summarizes the overall exposure-to-outcome association, including both the mediated pathway through the selected mediator PCs and the direct pathway."
      )
    )
  }

  mediation_dag_canvas <- function(mediator_content, x_detail, y_detail, alpha_label, beta_label, direct_label,
                                  alpha_tooltip = NULL, beta_tooltip = NULL, direct_tooltip = NULL) {
    tags$div(
      class = "mediation-dag-canvas",
      tags$svg(
        class = "mediation-arrow-svg",
        viewBox = "0 0 100 100",
        preserveAspectRatio = "none",
        xmlns = "http://www.w3.org/2000/svg",
        tags$defs(
          tags$marker(id = "dag-arrow-alpha", markerWidth = "8", markerHeight = "8", refX = "7", refY = "4", orient = "auto", markerUnits = "strokeWidth",
                      tags$polygon(points = "0 0, 8 4, 0 8", fill = "#2f5597")),
          tags$marker(id = "dag-arrow-beta", markerWidth = "8", markerHeight = "8", refX = "7", refY = "4", orient = "auto", markerUnits = "strokeWidth",
                      tags$polygon(points = "0 0, 8 4, 0 8", fill = "#7e57c2")),
          tags$marker(id = "dag-arrow-direct", markerWidth = "8", markerHeight = "8", refX = "7", refY = "4", orient = "auto", markerUnits = "strokeWidth",
                      tags$polygon(points = "0 0, 8 4, 0 8", fill = "#6b7b8f"))
        ),
        tags$line(class = "mediation-dag-line mediation-dag-line-alpha", x1 = "25", y1 = "78", x2 = "44", y2 = "30", `marker-end` = "url(#dag-arrow-alpha)"),
        tags$line(class = "mediation-dag-line mediation-dag-line-beta", x1 = "56", y1 = "30", x2 = "75", y2 = "78", `marker-end` = "url(#dag-arrow-beta)"),
        tags$line(class = "mediation-dag-line mediation-dag-line-direct", x1 = "25", y1 = "82", x2 = "75", y2 = "82", `marker-end` = "url(#dag-arrow-direct)")
      ),
      tags$div(
        class = "mediation-arrow-label mediation-arrow-label-xm",
        title = alpha_tooltip %||% "",
        alpha_label
      ),
      tags$div(
        class = "mediation-arrow-label mediation-arrow-label-my",
        title = beta_tooltip %||% "",
        beta_label
      ),
      tags$div(
        class = "mediation-arrow-label mediation-arrow-label-xy",
        title = direct_tooltip %||% "",
        direct_label
      ),
      tags$div(
        class = "mediation-node mediation-node-x",
        tags$div(class = "mediation-node-label", "Exposure"),
        tags$div(class = "mediation-node-value", "X"),
        tags$div(class = "mediation-node-detail", x_detail)
      ),
      tags$div(
        class = "mediation-node mediation-node-m",
        tags$div(class = "mediation-node-label", "Mediator"),
        tags$div(class = "mediation-node-value", "M"),
        mediator_content
      ),
      tags$div(
        class = "mediation-node mediation-node-y",
        tags$div(class = "mediation-node-label", "Outcome"),
        tags$div(class = "mediation-node-value", "Y"),
        tags$div(class = "mediation-node-detail", y_detail)
      )
    )
  }

  mediation_dag_placeholder <- function(framework) {
    app_panel_card(
      "Mediation Model",
      class = "mediation-dag-card",
      tags$div(
        class = "mediation-dag-layout",
        mediation_dag_canvas(
          mediator_content = tags$div(class = "mediation-node-detail", "Selected PCs"),
          x_detail = "Exposure contrast",
          y_detail = "Outcome",
          alpha_label = mediation_dag_path_label("a", path = "alpha", help_style = "label"),
          beta_label = mediation_dag_path_label("b", path = "beta", help_style = "label"),
          direct_label = tags$span(
            class = "mediation-path-symbol-wrap",
            app_label_help(
              tags$span(class = "mediation-path-symbol", "c'"),
              "c' is the direct exposure-to-outcome path. It represents the exposure effect that is not explained by the selected mediator PCs.",
              class = "mediation-path-concept"
            )
          )
        ),
        tags$div(
          class = "mediation-dag-summary-strip",
          tags$div(
            class = "mediation-dag-formula",
            tags$strong("Mediation model"),
            tags$p(class = "setup-meta", "This diagram separates the exposure effect on the outcome into a direct path and indirect paths through the selected mediator PCs.")
          )
        )
      )
    )
  }

  mediation_dag_card <- function(summary_rows, path_rows, framework, interval_label) {
    if (nrow(summary_rows) == 0 || nrow(path_rows) == 0) {
      return(mediation_dag_placeholder(framework))
    }
    contrast <- unique(path_rows$contrast)[[1]]
    path_rows <- dplyr::filter(path_rows, contrast == !!contrast)
    summary_row <- dplyr::filter(summary_rows, contrast == !!contrast) |> dplyr::slice_head(n = 1)
    if (nrow(summary_row) == 0) {
      summary_row <- summary_rows |> dplyr::slice_head(n = 1)
    }

    direct_est <- if (identical(framework, "bayesian")) summary_row$direct_est else summary_row$direct_effect
    indirect_est <- if (identical(framework, "bayesian")) summary_row$total_IIE_est else summary_row$total_indirect_effect
    indirect_lwr <- if (identical(framework, "bayesian")) {
      summary_row$total_IIE_lwr
    } else if ("NIE_lwr" %in% names(summary_row)) {
      summary_row$NIE_lwr
    } else {
      NA_real_
    }
    indirect_upr <- if (identical(framework, "bayesian")) {
      summary_row$total_IIE_upr
    } else if ("NIE_upr" %in% names(summary_row)) {
      summary_row$NIE_upr
    } else {
      NA_real_
    }
    total_est <- if (identical(framework, "bayesian")) summary_row$TE_est else summary_row$TE
    total_lwr <- if (identical(framework, "bayesian")) {
      summary_row$TE_lwr
    } else if ("TE_lwr" %in% names(summary_row)) {
      summary_row$TE_lwr
    } else {
      NA_real_
    }
    total_upr <- if (identical(framework, "bayesian")) {
      summary_row$TE_upr
    } else if ("TE_upr" %in% names(summary_row)) {
      summary_row$TE_upr
    } else {
      NA_real_
    }
    pm_est <- if (identical(framework, "bayesian")) summary_row$PM_est else summary_row$PM
    direct_lwr <- if (identical(framework, "bayesian")) {
      summary_row$direct_lwr
    } else if ("NDE_lwr" %in% names(summary_row)) {
      summary_row$NDE_lwr
    } else {
      NA_real_
    }
    direct_upr <- if (identical(framework, "bayesian")) {
      summary_row$direct_upr
    } else if ("NDE_upr" %in% names(summary_row)) {
      summary_row$NDE_upr
    } else {
      NA_real_
    }

    alpha_label <- mediation_dag_path_label("a", path_rows, framework, path = "alpha", help_style = "label")
    beta_label <- mediation_dag_path_label("b", path_rows, framework, path = "beta", help_style = "label")
    direct_tooltip <- mediation_dag_tooltip("Direct effect", direct_est, direct_lwr, direct_upr)
    direct_label <- tags$span(
      class = "mediation-path-symbol-wrap",
      app_label_help(
        tags$span(class = "mediation-path-symbol", "c'"),
        "c' is the direct exposure-to-outcome path. It represents the exposure effect that is not explained by the selected mediator PCs.",
        class = "mediation-path-concept"
      ),
      tags$span("="),
      tags$span(
        class = "interval-tooltip",
        `data-tooltip` = direct_tooltip,
        tabindex = "0",
        mediation_dag_number(direct_est, direct_lwr, direct_upr)
      )
    )
    alpha_tooltip <- paste(vapply(seq_len(nrow(path_rows)), function(i) {
      row <- path_rows[i, ]
      values <- mediation_dag_path_values(row, framework)
      paste0(row$mediator, " ", mediation_dag_tooltip("X -> M", values$alpha, values$alpha_lwr, values$alpha_upr))
    }, character(1)), collapse = "\n")
    beta_tooltip <- paste(vapply(seq_len(nrow(path_rows)), function(i) {
      row <- path_rows[i, ]
      values <- mediation_dag_path_values(row, framework)
      paste0(row$mediator, " ", mediation_dag_tooltip("M -> Y", values$beta, values$beta_lwr, values$beta_upr))
    }, character(1)), collapse = "\n")

    app_panel_card(
      "Mediation Model",
      class = "mediation-dag-card",
      tags$div(
        class = "mediation-dag-layout",
        mediation_dag_canvas(
          mediator_content = tags$div(class = "mediation-node-detail", paste(path_rows$mediator, collapse = ", ")),
          x_detail = contrast,
          y_detail = "Outcome",
          alpha_label = alpha_label,
          beta_label = beta_label,
          direct_label = direct_label,
          alpha_tooltip = alpha_tooltip,
          beta_tooltip = beta_tooltip,
          direct_tooltip = direct_tooltip
        ),
        tags$div(
          class = "mediation-dag-summary-strip",
          tags$div(
            class = "mediation-dag-formula",
            tags$div(
              class = "mediation-dag-formula-head",
              tags$strong("Effect decomposition")
            ),
            mediation_dag_effect_grid(
              path_rows,
              framework,
              indirect_est,
              indirect_lwr,
              indirect_upr,
              direct_est,
              direct_lwr,
              direct_upr,
              total_est,
              total_lwr,
              total_upr,
              help_style = "label"
            ),
            tags$div(
              class = "mediation-interval-guide",
              tags$span("95% intervals are available on hover. * marks intervals that do not include 0.")
            )
          )
        )
      )
    )
  }

  output$mediation_model_diagram <- renderUI({
    framework <- input$analysis_framework %||% "bayesian"
    if (identical(framework, "bayesian")) {
      return(mediation_dag_card(
        summary_rows = bayesian_mediation_summary(),
        path_rows = bayesian_mediation_path(),
        framework = framework,
        interval_label = "Values are posterior means with 95% credible intervals where shown."
      ))
    }
    summary_rows <- mediation_summary()
    summary_ci <- current_bootstrap_summary_ci()
    if (!is.null(summary_ci)) {
      summary_rows <- summary_rows |>
        left_join(
          summary_ci |> select(contrast, TE_lwr, TE_upr, NDE_lwr, NDE_upr, NIE_lwr, NIE_upr),
          by = "contrast"
        )
    }
    mediation_dag_card(
      summary_rows = summary_rows,
      path_rows = mediation_path_with_bootstrap_ci(),
      framework = framework,
      interval_label = if (is.null(current_bootstrap_path_ci())) {
        "Values are point estimates while bootstrap confidence intervals are prepared."
      } else {
        "Values include 95% bootstrap confidence intervals where shown."
      }
    )
  })

  output$mediation_summary_cards <- renderUI({
    framework <- input$analysis_framework %||% "bayesian"
    if (identical(framework, "bayesian")) {
      rows <- bayesian_mediation_summary()
      return(tagList(lapply(seq_len(nrow(rows)), function(i) {
        row <- rows[i, ]
        mediation_summary_panel(
          row$contrast,
          list(
            mediation_effect_card("Total Effect", row$TE_est, row$TE_lwr, row$TE_upr, row$TE_Pr_gt_0),
            mediation_effect_card("Direct Effect", row$direct_est, row$direct_lwr, row$direct_upr, row$direct_Pr_gt_0),
            mediation_effect_card("Indirect Effect", row$total_IIE_est, row$total_IIE_lwr, row$total_IIE_upr, row$total_IIE_Pr_gt_0),
            mediation_effect_card("Proportion Mediated", row$PM_est, row$PM_lwr, row$PM_upr, percent = TRUE)
          )
        )
      })))
    }

    rows <- mediation_summary()
    ci <- current_bootstrap_summary_ci()
    if (!is.null(ci)) {
      rows <- rows |>
        left_join(
          ci |> select(contrast, TE_lwr, TE_upr, NDE_lwr, NDE_upr, NIE_lwr, NIE_upr, PM_lwr, PM_upr),
          by = "contrast"
        )
    } else {
      rows <- rows |>
        mutate(
          TE_lwr = NA_real_, TE_upr = NA_real_,
          NDE_lwr = NA_real_, NDE_upr = NA_real_,
          NIE_lwr = NA_real_, NIE_upr = NA_real_,
          PM_lwr = NA_real_, PM_upr = NA_real_
        )
    }

    tagList(lapply(seq_len(nrow(rows)), function(i) {
      row <- rows[i, ]
      mediation_summary_panel(
        row$contrast,
        list(
          mediation_effect_card("Total Effect", row$TE, row$TE_lwr, row$TE_upr),
          mediation_effect_card("Direct Effect", row$direct_effect, row$NDE_lwr, row$NDE_upr),
          mediation_effect_card("Indirect Effect", row$total_indirect_effect, row$NIE_lwr, row$NIE_upr),
          mediation_effect_card("Proportion Mediated", row$PM, row$PM_lwr, row$PM_upr, percent = TRUE)
        )
      )
    }))
  })

  bayesian_path_evidence <- function(pr_gt, pr_lt) {
    if (is.na(pr_gt) && is.na(pr_lt)) {
      return("")
    }
    if (!is.na(pr_gt) && (is.na(pr_lt) || pr_gt >= pr_lt)) {
      return(paste0("Pr(positive) = ", format_number(100 * pr_gt, 1), "%"))
    }
    paste0("Pr(negative) = ", format_number(100 * pr_lt, 1), "%")
  }

  frequentist_path_cell <- function(est, lwr = NA_real_, upr = NA_real_, digits = 4) {
    if (is.na(lwr) || is.na(upr)) {
      return(format_number(est, digits = digits))
    }
    format_estimate_interval(est, lwr, upr, digits = digits)
  }

  mediation_path_table_ui <- function(df, framework) {
    if (nrow(df) == 0) {
      return(tags$p(class = "setup-meta", "No path rows match the current filters."))
    }
	    header <- if (identical(framework, "bayesian")) {
	      c("Contrast", "Mediator", "Exposure -> Mediator", "Mediator -> Outcome", "Indirect Effect", "Posterior Direction")
	    } else {
	      c("Contrast", "Mediator", "Exposure -> Mediator", "Mediator -> Outcome", "Indirect Effect")
	    }
    rows <- lapply(seq_len(nrow(df)), function(i) {
      if (identical(framework, "bayesian")) {
        tags$tr(
          tags$td(df$contrast[[i]]),
          tags$td(df$mediator[[i]]),
          tags$td(format_estimate_interval(df$alpha_est[[i]], df$alpha_lwr[[i]], df$alpha_upr[[i]])),
          tags$td(format_estimate_interval(df$beta_est[[i]], df$beta_lwr[[i]], df$beta_upr[[i]])),
          tags$td(format_estimate_interval(df$IIE_est[[i]], df$IIE_lwr[[i]], df$IIE_upr[[i]])),
          tags$td(bayesian_path_evidence(df$IIE_Pr_gt_0[[i]], df$IIE_Pr_lt_0[[i]]))
        )
      } else {
        value_or_na <- function(column) {
          if (!column %in% names(df)) {
            return(NA_real_)
          }
          df[[column]][[i]]
        }
        tags$tr(
          tags$td(df$contrast[[i]]),
          tags$td(df$mediator[[i]]),
          tags$td(frequentist_path_cell(df$alpha[[i]], value_or_na("alpha_lwr"), value_or_na("alpha_upr"), digits = 4)),
          tags$td(frequentist_path_cell(df$beta[[i]], value_or_na("beta_lwr"), value_or_na("beta_upr"), digits = 4)),
          tags$td(frequentist_path_cell(df$indirect_component[[i]], value_or_na("indirect_lwr"), value_or_na("indirect_upr"), digits = 4))
        )
      }
    })
    tags$table(
      class = "mediation-path-table",
      tags$thead(tags$tr(lapply(header, tags$th))),
      tags$tbody(rows)
    )
  }

  output$mediation_path_details <- renderUI({
    framework <- input$analysis_framework %||% "bayesian"
    if (identical(framework, "bayesian")) {
      return(mediation_path_table_ui(bayesian_mediation_path(), "bayesian"))
    }
    mediation_path_table_ui(mediation_path_with_bootstrap_ci(), "frequentist")
  })

  output$mediation_run_controls <- renderUI({
    if (isTRUE(mediation_result_available())) {
      return(NULL)
    }
    framework <- input$analysis_framework %||% "bayesian"
    button_id <- if (identical(framework, "bayesian")) "run_bayesian_mediation" else "run_frequentist_mediation"
    button_label <- if (identical(framework, "bayesian")) "Run Bayesian Mediation" else "Run Frequentist Mediation"
    running <- if (identical(framework, "bayesian")) isTRUE(bayesian_fit_running()) else isTRUE(frequentist_fit_running())
    tags$div(class = "mediation-run-action",
      tags$div(class = "workflow-action",
        actionButton(button_id, if (running) paste(button_label, "running") else button_label, class = "btn-primary")
      )
    )
  })

	  output$mediation_main_content <- renderUI({
	    framework <- input$analysis_framework %||% "bayesian"
	    if (!isTRUE(mediation_result_available())) {
		      return(tagList(
		        tags$div(
		          class = "mediation-run-strip",
		          if (identical(framework, "bayesian")) {
		            tags$div(
		              id = "bayesian-client-progress",
		              class = "bayes-client-progress hidden",
		              tags$div(
		                class = "bayes-client-progress-head",
		                tags$span(class = "bayes-client-progress-label", "Starting Bayesian model fitting"),
		                tags$span(class = "bayes-client-progress-percent", "0%")
		              ),
		              tags$div(
		                class = "bayes-client-progress-track",
		                tags$div(class = "bayes-client-progress-fill")
		              )
		            )
		          }
		        ),
		        mediation_dag_placeholder(framework)
		      ))
		    }

    framework_label <- if (identical(framework, "bayesian")) "Bayesian posterior estimates" else "Frequentist estimates"

	    tagList(
	      tags$div(
	        class = "mediation-header",
	        tags$div(
	          h3("Mediation Results"),
	          tags$p(class = "setup-meta", framework_label)
	        )
	      ),
	      uiOutput("mediation_model_diagram")
	    )
	  })

  order_mediators_for_forest_plot <- function(df) {
    mediators_in_order <- df$mediator |>
      unique() |>
      sort()
    df$mediator <- factor(df$mediator, levels = rev(mediators_in_order))
    df
  }

	  output$bayesian_mediation_summary <- renderDT({
	    bayesian_mediation_summary() |>
	      present_bayesian_mediation_summary() |>
	      datatable(rownames = FALSE, options = list(pageLength = 5, dom = "tip", scrollX = FALSE))
	  }, server = FALSE)

	  output$bayesian_mediation_probabilities <- renderDT({
	    bayesian_mediation_summary() |>
	      present_bayesian_mediation_probabilities() |>
	      datatable(rownames = FALSE, options = list(pageLength = 5, dom = "t", scrollX = FALSE))
	  }, server = FALSE)

	  output$bayesian_mediation_path <- renderDT({
	    bayesian_mediation_path() |>
	      present_bayesian_mediation_path() |>
	      datatable(rownames = FALSE, options = list(pageLength = 10, dom = "tip", scrollX = FALSE))
	  }, server = FALSE)

	  output$bayesian_indirect_plot <- renderPlot({
	    df <- bayesian_mediation_path()
      validate(need(nrow(df) > 0, "No mediator-specific effects match the current filters."))
      df <- order_mediators_for_forest_plot(df)
      dodge <- position_dodge(width = 0.55)
	    ggplot(df, aes(x = IIE_est, y = mediator, color = contrast)) +
	      geom_vline(xintercept = 0, linewidth = 0.35, color = "grey60") +
	      geom_errorbarh(
	        aes(xmin = IIE_lwr, xmax = IIE_upr),
	        position = dodge,
	        height = 0.18,
	        linewidth = 0.55
	      ) +
	      geom_point(position = dodge, size = 2.3) +
	      labs(x = "Indirect Effect, posterior mean with 95% CrI", y = "Mediator", color = "Contrast") +
	      theme_bw(base_size = 12) +
        theme(
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = if (dplyr::n_distinct(df$contrast) > 1) "bottom" else "none"
        )
	  })

	  output$mediation_summary <- renderDT({
	    mediation_summary() |>
      present_frequentist_mediation_summary() |>
      datatable(rownames = FALSE, options = list(pageLength = 5, dom = "tip", scrollX = FALSE))
  }, server = FALSE)

  output$mediation_path <- renderDT({
    mediation_path() |>
      present_frequentist_mediation_path() |>
      datatable(rownames = FALSE, options = list(pageLength = 10, dom = "tip", scrollX = FALSE))
  }, server = FALSE)

  output$indirect_plot <- renderPlot({
    df <- mediation_path_with_bootstrap_ci()
    validate(need(nrow(df) > 0, "No mediator-specific effects match the current filters."))
    df <- order_mediators_for_forest_plot(df)
    dodge <- position_dodge(width = 0.55)
    p <- ggplot(df, aes(x = indirect_component, y = mediator, color = contrast)) +
      geom_vline(xintercept = 0, linewidth = 0.35, color = "grey60") +
      labs(
        x = if (!is.null(current_bootstrap_path_ci())) "Indirect Effect, estimate with 95% bootstrap CI" else "Indirect Effect",
        y = "Mediator",
        color = "Contrast"
      ) +
      theme_bw(base_size = 12) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = if (dplyr::n_distinct(df$contrast) > 1) "bottom" else "none"
      )
    if (all(c("indirect_lwr", "indirect_upr") %in% names(df)) && any(!is.na(df$indirect_lwr) & !is.na(df$indirect_upr))) {
      p <- p +
        geom_errorbarh(
          aes(xmin = indirect_lwr, xmax = indirect_upr),
          position = dodge,
          height = 0.18,
          linewidth = 0.55
        )
    }
    p + geom_point(position = dodge, size = 2.4)
  })

	  output$bootstrap_status <- renderUI({
	    boot <- bootstrap_result_or_null()
	    if (is.null(boot)) {
	      return(tags$div(class = "setup-meta", "Bootstrap CI will run automatically after frequentist results are ready."))
	    } else if (!identical(boot$key, frequentist_current_cache_key())) {
	      return(tags$div(class = "setup-meta", style = "color:#8a5a00;", "Bootstrap CI is pending for the current setup."))
	    } else {
	      tags$div(class = "setup-meta", style = "color:#1b6e3c;", paste("95% CI ready with B =", boot$B))
	    }
	  })

  output$bootstrap_summary_ci <- renderDT({
    req(bootstrap_result())
    validate(need(identical(bootstrap_result()$key, frequentist_current_cache_key()), "Bootstrap CI is updating for the current setup."))
    bootstrap_result()$ci$summary_ci |>
      dplyr::filter((input$mediation_contrast_filter %||% "all") == "all" | contrast == (input$mediation_contrast_filter %||% "all")) |>
      present_bootstrap_summary_ci(contrast_labels = NULL) |>
      datatable(rownames = FALSE, options = list(pageLength = 5, dom = "tip", scrollX = FALSE))
  }, server = FALSE)

	  output$bootstrap_path_ci <- renderDT({
	    req(bootstrap_result())
	    validate(need(identical(bootstrap_result()$key, frequentist_current_cache_key()), "Bootstrap CI is updating for the current setup."))
	    bootstrap_result()$ci$path_ci |>
	      dplyr::filter(
	        (input$mediation_contrast_filter %||% "all") == "all" | contrast == (input$mediation_contrast_filter %||% "all"),
	        (input$mediation_mediator_filter %||% "all") == "all" | mediator == (input$mediation_mediator_filter %||% "all")
	      ) |>
	      present_bootstrap_path_ci(contrast_labels = NULL) |>
	      datatable(rownames = FALSE, options = list(pageLength = 10, dom = "tip", scrollX = FALSE))
	  }, server = FALSE)

  output$result_export_controls <- renderUI({
    if (!isTRUE(mediation_result_available())) return(NULL)
    tags$div(class = "setup-sidebar-section setup-sidebar-actions result-export-actions",
      tags$div(style = "display:flex;flex-direction:column;gap:8px;",
        downloadButton("download_result_tables", "Download Tables (.xlsx)", style = "width:100%;"),
        downloadButton("download_result_figures", "Download Figures (.zip)", style = "width:100%;")
      )
    )
  })

  result_export_bundle <- function() {
    req(mediation_result_available())
    framework <- input$analysis_framework %||% "bayesian"
    label <- if (framework == "bayesian") "Bayesian" else "Frequentist"
    a <- analysis()
    scenario <- input$sensitivity_scenario %||% "common"
    rho <- input$sensitivity_active_rho %||% 0
    mediators <- sensitivity_mediator_ids()
    target <- input$sensitivity_one_mediator %||% mediators[[1]]
    if (!target %in% mediators) target <- mediators[[1]]
    normalize_pair <- function(pair, sensitivity = FALSE) {
      dplyr::bind_rows(export_normalize(pair$summary, framework, "effects", sensitivity),
                       export_normalize(pair$path, framework, "paths", sensitivity))
    }
    if (framework == "bayesian") {
      result <- bayesian_current_result()
      artifact <- result$artifact
      effects <- export_normalize(get_bayesian_mediation_summary(artifact), framework, "effects")
      paths <- export_normalize(get_bayesian_mediation_path(artifact), framework, "paths")
      common <- bayesian_sensitivity_curve_cache(artifact, result$key, "common")
      one <- bayesian_sensitivity_all_mediator_curve_cache(artifact, result$key)
      curves <- bind_rows(normalize_pair(common, TRUE), normalize_pair(one, TRUE))
      selected <- normalize_pair(bayesian_sensitivity_selected_cache(artifact, result$key, rho, scenario,
                                      if (scenario == "one_at_a_time") target else NULL), TRUE)
    } else {
      fit <- frequentist_current_result()$result$fit
      effects <- export_normalize(export_join_ci(fit$summary, current_bootstrap_summary_ci(), "contrast"), framework, "effects")
      paths <- export_normalize(export_join_ci(fit$path, current_bootstrap_path_ci(), c("contrast", "mediator")), framework, "paths")
      sens <- sensitivity()
      ci <- current_bootstrap_sensitivity_ci()
      curve_rows <- lapply(c("common", "one_at_a_time"), function(mode) {
        prefix <- if (mode == "common") "common" else "one"
        summary <- sens[[mode]]$sensitivity_summary |> mutate(contrast = contrast_label_current(exposure))
        path <- sens[[mode]]$sensitivity_path |> mutate(contrast = contrast_label_current(exposure))
        keys <- c("contrast", "rho", if (mode == "one_at_a_time") "varied_mediator")
        summary <- export_join_ci(summary, ci[[paste0(prefix, "_summary_ci")]], keys)
        path <- export_join_ci(path, ci[[paste0(prefix, "_path_ci")]], c(keys, "mediator"))
        normalize_pair(list(summary = summary, path = path), TRUE)
      })
      curves <- bind_rows(curve_rows)
      selected <- curves |> filter(.data$scenario == .env$scenario, abs(.data$rho - .env$rho) < 1e-10)
      if (scenario == "one_at_a_time") selected <- selected |> filter(varied_mediator == mediator_label(target))
    }
    contrasts <- unique(effects$contrast)
    choose_contrast <- function(value) if (!is.null(value) && value %in% contrasts) value else contrasts[[1]]
    selection <- list(mediation_contrast = choose_contrast(input$mediation_contrast_filter),
                      sensitivity_contrast = choose_contrast(input$sensitivity_contrast_filter),
                      mediator = input$mediation_mediator_filter %||% "all", scenario = scenario,
                      varied_mediator = mediator_label(target), rho = rho)
    if (selection$mediator != "all") selection$mediator <- mediator_label(selection$mediator)
    if (!selection$mediator %in% c("all", paths$mediator)) selection$mediator <- "all"
    settings <- data.frame(
      Setting = c("Common correlation", "One mediator", "Computed correlation grid", "Selected scenario", "Selected correlation", "Selected mediator",
                  "Interval interpretation", "Unavailable intervals", "Feasibility", "Zero crossings", "Table coverage", "Figure selection"),
      Value = c("All mediator-outcome residual correlations change together.",
                "One mediator-outcome residual correlation changes; the others remain zero.",
                paste(sort(unique(curves$rho)), collapse = ", "),
                if (scenario == "common") "Common correlation" else "One mediator", as.character(rho),
                if (scenario == "common") "All mediators" else mediator_label(target),
                paste(export_interval_type(framework), "; Interval excludes zero marks intervals entirely above or below zero."),
                "Blank bounds mean an interval is unavailable; no model or bootstrap is rerun during export.",
                "Sensitivity intervals use feasible draws or bootstrap replicates. Retained fractions are not included in this workbook; consult application warnings before interpretation.",
                "Linear interpolation within the computed grid; no crossing does not rule out a crossing outside this range.",
                "All contrasts and confirmed mediator PCs. Selected Rho contains the selected scenario and correlation for all contrasts.",
                paste("Mediation:", selection$mediation_contrast, "; Sensitivity:", selection$sensitivity_contrast)),
      stringsAsFactors = FALSE)
    setup <- data.frame(
      Setting = c("Framework", "Exported at", "Data source", "Observations", "Mediator features", "Exposure", "Reference", "Outcome",
                  "Coordinates", "Additional covariates", "Confirmed PCs", "PC signs", "PCA method", "Estimate type", "Bootstrap replicates", "Loading features per direction"),
      Value = c(label, format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), a$spec$data_source, nrow(a$dat_full_allpc),
                a$spec$mediator_features$n_features, a$spec$exposure$variable, a$spec$exposure$reference, a$spec$outcome,
                paste(unlist(a$spec$coordinates), collapse = ", "), paste(unname(a$spec$covariates), collapse = ", "),
                paste(confirmed_selected_pc_ids(), collapse = ", "),
                paste(confirmed_selected_pc_ids(), confirmed_pc_sign_vector()[confirmed_selected_pc_ids()], sep = "=", collapse = ", "),
                a$pca$method %||% "exact", if (framework == "bayesian") "Posterior means" else "Point estimates",
                if (framework == "frequentist") as.character(bootstrap_result_or_null()$B %||% "Not available") else "Not applicable",
                as.character(input$loading_n %||% 10)),
      stringsAsFactors = FALSE)
    list(framework = framework, framework_label = label, outcome = a$spec$outcome,
         setup = setup, settings = settings, effects = effects, paths = paths,
         curves = curves, selected = selected, selection = selection)
  }

  result_context_figures <- function() {
    a <- analysis()
    pcs <- confirmed_selected_pc_ids()
    signs <- confirmed_pc_sign_vector()
    spatial <- orient_spatial_scores(a$spatial_inputs$plot_df, signs, pca_fit = a$pca$pca_fit)
    loadings <- orient_loading_data(a$loading_df, signs, pca_fit = a$pca$pca_fit)
    n <- as.integer(input$loading_n %||% 10)
    plots <- list(Spatial_Exposure.png = overview_exposure_plot(),
                  Spatial_Outcome.png = overview_outcome_plot(),
                  PCA_Explained_Variance.png = explained_variance_plot())
    # Export every confirmed mediator, not only the current inspection filter.
    for (pc in pcs) {
      index <- as.integer(sub("^PC", "", pc))
      plots[[paste0(pc, "_Spatial.png")]] <- plot_pc_spatial(
        spatial, pc, a$pca$var_explained[[index]], pc,
        a$spatial_inputs$category_polygon_sf, flip_y = a$spatial_inputs$flip_y
      )
      plots[[paste0(pc, "_Loadings.png")]] <- export_loading_plot(
        top_loading_table(loadings, pc, "positive", n),
        top_loading_table(loadings, pc, "negative", n), pc
      )
    }
    plots
  }

  export_filename <- function(extension) {
    paste0(if ((input$analysis_framework %||% "bayesian") == "bayesian") "Bayesian" else "Frequentist",
           "_", format(Sys.Date(), "%Y-%m-%d"), ".", extension)
  }
  output$download_result_tables <- downloadHandler(
    filename = function() export_filename("xlsx"),
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    content = function(file) {
      withProgress(message = "Preparing result tables", value = 0.1, {
        bundle <- isolate(result_export_bundle())
        setProgress(0.8)
        write_result_workbook(bundle, file)
      })
    }
  )
  output$download_result_figures <- downloadHandler(
    filename = function() export_filename("zip"), contentType = "application/zip",
    content = function(file) {
      withProgress(message = "Preparing figures", value = 0.1, {
        bundle <- isolate(result_export_bundle())
        setProgress(0.6)
        write_result_figures(bundle, file, context_figures = isolate(result_context_figures()))
      })
    }
  )

  output$sensitivity_workflow_tabs <- renderUI({
    scenario <- input$sensitivity_scenario %||% "common"
    selected <- isolate(input$sensitivity_result_view) %||% "effect"
    if (identical(selected, "comparison") && !identical(scenario, "one_at_a_time")) {
      selected <- "effect"
    }

    tabsetPanel(
      id = "sensitivity_result_view",
      type = "tabs",
      selected = selected,
      tabPanel("Indirect Effect Across ρ", value = "effect", uiOutput("sensitivity_effect_content")),
      tabPanel("Effects at Selected ρ", value = "decomposition", uiOutput("sensitivity_decomposition_content")),
      if (identical(scenario, "one_at_a_time")) {
        tabPanel("Compare Mediators", value = "comparison", uiOutput("sensitivity_comparison_content"))
      }
    )
  })

  sensitivity_unavailable_card <- function(framework) {
    app_panel_card(
      if (identical(framework, "bayesian")) "Bayesian mediation is not fitted yet" else "Frequentist mediation is not fitted yet",
      class = "setup-empty",
      tags$p("Run mediation for the current data and selected PCs to view sensitivity results.")
    )
  }

  output$sensitivity_effect_content <- renderUI({
    framework <- input$analysis_framework %||% "bayesian"
    if (!isTRUE(mediation_result_available())) {
      return(sensitivity_unavailable_card(framework))
    }
    scenario <- input$sensitivity_scenario %||% "common"
    primary_plot_id <- if (identical(scenario, "common")) {
      if (identical(framework, "bayesian")) "bayesian_sensitivity_curve" else "sensitivity_curve"
    } else {
      if (identical(framework, "bayesian")) "bayesian_sensitivity_target_curve" else "sensitivity_target_curve"
    }
    primary_plot_title <- if (identical(scenario, "common")) {
      "Total Indirect Effect"
    } else {
      mediator_ids <- sensitivity_mediator_ids()
      mediator <- input$sensitivity_one_mediator %||% if (length(mediator_ids) > 0) mediator_ids[[1]] else "PC1"
      paste(mediator_label(mediator), "Indirect Effect")
    }

    tagList(
      app_panel_card("Zero Crossing", uiOutput("sensitivity_tipping_summary")),
      app_panel_card(primary_plot_title, plotOutput(primary_plot_id, height = 430))
    )
  })

  output$sensitivity_decomposition_content <- renderUI({
    framework <- input$analysis_framework %||% "bayesian"
    if (!isTRUE(mediation_result_available())) {
      return(sensitivity_unavailable_card(framework))
    }
    warning_id <- if (identical(framework, "bayesian")) "bayesian_sensitivity_warnings" else "frequentist_sensitivity_warnings"
    tagList(
      tags$p(
        class = "setup-meta",
        if (identical(framework, "bayesian")) {
          paste("Choose ρ to see the corresponding paths and effects. Posterior means with 95% credible intervals.", interval_star_note())
        } else {
          paste("Choose ρ to see the corresponding paths and effects. Point estimates with 95% bootstrap intervals where available.", interval_star_note())
        }
      ),
      uiOutput("sensitivity_rho_control"),
      uiOutput("sensitivity_selected_decomposition"),
      uiOutput(warning_id)
    )
  })

  output$sensitivity_comparison_content <- renderUI({
    framework <- input$analysis_framework %||% "bayesian"
    if (!isTRUE(mediation_result_available())) {
      return(sensitivity_unavailable_card(framework))
    }
    app_panel_card(
      "All Mediators Comparison",
      tags$p(class = "setup-meta", "Each curve varies one mediator at a time and shows that mediator's own indirect effect."),
      plotOutput(
        if (identical(framework, "bayesian")) "bayesian_sensitivity_all_mediator_comparison" else "sensitivity_all_mediator_comparison",
        height = 430
      )
    )
  })

	  bayesian_sensitivity_selected <- reactive({
	    req(identical(input$analysis_framework, "bayesian"))
	    artifact <- bayesian_current_artifact()
	    validate(need(!is.null(artifact), "Run Bayesian Mediation first for the current data and selected PCs."))
      result <- bayesian_current_result()
	    scenario <- input$sensitivity_scenario %||% "common"
	    rho <- input$sensitivity_active_rho %||% 0

	    if (identical(scenario, "common")) {
	      bayesian_sensitivity_selected_cache(
	        artifact = artifact,
          artifact_key = result$key,
	        rho = rho,
	        mode = "common"
	      )
	    } else {
	      mediator <- input$sensitivity_one_mediator %||% bayes_A_bundle_mediators(artifact$sensitivity_A$residual_draws)[[1]]
	      bayesian_sensitivity_selected_cache(
	        artifact = artifact,
          artifact_key = result$key,
	        rho = rho,
	        mode = "one_at_a_time",
	        mediator = mediator
	      )
	    }
	  })

  sensitivity_path_selected <- reactive({
    req(identical(input$analysis_framework, "frequentist"))
    req(input$sensitivity_scenario, input$sensitivity_active_rho)

    path_rows <- if (identical(input$sensitivity_scenario, "common")) {
      sensitivity()$common$sensitivity_path |>
        filter(valid_parameter, abs(rho - input$sensitivity_active_rho) < 1e-10) |>
        mutate(contrast = contrast_label_current(exposure), mediator = mediator_label(mediator)) |>
        filter_sensitivity_contrast() |>
        transmute(
          contrast,
          rho,
          mediator,
          alpha = alpha_X_to_M,
          beta = beta_adjusted,
          indirect_component = indirect_adjusted
        )
    } else {
      req(input$sensitivity_one_mediator)
      sensitivity()$one_at_a_time$sensitivity_path |>
        filter(
          valid_parameter,
          varied_mediator == input$sensitivity_one_mediator,
          mediator == input$sensitivity_one_mediator,
          abs(rho - input$sensitivity_active_rho) < 1e-10
        ) |>
        mutate(contrast = contrast_label_current(exposure), mediator = mediator_label(mediator)) |>
        filter_sensitivity_contrast() |>
        transmute(
          contrast,
          rho,
          mediator,
          alpha = alpha_X_to_M,
          beta = beta_adjusted,
          indirect_component = indirect_adjusted
        )
    }
    ci <- current_bootstrap_sensitivity_ci()
    if (is.null(ci)) {
      return(path_rows)
    }
    path_ci <- if (identical(input$sensitivity_scenario, "common")) {
      ci$common_path_ci
    } else {
      ci$one_path_ci
    }
    path_rows |>
      left_join(
        path_ci |>
          mutate(mediator = mediator_label(mediator)) |>
          select(
            contrast,
            rho,
            mediator,
            alpha_lwr,
            alpha_upr,
            beta_lwr,
            beta_upr,
            indirect_lwr,
            indirect_upr
          ),
        by = c("contrast", "rho", "mediator")
      )
  })

  sensitivity_value_cell <- function(est, lwr = NA_real_, upr = NA_real_, percent = FALSE) {
    if (isTRUE(percent)) {
      est_txt <- paste0(format_number(100 * est, 1), interval_star(lwr, upr), "%")
      if (is.na(lwr) || is.na(upr)) {
        return(est_txt)
      }
      return(paste0(est_txt, " [", format_number(100 * lwr, 1), "%, ", format_number(100 * upr, 1), "%]"))
    }
    format_estimate_interval(est, lwr, upr)
  }

  sensitivity_mediator_table <- function(path_df, framework = c("bayesian", "frequentist")) {
    framework <- match.arg(framework)
    if (nrow(path_df) == 0) {
      return(tags$p(class = "setup-meta", "No mediator-specific values are feasible at the selected ρ."))
    }
    rows <- lapply(seq_len(nrow(path_df)), function(i) {
      if (identical(framework, "bayesian")) {
        tags$tr(
          tags$td(mediator_label(path_df$mediator[[i]])),
          tags$td(format_estimate_interval(path_df$IIE_est[[i]], path_df$IIE_lwr[[i]], path_df$IIE_upr[[i]]))
        )
      } else {
        indirect_value <- if ("indirect_component" %in% names(path_df)) {
          path_df$indirect_component[[i]]
        } else {
          path_df$indirect_est[[i]]
        }
        tags$tr(
          tags$td(path_df$mediator[[i]]),
          tags$td(format_number(indirect_value))
        )
      }
    })
    tags$table(
      class = "sensitivity-mini-table",
      tags$thead(tags$tr(tags$th("Mediator"), tags$th("Indirect Effect"))),
      tags$tbody(rows)
    )
  }

  first_or_empty <- function(x) {
    if (length(x) == 0 || is.null(x[[1]]) || is.na(x[[1]])) "" else x[[1]]
  }

  zero_crossing_rho <- function(rho, value) {
    ok <- is.finite(rho) & is.finite(value)
    rho <- rho[ok]
    value <- value[ok]
    if (length(rho) == 0) {
      return(NA_real_)
    }
    ord <- order(rho)
    rho <- rho[ord]
    value <- value[ord]

    exact <- which(abs(value) < .Machine$double.eps^0.5)
    if (length(exact) > 0) {
      return(rho[[exact[[1]]]])
    }
    if (length(rho) < 2) {
      return(NA_real_)
    }

    changes <- which(value[-length(value)] * value[-1] < 0)
    if (length(changes) == 0) {
      return(NA_real_)
    }
    i <- changes[[1]]
    rho[[i]] - value[[i]] * (rho[[i + 1]] - rho[[i]]) / (value[[i + 1]] - value[[i]])
  }

  sensitivity_primary_curve_data <- reactive({
    framework <- input$analysis_framework %||% "bayesian"
    scenario <- input$sensitivity_scenario %||% "common"

    if (identical(framework, "bayesian")) {
      if (identical(scenario, "common")) {
        return(bayesian_sensitivity_curve()$summary |>
          filter_sensitivity_contrast() |>
          transmute(contrast, rho, value = total_IIE_est))
      }
      target <- input$sensitivity_one_mediator %||% ""
      return(bayesian_sensitivity_curve()$path |>
        dplyr::filter(mediator == target) |>
        filter_sensitivity_contrast() |>
        transmute(contrast, rho, value = IIE_est))
    }

    if (identical(scenario, "common")) {
      return(sensitivity()$common$sensitivity_summary |>
        filter(valid_parameter) |>
        mutate(contrast = contrast_label_current(exposure)) |>
        filter_sensitivity_contrast() |>
        transmute(contrast, rho, value = NIE_adjusted))
    }

    target <- input$sensitivity_one_mediator %||% ""
    sensitivity()$one_at_a_time$sensitivity_path |>
      filter(valid_parameter, varied_mediator == target, mediator == target) |>
      mutate(contrast = contrast_label_current(exposure)) |>
      filter_sensitivity_contrast() |>
      transmute(contrast, rho, value = indirect_adjusted)
  })

  output$sensitivity_tipping_summary <- renderUI({
    curve <- sensitivity_primary_curve_data()
    validate(need(nrow(curve) > 0, "No sensitivity curve rows match the current filters."))

    tipping <- curve |>
      group_by(contrast) |>
      summarise(
        rho_zero_crossing = zero_crossing_rho(rho, value),
        crossing_in_grid = is.finite(rho_zero_crossing),
        .groups = "drop"
      )

    scenario <- input$sensitivity_scenario %||% "common"
    quantity <- if (identical(scenario, "common")) {
      "Total indirect effect = 0"
    } else {
      paste(mediator_label(input$sensitivity_one_mediator %||% ""), "indirect effect = 0")
    }

    tags$div(
      tags$p(class = "setup-meta", paste("Zero crossing for the primary curve:", quantity)),
      tags$div(
        class = "sensitivity-tipping-grid",
        lapply(seq_len(nrow(tipping)), function(i) {
          tags$div(
            class = "sensitivity-tipping-card",
            tags$div(class = "sensitivity-tipping-label", tipping$contrast[[i]]),
            tags$div(
              class = "sensitivity-tipping-value",
              if (isTRUE(tipping$crossing_in_grid[[i]])) {
                paste0("ρ = ", format_number(tipping$rho_zero_crossing[[i]], digits = 3))
              } else {
                "No crossing in grid"
              }
            )
          )
        })
      )
    )
  })

  sensitivity_rho_label <- function() {
    rho <- input$sensitivity_active_rho %||% 0
    paste0("ρ = ", format_number(as.numeric(rho), digits = 2))
  }

  sensitivity_summary_values <- function(row, framework) {
    if (identical(framework, "bayesian")) {
      return(list(
        direct_est = row$direct_est,
        direct_lwr = row$direct_lwr,
        direct_upr = row$direct_upr,
        indirect_est = row$total_IIE_est,
        indirect_lwr = row$total_IIE_lwr,
        indirect_upr = row$total_IIE_upr,
        total_est = row$TE_est,
        total_lwr = row$TE_lwr,
        total_upr = row$TE_upr
      ))
    }
    value_or_na <- function(name) {
      if (!name %in% names(row)) {
        return(NA_real_)
      }
      row[[name]]
    }
    list(
      direct_est = row$adjusted_direct_effect,
      direct_lwr = value_or_na("adjusted_direct_lwr"),
      direct_upr = value_or_na("adjusted_direct_upr"),
      indirect_est = row$adjusted_total_indirect,
      indirect_lwr = value_or_na("adjusted_indirect_lwr"),
      indirect_upr = value_or_na("adjusted_indirect_upr"),
      total_est = row$adjusted_total_effect,
      total_lwr = value_or_na("adjusted_total_lwr"),
      total_upr = value_or_na("adjusted_total_upr")
    )
  }

  sensitivity_dag_panel <- function(summary_row, path_rows, framework) {
    if (nrow(path_rows) == 0) {
      return(tags$p(class = "setup-meta", "No mediator component rows match the selected sensitivity setting."))
    }
    contrast <- summary_row$contrast[[1]]
    values <- sensitivity_summary_values(summary_row, framework)
    alpha_label <- mediation_dag_path_label("a", path_rows, framework, path = "alpha")
    beta_label <- mediation_dag_path_label("b", path_rows, framework, path = "beta")
    direct_tooltip <- mediation_dag_tooltip("Direct effect", values$direct_est, values$direct_lwr, values$direct_upr)
    direct_label <- tags$span(
      class = "mediation-path-symbol-wrap",
      tags$span(
        class = "interval-tooltip",
        `data-tooltip` = direct_tooltip,
        tabindex = "0",
        paste0("c' = ", mediation_dag_number(values$direct_est, values$direct_lwr, values$direct_upr))
      ),
      app_help_tooltip("c' is the direct exposure-to-outcome path after accounting for the selected mediator PCs under this ρ scenario.")
    )

    tags$div(
      class = "sensitivity-dag-block",
      tags$div(
        class = "sensitivity-dag-canvas",
        tags$svg(
          class = "mediation-arrow-svg",
          viewBox = "0 0 100 100",
          preserveAspectRatio = "none",
          xmlns = "http://www.w3.org/2000/svg",
          tags$defs(
            tags$marker(id = "sensitivity-arrow-alpha", markerWidth = "8", markerHeight = "8", refX = "7", refY = "4", orient = "auto", markerUnits = "strokeWidth",
                        tags$polygon(points = "0 0, 8 4, 0 8", fill = "#2f5597")),
            tags$marker(id = "sensitivity-arrow-beta", markerWidth = "8", markerHeight = "8", refX = "7", refY = "4", orient = "auto", markerUnits = "strokeWidth",
                        tags$polygon(points = "0 0, 8 4, 0 8", fill = "#7e57c2")),
            tags$marker(id = "sensitivity-arrow-direct", markerWidth = "8", markerHeight = "8", refX = "7", refY = "4", orient = "auto", markerUnits = "strokeWidth",
                        tags$polygon(points = "0 0, 8 4, 0 8", fill = "#6b7b8f"))
          ),
          tags$line(class = "mediation-dag-line mediation-dag-line-alpha", x1 = "25", y1 = "78", x2 = "44", y2 = "30", `marker-end` = "url(#sensitivity-arrow-alpha)"),
          tags$line(class = "mediation-dag-line mediation-dag-line-beta", x1 = "56", y1 = "30", x2 = "75", y2 = "78", `marker-end` = "url(#sensitivity-arrow-beta)"),
          tags$line(class = "mediation-dag-line mediation-dag-line-direct", x1 = "25", y1 = "82", x2 = "75", y2 = "82", `marker-end` = "url(#sensitivity-arrow-direct)")
        ),
        tags$div(
          class = "mediation-arrow-label mediation-arrow-label-xm",
          alpha_label
        ),
        tags$div(
          class = "mediation-arrow-label mediation-arrow-label-my",
          beta_label
        ),
        tags$div(
          class = "mediation-arrow-label mediation-arrow-label-xy",
          title = direct_tooltip,
          direct_label
        ),
        tags$div(class = "sensitivity-rho-box", sensitivity_rho_label()),
        tags$div(
          class = "mediation-node mediation-node-x",
          tags$div(class = "mediation-node-label", "Exposure"),
          tags$div(class = "mediation-node-value", "X"),
          tags$div(class = "mediation-node-detail", contrast)
        ),
        tags$div(
          class = "mediation-node mediation-node-m",
          tags$div(class = "mediation-node-label", "Mediator"),
          tags$div(class = "mediation-node-value", "M"),
          tags$div(class = "mediation-node-detail", paste(unique(path_rows$mediator), collapse = ", "))
        ),
        tags$div(
          class = "mediation-node mediation-node-y",
          tags$div(class = "mediation-node-label", "Outcome"),
          tags$div(class = "mediation-node-value", "Y"),
          tags$div(class = "mediation-node-detail", "Outcome")
        )
      ),
      mediation_dag_effect_grid(
        path_rows,
        framework,
        values$indirect_est,
        values$indirect_lwr,
        values$indirect_upr,
        values$direct_est,
        values$direct_lwr,
        values$direct_upr,
        values$total_est,
        values$total_lwr,
        values$total_upr
      )
    )
  }

  output$sensitivity_selected_decomposition <- renderUI({
    framework <- input$analysis_framework %||% "bayesian"
    scenario <- input$sensitivity_scenario %||% "common"

    if (identical(framework, "bayesian")) {
      selected <- bayesian_sensitivity_selected()
      summary <- filter_sensitivity_contrast(selected$summary)
      path <- filter_sensitivity_contrast(selected$path)
      if (identical(scenario, "one_at_a_time")) {
        target <- input$sensitivity_one_mediator %||% first_or_empty(selected$path$mediator_varied)
        path <- dplyr::filter(path, mediator == target)
      }
      path <- path |>
        dplyr::mutate(
          mediator_id = mediator,
          mediator_varied_id = mediator_varied,
          mediator = mediator_label(mediator),
          mediator_varied = mediator_label(mediator_varied)
        )

      return(tagList(lapply(seq_len(nrow(summary)), function(i) {
        row <- summary[i, ]
        path_rows <- dplyr::filter(path, contrast == row$contrast)
        app_panel_card(
          row$contrast,
          class = "sensitivity-decomp-panel",
          sensitivity_dag_panel(row, path_rows, framework = "bayesian")
        )
      })))
    }

    summary <- sensitivity_selected()
    path <- sensitivity_path_selected()
    tagList(lapply(seq_len(nrow(summary)), function(i) {
      row <- summary[i, ]
      path_rows <- dplyr::filter(path, contrast == row$contrast)
      app_panel_card(
        row$contrast,
        class = "sensitivity-decomp-panel",
        sensitivity_dag_panel(row, path_rows, framework = "frequentist")
      )
    }))
  })

	  bayesian_sensitivity_curve <- reactive({
	    req(identical(input$analysis_framework, "bayesian"))
	    artifact <- bayesian_current_artifact()
	    validate(need(!is.null(artifact), "Run Bayesian Mediation first for the current data and selected PCs."))
	    result <- bayesian_current_result()
	    scenario <- input$sensitivity_scenario %||% "common"

	    if (identical(scenario, "common")) {
	      bayesian_sensitivity_curve_cache(artifact = artifact, artifact_key = result$key, mode = "common")
	    } else {
	      mediator <- input$sensitivity_one_mediator %||% bayes_A_bundle_mediators(artifact$sensitivity_A$residual_draws)[[1]]
	      bayesian_sensitivity_curve_cache(artifact = artifact, artifact_key = result$key, mode = "one_at_a_time", mediator = mediator)
	    }
	  })

	  output$bayesian_sensitivity_summary <- renderDT({
	    bayesian_sensitivity_selected()$summary |>
	      present_bayesian_sensitivity_summary() |>
	      datatable(rownames = FALSE, options = list(pageLength = 5, dom = "tip", scrollX = FALSE))
	  }, server = FALSE)

	  output$bayesian_sensitivity_path <- renderDT({
	    bayesian_sensitivity_selected()$path |>
	      present_bayesian_sensitivity_path() |>
	      datatable(rownames = FALSE, options = list(pageLength = 10, dom = "tip", scrollX = FALSE))
	  }, server = FALSE)

	  output$bayesian_sensitivity_warnings <- renderUI({
	    selected <- bayesian_sensitivity_selected()
	    warnings <- list()
	    feasible_fraction <- min(selected$summary$valid_draw_fraction, na.rm = TRUE)

	    if (is.finite(feasible_fraction) && feasible_fraction <= 0.25) {
	      warnings <- c(
	        warnings,
	        list(tags$p("This sensitivity setting is poorly feasible for the posterior covariance structure; do not interpret it as an ordinary valid result."))
	      )
	    } else if (is.finite(feasible_fraction) && feasible_fraction < 1) {
	      warnings <- c(
	        warnings,
	        list(tags$p("This sensitivity setting is not feasible for all posterior draws; interpret the result with caution."))
	      )
	    }

	    if (length(warnings) == 0) {
	      return(NULL)
	    }

	    tags$div(style = "color:#8a5a00;", warnings)
	  })

	  output$bayesian_sensitivity_curve <- renderPlot({
	    curve <- bayesian_sensitivity_curve()$summary |>
        filter_sensitivity_contrast()
	    active_rho <- input$sensitivity_active_rho %||% 0
      validate(need(nrow(curve) > 0, "No sensitivity curve rows match the current contrast filter."))

	    ggplot(curve, aes(x = rho, y = total_IIE_est, color = contrast, fill = contrast)) +
	      geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
	      geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
	      geom_vline(xintercept = active_rho, linetype = "dashed", color = "grey45") +
	      geom_ribbon(aes(ymin = total_IIE_lwr, ymax = total_IIE_upr), alpha = 0.12, color = NA) +
	      geom_line(linewidth = 0.9) +
	      geom_point(size = 1.4) +
	      labs(x = "Residual correlation (ρ)", y = "Total Indirect Effect, posterior mean with 95% CrI", color = "Contrast", fill = "Contrast") +
	      theme_bw(base_size = 12)
	  })

    output$bayesian_sensitivity_target_curve <- renderPlot({
      req(identical(input$sensitivity_scenario, "one_at_a_time"))
      target <- input$sensitivity_one_mediator %||% ""
      curve <- bayesian_sensitivity_curve()$path |>
        dplyr::filter(mediator == target) |>
        filter_sensitivity_contrast()
      active_rho <- input$sensitivity_active_rho %||% 0
      validate(need(nrow(curve) > 0, "No target mediator sensitivity curve rows match the current filters."))

      ggplot(curve, aes(x = rho, y = IIE_est, color = contrast, fill = contrast)) +
        geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
        geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
        geom_vline(xintercept = active_rho, linetype = "dashed", color = "grey45") +
        geom_ribbon(aes(ymin = IIE_lwr, ymax = IIE_upr), alpha = 0.12, color = NA) +
        geom_line(linewidth = 0.9) +
        geom_point(size = 1.4) +
        labs(
          x = paste("Residual correlation for", mediator_label(target)),
          y = paste(mediator_label(target), "Indirect Effect, posterior mean with 95% CrI"),
          color = "Contrast",
          fill = "Contrast"
        ) +
        theme_bw(base_size = 12)
    })

		  output$bayesian_sensitivity_path_curve <- renderPlot({
		    curve <- bayesian_sensitivity_curve()$path |>
	        filter_sensitivity_contrast()
	    selected_varied <- input$sensitivity_one_mediator %||% (unique(curve$mediator)[[1]] %||% "")
	    curve <- curve |>
	      mutate(
	        mediator_id = mediator,
	        mediator = mediator_label(mediator),
	        role = if_else(
	          identical(input$sensitivity_scenario, "one_at_a_time") & mediator_id == selected_varied,
	          "Varied mediator",
	          "Other mediators"
	        )
	      )
	    active_rho <- input$sensitivity_active_rho %||% 0
      validate(need(nrow(curve) > 0, "No mediator-specific sensitivity curve rows match the current filters."))

	    ggplot(curve, aes(x = rho, y = IIE_est, color = mediator, fill = mediator)) +
	      geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
	      geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
	      geom_vline(xintercept = active_rho, linetype = "dashed", color = "grey45") +
	      geom_ribbon(aes(ymin = IIE_lwr, ymax = IIE_upr), alpha = 0.10, color = NA) +
	      geom_line(aes(linetype = role, alpha = role), linewidth = 0.8) +
	      scale_alpha_manual(values = c("Varied mediator" = 1, "Other mediators" = 0.55), guide = "none") +
	      facet_wrap(~ contrast, scales = "free_y") +
		      labs(x = "Residual correlation (ρ)", y = "Indirect Effect, posterior mean with 95% CrI", color = "Mediator", fill = "Mediator", linetype = "Role") +
		      theme_bw(base_size = 12)
		  })

	  output$bayesian_sensitivity_all_mediator_comparison <- renderPlot({
	    req(identical(input$analysis_framework, "bayesian"))
	    req(identical(input$sensitivity_scenario, "one_at_a_time"))
	    artifact <- bayesian_current_artifact()
	    validate(need(!is.null(artifact), "Run Bayesian Mediation first for the current data and selected PCs."))
	    result <- bayesian_current_result()
	    curve <- bayesian_sensitivity_all_mediator_curve_cache(
	      artifact = artifact,
	      artifact_key = result$key
	    )$path |>
	      dplyr::filter(mediator == mediator_varied) |>
	      filter_sensitivity_contrast() |>
	      mutate(mediator = mediator_label(mediator))
	    validate(need(nrow(curve) > 0, "No mediator comparison rows match the current contrast filter."))

	    ggplot(curve, aes(x = rho, y = IIE_est, color = mediator, fill = mediator)) +
	      geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
	      geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
	      geom_vline(xintercept = input$sensitivity_active_rho %||% 0, linetype = "dashed", color = "grey45") +
	      geom_ribbon(aes(ymin = IIE_lwr, ymax = IIE_upr), alpha = 0.08, color = NA) +
	      geom_line(linewidth = 0.85) +
	      geom_point(size = 1.2) +
	      facet_wrap(~ contrast, scales = "free_y") +
	      labs(
	        x = "Residual correlation for each mediator varied one at a time",
	        y = "Mediator-specific Indirect Effect, posterior mean with 95% CrI",
	        color = "Mediator",
	        fill = "Mediator"
	      ) +
	      theme_bw(base_size = 12)
	  })

		  output$rho0_validation <- renderPrint({
	    s <- sensitivity()
    list(
      common_rho0_matches_main = s$rho0_common$passed,
      one_at_a_time_rho0_matches_main = s$rho0_one_at_a_time$passed
    )
  })

  sensitivity_selected <- reactive({
    req(identical(input$analysis_framework, "frequentist"))
    req(input$sensitivity_scenario, input$sensitivity_active_rho)

    rows <- if (identical(input$sensitivity_scenario, "common")) {
      sensitivity()$common$sensitivity_summary |>
        filter(valid_parameter, abs(rho - input$sensitivity_active_rho) < 1e-10) |>
        mutate(contrast = contrast_label_current(exposure)) |>
        filter_sensitivity_contrast() |>
        transmute(
          contrast,
          rho,
          adjusted_total_indirect = NIE_adjusted,
          adjusted_direct_effect = NDE_adjusted,
          adjusted_total_effect = NDE_adjusted + NIE_adjusted,
          PM_adjusted
        )
    } else {
      req(input$sensitivity_one_mediator)
      sensitivity()$one_at_a_time$sensitivity_summary |>
        filter(
          valid_parameter,
          varied_mediator == input$sensitivity_one_mediator,
          abs(rho - input$sensitivity_active_rho) < 1e-10
        ) |>
        mutate(contrast = contrast_label_current(exposure)) |>
        filter_sensitivity_contrast() |>
        transmute(
          contrast,
          varied_mediator,
          rho,
          adjusted_total_indirect = NIE_adjusted,
          adjusted_direct_effect = NDE_adjusted,
          adjusted_total_effect = NDE_adjusted + NIE_adjusted,
          PM_adjusted
        )
    }
    ci <- current_bootstrap_sensitivity_ci()
    if (is.null(ci)) {
      return(rows)
    }
    summary_ci <- if (identical(input$sensitivity_scenario, "common")) {
      ci$common_summary_ci
    } else {
      ci$one_summary_ci
    }
    join_cols <- c("contrast", "rho")
    if ("varied_mediator" %in% names(rows)) {
      join_cols <- c(join_cols, "varied_mediator")
    }
    rows |>
      left_join(
        summary_ci |>
          select(
            contrast,
            varied_mediator,
            rho,
            adjusted_direct_lwr,
            adjusted_direct_upr,
            adjusted_indirect_lwr,
            adjusted_indirect_upr,
            adjusted_total_lwr,
            adjusted_total_upr
          ),
        by = join_cols
      )
  })

  output$sensitivity_selected <- renderDT({
    sensitivity_selected() |>
      present_frequentist_sensitivity_selected() |>
      datatable(rownames = FALSE, options = list(pageLength = 5, dom = "t"))
  }, server = FALSE)

  output$frequentist_sensitivity_warnings <- renderUI({
    selected <- sensitivity_selected()
    warnings <- list()

    if (nrow(selected) == 0) {
      warnings <- c(warnings, list(tags$p("The selected ρ is not feasible for the current frequentist residual covariance structure.")))
    }
    if (length(warnings) == 0) {
      return(NULL)
    }
    tags$div(style = "color:#8a5a00;", warnings)
  })

  output$sensitivity_curve <- renderPlot({
    req(identical(input$analysis_framework, "frequentist"))
    req(input$sensitivity_scenario)

    if (identical(input$sensitivity_scenario, "common")) {
      curve <- sensitivity()$common$sensitivity_summary |>
        filter(valid_parameter) |>
        mutate(contrast = contrast_label_current(exposure)) |>
        filter_sensitivity_contrast() |>
        transmute(contrast, rho, NIE_adjusted)
      ci <- current_bootstrap_sensitivity_ci()
      if (!is.null(ci)) {
        curve <- curve |>
          left_join(
            ci$common_summary_ci |>
              select(contrast, rho, adjusted_indirect_lwr, adjusted_indirect_upr),
            by = c("contrast", "rho")
          )
      }
      p <- ggplot(curve, aes(x = rho, y = NIE_adjusted, color = contrast)) +
        geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
        geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
        geom_vline(xintercept = input$sensitivity_active_rho %||% 0, linetype = "dashed", color = "grey45")
      if (all(c("adjusted_indirect_lwr", "adjusted_indirect_upr") %in% names(curve))) {
        p <- p + geom_ribbon(aes(ymin = adjusted_indirect_lwr, ymax = adjusted_indirect_upr, fill = contrast), alpha = 0.10, color = NA)
      }
      p +
        geom_line(linewidth = 0.9) +
        geom_point(size = 1.4) +
        labs(x = "Residual correlation (ρ)", y = "Total Indirect Effect", color = "Contrast", fill = "Contrast") +
        theme_bw(base_size = 12)
    } else {
      req(input$sensitivity_one_mediator)
      curve <- sensitivity()$one_at_a_time$sensitivity_summary |>
        filter(valid_parameter, varied_mediator == input$sensitivity_one_mediator) |>
        mutate(contrast = contrast_label_current(exposure)) |>
        filter_sensitivity_contrast() |>
        transmute(contrast, varied_mediator, rho, NIE_adjusted)
      ci <- current_bootstrap_sensitivity_ci()
      if (!is.null(ci)) {
        curve <- curve |>
          left_join(
            ci$one_summary_ci |>
              select(contrast, varied_mediator, rho, adjusted_indirect_lwr, adjusted_indirect_upr),
            by = c("contrast", "varied_mediator", "rho")
          )
      }
      p <- ggplot(curve, aes(x = rho, y = NIE_adjusted, color = contrast)) +
        geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
        geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
        geom_vline(xintercept = input$sensitivity_active_rho %||% 0, linetype = "dashed", color = "grey45")
      if (all(c("adjusted_indirect_lwr", "adjusted_indirect_upr") %in% names(curve))) {
        p <- p + geom_ribbon(aes(ymin = adjusted_indirect_lwr, ymax = adjusted_indirect_upr, fill = contrast), alpha = 0.10, color = NA)
      }
      p +
        geom_line(linewidth = 0.9) +
        geom_point(size = 1.4) +
        labs(x = paste("Residual correlation for", mediator_label(input$sensitivity_one_mediator)), y = "Total Indirect Effect", color = "Contrast", fill = "Contrast") +
        theme_bw(base_size = 12)
    }
  })

  output$sensitivity_target_curve <- renderPlot({
    req(identical(input$analysis_framework, "frequentist"))
    req(identical(input$sensitivity_scenario, "one_at_a_time"))
    req(input$sensitivity_one_mediator)

    target <- input$sensitivity_one_mediator
    curve <- sensitivity()$one_at_a_time$sensitivity_path |>
      filter(
        valid_parameter,
        varied_mediator == target,
        mediator == target
      ) |>
      mutate(contrast = contrast_label_current(exposure)) |>
      filter_sensitivity_contrast()
    ci <- current_bootstrap_sensitivity_ci()
    if (!is.null(ci)) {
      curve <- curve |>
        mutate(mediator = mediator_label(mediator), varied_mediator = mediator_label(varied_mediator)) |>
        left_join(
          ci$one_path_ci |>
            mutate(mediator = mediator_label(mediator), varied_mediator = mediator_label(varied_mediator)) |>
            select(contrast, varied_mediator, mediator, rho, indirect_lwr, indirect_upr),
          by = c("contrast", "varied_mediator", "mediator", "rho")
        )
    }
    validate(need(nrow(curve) > 0, "No target mediator sensitivity curve rows match the current filters."))

    p <- ggplot(curve, aes(x = rho, y = indirect_adjusted, color = contrast)) +
      geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
      geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
      geom_vline(xintercept = input$sensitivity_active_rho %||% 0, linetype = "dashed", color = "grey45")
    if (all(c("indirect_lwr", "indirect_upr") %in% names(curve))) {
      p <- p + geom_ribbon(aes(ymin = indirect_lwr, ymax = indirect_upr, fill = contrast), alpha = 0.10, color = NA)
    }
    p +
      geom_line(linewidth = 0.9) +
      geom_point(size = 1.4) +
      labs(
        x = paste("Residual correlation for", mediator_label(target)),
        y = paste(mediator_label(target), "Indirect Effect"),
        color = "Contrast",
        fill = "Contrast"
      ) +
      theme_bw(base_size = 12)
  })

	  output$sensitivity_path_curve <- renderPlot({
	    req(identical(input$analysis_framework, "frequentist"))
	    req(input$sensitivity_scenario)

    if (identical(input$sensitivity_scenario, "common")) {
      sensitivity()$common$sensitivity_path |>
        filter(valid_parameter) |>
        mutate(contrast = contrast_label_current(exposure), mediator = mediator_label(mediator)) |>
        filter_sensitivity_contrast() |>
        ggplot(aes(x = rho, y = indirect_adjusted, color = mediator)) +
        geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
        geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
        geom_vline(xintercept = input$sensitivity_active_rho %||% 0, linetype = "dashed", color = "grey45") +
        geom_line(linewidth = 0.8) +
        facet_wrap(~ contrast, scales = "free_y") +
        labs(x = "Residual correlation (ρ)", y = "Indirect Effect", color = "Mediator") +
        theme_bw(base_size = 12)
    } else {
      req(input$sensitivity_one_mediator)
      sensitivity()$one_at_a_time$sensitivity_path |>
        filter(valid_parameter, varied_mediator == input$sensitivity_one_mediator) |>
        mutate(
          contrast = contrast_label_current(exposure),
          mediator_id = mediator,
          mediator = mediator_label(mediator),
          role = if_else(mediator_id == input$sensitivity_one_mediator, "Varied mediator", "Other mediators")
        ) |>
        filter_sensitivity_contrast() |>
        ggplot(aes(x = rho, y = indirect_adjusted, color = mediator)) +
        geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
        geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
        geom_vline(xintercept = input$sensitivity_active_rho %||% 0, linetype = "dashed", color = "grey45") +
        geom_line(aes(linetype = role, alpha = role), linewidth = 0.8) +
        scale_alpha_manual(values = c("Varied mediator" = 1, "Other mediators" = 0.55), guide = "none") +
        facet_wrap(~ contrast, scales = "free_y") +
	        labs(x = paste("Residual correlation for", mediator_label(input$sensitivity_one_mediator)), y = "Indirect Effect", color = "Mediator", linetype = "Role") +
	        theme_bw(base_size = 12)
	    }
	  })

	  output$sensitivity_all_mediator_comparison <- renderPlot({
	    req(identical(input$analysis_framework, "frequentist"))
	    req(identical(input$sensitivity_scenario, "one_at_a_time"))

	    curve <- sensitivity()$one_at_a_time$sensitivity_path |>
	      filter(valid_parameter, mediator == varied_mediator) |>
	      mutate(
	        contrast = contrast_label_current(exposure),
	        mediator = mediator_label(mediator),
	        varied_mediator = mediator_label(varied_mediator)
	      ) |>
	      filter_sensitivity_contrast()
	    ci <- current_bootstrap_sensitivity_ci()
	    if (!is.null(ci)) {
	      curve <- curve |>
	        left_join(
	          ci$one_path_ci |>
	            mutate(mediator = mediator_label(mediator), varied_mediator = mediator_label(varied_mediator)) |>
	            select(contrast, mediator, varied_mediator, rho, indirect_lwr, indirect_upr),
	          by = c("contrast", "mediator", "varied_mediator", "rho")
	        )
	    }
	    validate(need(nrow(curve) > 0, "No mediator comparison rows match the current contrast filter."))

	    p <- ggplot(curve, aes(x = rho, y = indirect_adjusted, color = mediator)) +
	      geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
	      geom_vline(xintercept = 0, linewidth = 0.3, color = "grey70") +
	      geom_vline(xintercept = input$sensitivity_active_rho %||% 0, linetype = "dashed", color = "grey45")
	    if (all(c("indirect_lwr", "indirect_upr") %in% names(curve))) {
	      p <- p + geom_ribbon(aes(ymin = indirect_lwr, ymax = indirect_upr, fill = mediator), alpha = 0.08, color = NA)
	    }
	    p +
	      geom_line(linewidth = 0.85) +
	      geom_point(size = 1.2) +
	      facet_wrap(~ contrast, scales = "free_y") +
	      labs(
	        x = "Residual correlation for each mediator varied one at a time",
	        y = "Mediator-specific Indirect Effect",
	        color = "Mediator",
	        fill = "Mediator"
	      ) +
	      theme_bw(base_size = 12)
	  })

	  output$sensitivity_tipping <- renderDT({
	    req(identical(input$analysis_framework, "frequentist"))
    req(input$sensitivity_scenario)

    tipping_tbl <- if (identical(input$sensitivity_scenario, "common")) {
      sensitivity()$tipping_common$total |>
        mutate(contrast = contrast_label_current(exposure)) |>
        filter_sensitivity_contrast() |>
        select(contrast, rho_zero_crossing, crossing_in_grid)
    } else {
      req(input$sensitivity_one_mediator)
      sensitivity()$tipping_one_at_a_time$total |>
        filter(varied_mediator == input$sensitivity_one_mediator) |>
        mutate(contrast = contrast_label_current(exposure), varied_mediator = mediator_label(varied_mediator)) |>
        filter_sensitivity_contrast() |>
        select(contrast, varied_mediator, rho_zero_crossing, crossing_in_grid)
    }

    tipping_tbl |>
      rename(
	        `ρ zero crossing` = rho_zero_crossing,
        `crossing in grid` = crossing_in_grid
      ) |>
      round_numeric(4) |>
      datatable(rownames = FALSE, options = list(pageLength = 5, dom = "t"))
  }, server = FALSE)
}

shinyApp(ui = ui, server = server)
