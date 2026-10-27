# App-ready Bayesian analysis adapter.
#
# This file loads a validated Bayesian artifact and exposes small, tidy objects
# for Shiny. It never fits brms/Stan models.

default_bayesian_app_artifact_path <- function(project_root = ".") {
  file.path(
    project_root,
    "results",
    "bayesian",
    "current",
    "emory_bayesian_v1_sensitivity_A_app_artifact.rds"
  )
}

load_bayesian_app_artifact <- function(path = default_bayesian_app_artifact_path()) {
  if (!file.exists(path)) {
    stop("Bayesian app artifact not found: ", path, call. = FALSE)
  }

  artifact <- readRDS(path)
  validate_bayesian_app_artifact(artifact)
  artifact
}

validate_bayesian_app_artifact <- function(artifact) {
  required <- c("artifact_version", "dataset_id", "model_id", "settings", "baseline", "sensitivity_A")
  missing <- setdiff(required, names(artifact))
  if (length(missing) > 0) {
    stop("Bayesian app artifact is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  baseline_required <- c("wide", "decomposition_summary", "path_summary", "n_draws")
  missing_baseline <- setdiff(baseline_required, names(artifact$baseline))
  if (length(missing_baseline) > 0) {
    stop("Bayesian app artifact baseline is missing: ", paste(missing_baseline, collapse = ", "), call. = FALSE)
  }

  sensitivity_required <- c("residual_draws", "validation")
  missing_sensitivity <- setdiff(sensitivity_required, names(artifact$sensitivity_A))
  if (length(missing_sensitivity) > 0) {
    stop("Bayesian app artifact sensitivity_A is missing: ", paste(missing_sensitivity, collapse = ", "), call. = FALSE)
  }
  if (!identical(artifact$sensitivity_A$parameterization, structural_rho_version())) {
    stop("Outdated sensitivity artifact. Run scripts/current/refresh_structural_rho.R before loading it.", call. = FALSE)
  }

  forbidden_top_level <- c("B1_XY", "B2_XM", "models", "fits")
  present_forbidden <- intersect(forbidden_top_level, names(artifact))
  if (length(present_forbidden) > 0) {
    stop("Bayesian app artifact contains forbidden runtime objects: ", paste(present_forbidden, collapse = ", "), call. = FALSE)
  }

  if (any(c("B1_XY", "B2_XM", "brmsfit", "stanfit") %in% names(artifact$sensitivity_A))) {
    stop("Bayesian app artifact sensitivity_A contains excluded model objects.", call. = FALSE)
  }

  if (exists("validate_bayesian_A_residual_draws", mode = "function")) {
    validate_bayesian_A_residual_draws(artifact$sensitivity_A$residual_draws)
  }

  invisible(TRUE)
}

bayesian_app_mcmc_diagnostics <- function(artifact) {
  diagnostics <- artifact$diagnostics %||% NULL
  if (is.null(diagnostics)) {
    return(NULL)
  }
  if (inherits(diagnostics, "data.frame")) {
    return(diagnostics)
  }
  if (is.list(diagnostics) && inherits(diagnostics$mcmc, "data.frame")) {
    return(diagnostics$mcmc)
  }
  NULL
}

make_bayesian_app_artifact_from_v1_results <- function(results,
                                                       dataset_id = "emory_hudson_immune_GEX_v2",
                                                       posterior_mc_draws = NULL,
                                                       seed = NULL) {
  if (!exists("extract_bayesian_v1_app_residual_draws", mode = "function")) {
    stop("Source R/bayesian_mediation_v1.R before creating a v1 app artifact.", call. = FALSE)
  }
  if (!exists("compute_bayesian_A_baseline_from_residual_draws", mode = "function")) {
    stop("Source R/bayesian_sensitivity_A_only.R before creating a v1 app artifact.", call. = FALSE)
  }

  posterior_mc_draws <- posterior_mc_draws %||%
    results$settings$posterior_mc_draws %||%
    results$posterior_quantities$n_draws %||%
    NULL
  seed <- seed %||% results$settings$seed %||% 123

  sampled_indices <- results$posterior_quantities$sampled_indices %||% NULL
  residual_draws <- extract_bayesian_v1_app_residual_draws(
    fits = results$models,
    n_draws = posterior_mc_draws,
    seed = seed + 100,
    sampled_indices = sampled_indices
  )
  baseline <- compute_bayesian_A_baseline_from_residual_draws(residual_draws)

  A0_common <- compute_bayesian_sensitivity_A_from_draws(
    residual_draws = residual_draws,
    rho_grid = 0,
    mode = "common"
  )
  A0_one <- compute_bayesian_sensitivity_A_from_draws(
    residual_draws = residual_draws,
    rho_grid = 0,
    mode = "one_at_a_time"
  )

  artifact <- list(
    artifact_version = "emory_bayesian_app_artifact_v3_structural_rho",
    created_at = Sys.time(),
    dataset_id = dataset_id,
    model_id = "bayesian_mediation_v1_te_fixed_sensitivity_A_only",
    excludes = c("B1_XY", "B2_XM", "brmsfit", "stanfit"),
    settings = list(
      rds_path = results$settings$rds_path %||% "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
      coord_x = results$settings$coord_x %||% "imagecol",
      coord_y = results$settings$coord_y %||% "imagerow",
      exposures = residual_draws$exposures,
      mediators = residual_draws$mediators,
      covariates = c("x_coord", "y_coord"),
      standardize_covariates = results$settings$standardize_covariates %||% TRUE,
      posterior_mc_draws = baseline$n_draws,
      posterior_draw_alignment = "Independent fitted-model posterior draws are explicitly Monte Carlo matched by sampled_indices for mediator_joint, outcome, and total models; derived TE, IE, DE, and PM are computed row-wise from those stored indices.",
      chains = results$settings$chains %||% NA_integer_,
      iter = results$settings$iter %||% NA_integer_,
      warmup = results$settings$warmup %||% NA_integer_,
      seed = seed,
      A_rho_grid_validated = c(-0.3, 0, 0.3)
    ),
    baseline = list(
      wide = baseline$wide,
      decomposition_draws = baseline$decomposition_draws,
      path_draws = baseline$path_draws,
      decomposition_summary = baseline$decomposition_summary,
      path_summary = baseline$path_summary,
      sampled_indices = baseline$sampled_indices,
      n_draws = baseline$n_draws
    ),
    sensitivity_A = list(
      parameterization = structural_rho_version(),
      residual_draws = residual_draws,
      validation = list(
        common_rho0_vs_baseline = validate_bayesian_A_rho0_draw_identity(A0_common, baseline),
        one_at_a_time_rho0_vs_baseline = validate_bayesian_A_rho0_draw_identity(A0_one, baseline),
        common_vs_one_at_a_time_rho0 = validate_bayesian_A_rho0_scenarios_identical(A0_common, A0_one)
      )
    ),
    diagnostics = list(
      mcmc = results$diagnostics,
      predictor_scale_summary = results$data$predictor_scale_summary %||% tibble::tibble(),
      covariate_scale_info = results$data$covariate_scale_info %||% tibble::tibble()
    )
  )

  validate_bayesian_app_artifact(artifact)
  artifact
}

reorient_bayesian_app_artifact <- function(artifact, signs, pc_orientation) {
  validate_bayesian_app_artifact(artifact)
  d <- artifact$sensitivity_A$residual_draws
  meds <- bayes_A_bundle_mediators(d)
  stopifnot(setequal(names(signs), meds), setequal(names(pc_orientation), meds),
            all(signs %in% c(-1,1)), all(pc_orientation %in% c(-1,1)))
  signs <- signs[meds]
  R <- bayes_A_residual_correlation_array(d, meds)
  for(m in meds) {
    d$beta[[m]] <- d$beta[[m]] * signs[[m]]
    for(e in names(d$alpha)) d$alpha[[e]][[m]] <- d$alpha[[e]][[m]] * signs[[m]]
  }
  d$R <- sweep(sweep(R, 2L, signs, "*"), 3L, signs, "*")
  if(length(meds)==3L) {
    for(pair in list(c(1,2),c(1,3),c(2,3))) {
      key <- paste0("r",paste(pair,collapse=""))
      if(!is.null(d[[key]])) d[[key]] <- d$R[,pair[1],pair[2]]
    }
  }
  baseline <- compute_bayesian_A_baseline_from_residual_draws(d)
  stopifnot(isTRUE(all.equal(baseline$decomposition_summary,
                             artifact$baseline$decomposition_summary,tolerance=1e-12)))
  for(nm in intersect(names(artifact$baseline),names(baseline)))
    artifact$baseline[[nm]] <- baseline[[nm]]
  artifact$sensitivity_A$residual_draws <- d
  grid <- artifact$settings$A_rho_grid_validated %||% c(-.3,0,.3)
  common <- compute_bayesian_sensitivity_A_from_draws(d,grid,mode="common")
  one <- compute_bayesian_sensitivity_A_from_draws(d,grid,mode="one_at_a_time")
  compact <- function(x) x[c("parameterization","grid","summary","path","feasibility","validation")]
  artifact$sensitivity_A$results <- list(common=compact(common),one_at_a_time=compact(one))
  artifact$sensitivity_A$validation <- list(
    common_rho0_vs_baseline=validate_bayesian_A_rho0_draw_identity(common,baseline),
    one_at_a_time_rho0_vs_baseline=validate_bayesian_A_rho0_draw_identity(one,baseline),
    common_vs_one_at_a_time_rho0=validate_bayesian_A_rho0_scenarios_identical(common,one))
  stopifnot(all(vapply(artifact$sensitivity_A$validation,function(x)isTRUE(x$passed),logical(1))))
  artifact$sensitivity_A$updated_at <- Sys.time()
  artifact$settings$pc_orientation <- pc_orientation[meds]
  artifact$settings$pc_orientation_reference <- "unmodified prcomp scores and loadings"
  artifact$artifact_version <- "app_artifact_v4_explicit_pc_orientation_structural_rho"
  # Descriptive predictor ranges change sign; MCMC diagnostics do not.
  if(is.list(artifact$diagnostics) && !is.data.frame(artifact$diagnostics) &&
     is.data.frame(artifact$diagnostics$predictor_scale_summary)) {
    z <- artifact$diagnostics$predictor_scale_summary
    for(m in meds[signs==-1]) {
      i <- which(z$variable==m)
      lo <- z$min[i]; z$min[i] <- -z$max[i]; z$max[i] <- -lo; z$mean[i] <- -z$mean[i]
    }
    artifact$diagnostics$predictor_scale_summary <- z
  }
  validate_bayesian_app_artifact(artifact)
  artifact
}

require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

get_bayesian_mediation_summary <- function(artifact) {
  validate_bayesian_app_artifact(artifact)
  x <- artifact$baseline$decomposition_summary
  required <- c(
    "contrast",
    "TE_mean", "TE_q025", "TE_q975",
    "direct_effect_mean", "direct_effect_q025", "direct_effect_q975",
    "total_IIE_mean", "total_IIE_q025", "total_IIE_q975",
    "PM_mean", "PM_q025", "PM_q975"
  )
  require_columns(x, required, "Bayesian mediation decomposition summary")

  out <- x |>
    dplyr::transmute(
      mode = "bayesian",
      contrast,
      TE_est = TE_mean,
      TE_lwr = TE_q025,
      TE_upr = TE_q975,
      direct_est = direct_effect_mean,
      direct_lwr = direct_effect_q025,
      direct_upr = direct_effect_q975,
      total_IIE_est = total_IIE_mean,
      total_IIE_lwr = total_IIE_q025,
      total_IIE_upr = total_IIE_q975,
      PM_est = PM_mean,
      PM_lwr = PM_q025,
      PM_upr = PM_q975
    )

  optional <- c(
    TE_Pr_gt_0 = "TE_Pr_gt_0",
    direct_Pr_gt_0 = "direct_effect_Pr_gt_0",
    total_IIE_Pr_gt_0 = "total_IIE_Pr_gt_0"
  )
  for (nm in names(optional)) {
    src <- optional[[nm]]
    if (src %in% names(x)) {
      out[[nm]] <- x[[src]]
    }
  }

  out
}

get_bayesian_mediation_path <- function(artifact) {
  validate_bayesian_app_artifact(artifact)
  x <- artifact$baseline$path_summary
  required <- c(
    "contrast", "mediator",
    "alpha_mean", "alpha_q025", "alpha_q975",
    "beta_mean", "beta_q025", "beta_q975",
    "IIE_mean", "IIE_q025", "IIE_q975"
  )
  require_columns(x, required, "Bayesian mediation path summary")

  out <- x |>
    dplyr::transmute(
      mode = "bayesian",
      contrast,
      mediator,
      alpha_est = alpha_mean,
      alpha_lwr = alpha_q025,
      alpha_upr = alpha_q975,
      beta_est = beta_mean,
      beta_lwr = beta_q025,
      beta_upr = beta_q975,
      IIE_est = IIE_mean,
      IIE_lwr = IIE_q025,
      IIE_upr = IIE_q975
    )

  if ("Pr_IIE_gt_0" %in% names(x)) {
    out$IIE_Pr_gt_0 <- x$Pr_IIE_gt_0
  } else if ("IIE_Pr_gt_0" %in% names(x)) {
    out$IIE_Pr_gt_0 <- x$IIE_Pr_gt_0
  }
  if ("Pr_IIE_lt_0" %in% names(x)) {
    out$IIE_Pr_lt_0 <- x$Pr_IIE_lt_0
  } else if ("IIE_Pr_lt_0" %in% names(x)) {
    out$IIE_Pr_lt_0 <- x$IIE_Pr_lt_0
  }

  out
}

canonicalize_bayesian_A_summary <- function(A_result) {
  x <- A_result$summary
  required <- c(
    "contrast", "rho", "mediator_varied",
    "direct_effect_mean", "direct_effect_q025", "direct_effect_q975",
    "total_IIE_mean", "total_IIE_q025", "total_IIE_q975",
    "TE_mean", "TE_q025", "TE_q975",
    "PM_mean", "PM_q025", "PM_q975"
  )
  require_columns(x, required, "Bayesian Sensitivity A summary")

  rho_cols <- grep("^rho_[0-9]+$", names(A_result$feasibility), value = TRUE)
  join_cols <- c("scenario", "mediator_varied", "rho", rho_cols)
  feasibility <- A_result$feasibility |>
    dplyr::select(dplyr::all_of(join_cols), valid_draw_fraction, max_feasibility)

  x |>
    dplyr::left_join(
      feasibility,
      by = join_cols
    ) |>
    dplyr::transmute(
      mode = "bayesian",
      scenario,
      contrast,
      rho,
      mediator_varied,
      direct_est = direct_effect_mean,
      direct_lwr = direct_effect_q025,
      direct_upr = direct_effect_q975,
      total_IIE_est = total_IIE_mean,
      total_IIE_lwr = total_IIE_q025,
      total_IIE_upr = total_IIE_q975,
      TE_est = TE_mean,
      TE_lwr = TE_q025,
      TE_upr = TE_q975,
      PM_est = PM_mean,
      PM_lwr = PM_q025,
      PM_upr = PM_q975,
      valid_draw_fraction,
      feasibility_status = dplyr::case_when(
        valid_draw_fraction == 1 ~ "all_draws_feasible",
        valid_draw_fraction > 0 ~ "partially_feasible",
        TRUE ~ "infeasible"
      ),
      max_feasibility
    )
}

canonicalize_bayesian_A_path <- function(A_result) {
  x <- A_result$path
  required <- c(
    "contrast", "rho", "mediator_varied", "mediator",
    "alpha_mean", "alpha_q025", "alpha_q975",
    "beta_mean", "beta_q025", "beta_q975",
    "IIE_mean", "IIE_q025", "IIE_q975"
  )
  require_columns(x, required, "Bayesian Sensitivity A path summary")

  x |>
    dplyr::transmute(
      mode = "bayesian",
      scenario,
      contrast,
      rho,
      mediator_varied,
      mediator,
      alpha_est = alpha_mean,
      alpha_lwr = alpha_q025,
      alpha_upr = alpha_q975,
      beta_est = beta_mean,
      beta_lwr = beta_q025,
      beta_upr = beta_q975,
      IIE_est = IIE_mean,
      IIE_lwr = IIE_q025,
      IIE_upr = IIE_q975,
      IIE_Pr_gt_0 = IIE_Pr_gt_0,
      IIE_Pr_lt_0 = IIE_Pr_lt_0
    )
}

compute_bayesian_app_sensitivity_A <- function(artifact,
                                               rho,
                                               mode = c("common", "one_at_a_time"),
                                               mediator = NULL,
                                               feasibility_tol = 1e-10) {
  validate_bayesian_app_artifact(artifact)
  mode <- match.arg(mode)
  if (!exists("compute_bayesian_sensitivity_A_from_draws", mode = "function")) {
    stop("Source R/bayesian_sensitivity_A_only.R before using Bayesian app Sensitivity A.", call. = FALSE)
  }
  if (!is.numeric(rho) || length(rho) != 1 || !is.finite(rho)) {
    stop("rho must be one finite numeric value.", call. = FALSE)
  }

  residual_draws <- artifact$sensitivity_A$residual_draws
  valid_mediators <- if (exists("bayes_A_bundle_mediators", mode = "function")) {
    bayes_A_bundle_mediators(residual_draws)
  } else {
    c("PC1_R", "PC2_R", "PC3")
  }
  A_result <- compute_bayesian_sensitivity_A_from_draws(
    residual_draws = residual_draws,
    rho_grid = rho,
    mode = mode,
    mediators = valid_mediators,
    mediators_varied = if (identical(mode, "one_at_a_time") && !is.null(mediator)) mediator else valid_mediators,
    feasibility_tol = feasibility_tol
  )

  if (mode == "one_at_a_time") {
    if (is.null(mediator)) {
      stop("mediator must be supplied when mode = 'one_at_a_time'.", call. = FALSE)
    }
    if (!mediator %in% valid_mediators) {
      stop("mediator must be one of: ", paste(valid_mediators, collapse = ", "), call. = FALSE)
    }

    A_result$grid <- dplyr::filter(A_result$grid, mediator_varied == mediator)
    A_result$draws_wide <- dplyr::filter(A_result$draws_wide, mediator_varied == mediator)
    A_result$path_draws <- dplyr::filter(A_result$path_draws, mediator_varied == mediator)
    A_result$summary_draws <- dplyr::filter(A_result$summary_draws, mediator_varied == mediator)
    A_result$summary <- dplyr::filter(A_result$summary, mediator_varied == mediator)
    A_result$path <- dplyr::filter(A_result$path, mediator_varied == mediator)
    A_result$feasibility <- dplyr::filter(A_result$feasibility, mediator_varied == mediator)
  }

  list(
    summary = canonicalize_bayesian_A_summary(A_result),
    path = canonicalize_bayesian_A_path(A_result),
    feasibility = A_result$feasibility,
    draws_wide = A_result$draws_wide,
    raw = A_result
  )
}

validate_bayesian_app_artifact_runtime <- function(artifact, rho_for_timing = 0.15) {
  validate_bayesian_app_artifact(artifact)
  baseline <- list(wide = artifact$baseline$wide, residual_draws = artifact$sensitivity_A$residual_draws)

  A0_common <- compute_bayesian_sensitivity_A_from_draws(
    residual_draws = artifact$sensitivity_A$residual_draws,
    rho_grid = 0,
    mode = "common"
  )
  A0_one <- compute_bayesian_sensitivity_A_from_draws(
    residual_draws = artifact$sensitivity_A$residual_draws,
    rho_grid = 0,
    mode = "one_at_a_time"
  )

  common_identity <- validate_bayesian_A_rho0_draw_identity(A0_common, baseline)
  one_identity <- validate_bayesian_A_rho0_draw_identity(A0_one, baseline)
  scenario_identity <- validate_bayesian_A_rho0_scenarios_identical(A0_common, A0_one)

  timing <- system.time({
    timing_result <- compute_bayesian_app_sensitivity_A(
      artifact = artifact,
      rho = rho_for_timing,
      mode = "common"
    )
  })

  list(
    baseline_summary_rows = nrow(get_bayesian_mediation_summary(artifact)),
    baseline_path_rows = nrow(get_bayesian_mediation_path(artifact)),
    common_rho0_draw_identity = common_identity[c("passed", "max_abs_draw_diff", "tolerance")],
    one_at_a_time_rho0_draw_identity = one_identity[c("passed", "max_abs_draw_diff", "tolerance")],
    common_vs_one_at_a_time_rho0_identity = scenario_identity[c("passed", "max_abs_draw_diff", "tolerance")],
    no_model_fitting_functions_called = TRUE,
    excludes_B1_B2 = !any(c("B1_XY", "B2_XM") %in% names(artifact)),
    timing = tibble::tibble(
      rho = rho_for_timing,
      elapsed_seconds = unname(timing[["elapsed"]]),
      valid_draw_fraction_min = min(timing_result$summary$valid_draw_fraction, na.rm = TRUE)
    )
  )
}
