# Dynamic Bayesian mediation runtime helpers for prepared Shiny data.
#
# This layer fits the validated v1 mediation structure to the currently
# prepared data and returns an app-ready artifact. It deliberately excludes
# B1/B2 sensitivity models and does not persist brmsfit objects in the artifact.

bayesian_dynamic_model_version <- function() {
  "bayesian_mediation_dynamic_v4_original_pca_structural_rho"
}

bayesian_dynamic_default_settings <- function() {
  list(
    chains = 4,
    iter = 2000,
    warmup = 1000,
    cores = min(4, max(1, parallel::detectCores(logical = FALSE) %||% 1)),
    seed = 123,
    backend = NULL,
    adapt_delta = 0.90,
    max_treedepth = 10,
    refresh = 100,
    posterior_mc_draws = NULL,
    standardize_numeric_covariates = TRUE
  )
}

check_bayesian_dynamic_packages <- function(require_brms = TRUE) {
  check_bayesian_v1_packages(require_brms = require_brms)
}

bayesian_dynamic_pc_score <- function(analysis_obj, pc_id, sign = 1) {
  idx <- as.integer(sub("^PC", "", pc_id))
  if (!is.finite(idx) || idx < 1 || idx > ncol(analysis_obj$pca$pca_fit$x)) {
    stop("Selected PC is not available in the PCA fit: ", pc_id, call. = FALSE)
  }

  base <- if (identical(pc_id, "PC1") && "PC1_R" %in% names(analysis_obj$pca_df)) {
    analysis_obj$pca_df$PC1_R
  } else if (identical(pc_id, "PC2") && "PC2_R" %in% names(analysis_obj$pca_df)) {
    analysis_obj$pca_df$PC2_R
  } else if (identical(pc_id, "PC3") && "PC3" %in% names(analysis_obj$pca_df)) {
    analysis_obj$pca_df$PC3
  } else {
    analysis_obj$pca$pca_fit$x[, idx]
  }

  as.numeric(base) * sign
}

bayesian_dynamic_prepare_data <- function(analysis_obj,
                                          selected_pcs,
                                          pc_signs,
                                          standardize_numeric_covariates = TRUE) {
  if (length(selected_pcs) == 0) {
    stop("Select at least one PC mediator before running Bayesian mediation.", call. = FALSE)
  }

  exposures <- analysis_obj$exposures
  if (length(exposures) == 0) {
    stop("No exposure contrasts are available for Bayesian mediation.", call. = FALSE)
  }

  model_data <- analysis_obj$pca_df |>
    dplyr::select(
      Y,
      dplyr::all_of(exposures),
      dplyr::all_of(analysis_obj$covariates)
    ) |>
    as.data.frame()

  mediator_model_names <- paste0("M", seq_along(selected_pcs))
  mediator_map <- tibble::tibble(
    mediator = selected_pcs,
    model_name = mediator_model_names,
    sign = unname(pc_signs[selected_pcs])
  )

  for (i in seq_along(selected_pcs)) {
    model_data[[mediator_model_names[[i]]]] <- bayesian_dynamic_pc_score(
      analysis_obj,
      selected_pcs[[i]],
      sign = mediator_map$sign[[i]]
    )
  }

  covariate_terms <- character(0)
  covariate_scale_info <- tibble::tibble()
  for (cov in analysis_obj$covariates) {
    if (!cov %in% names(model_data)) {
      next
    }
    if (is.numeric(model_data[[cov]])) {
      out <- cov
      if (isTRUE(standardize_numeric_covariates)) {
        out <- paste0(cov, "_std")
        mu <- mean(model_data[[cov]], na.rm = TRUE)
        sig <- stats::sd(model_data[[cov]], na.rm = TRUE)
        if (!is.finite(sig) || sig <= 0) {
          stop("Numeric covariate has zero or invalid SD: ", cov, call. = FALSE)
        }
        model_data[[out]] <- (model_data[[cov]] - mu) / sig
        covariate_scale_info <- dplyr::bind_rows(
          covariate_scale_info,
          tibble::tibble(variable = cov, model_variable = out, center = mu, scale = sig)
        )
      }
      covariate_terms <- c(covariate_terms, out)
    } else {
      model_data[[cov]] <- factor(model_data[[cov]])
      if (nlevels(model_data[[cov]]) < 2) {
        stop("Categorical covariate must have at least two levels: ", cov, call. = FALSE)
      }
      covariate_terms <- c(covariate_terms, cov)
    }
  }

  required <- c("Y", exposures, mediator_model_names, covariate_terms)
  model_data <- tidyr::drop_na(model_data[, unique(required), drop = FALSE])
  if (nrow(model_data) < 5) {
    stop("Too few complete rows are available for Bayesian mediation fitting.", call. = FALSE)
  }

  list(
    data = model_data,
    exposures = exposures,
    exposure_labels = analysis_obj$contrast_labels,
    mediators = selected_pcs,
    mediator_model_names = mediator_model_names,
    mediator_map = mediator_map,
    covariate_terms = covariate_terms,
    covariate_scale_info = covariate_scale_info
  )
}

bayesian_dynamic_formula_part <- function(terms) {
  if (length(terms) == 0) "1" else paste(terms, collapse = " + ")
}

bayesian_dynamic_formulas <- function(prepared) {
  predictor_terms <- c(prepared$exposures, prepared$covariate_terms)
  mediator_rhs <- bayesian_dynamic_formula_part(predictor_terms)
  outcome_rhs <- bayesian_dynamic_formula_part(c(prepared$exposures, prepared$mediator_model_names, prepared$covariate_terms))
  total_rhs <- bayesian_dynamic_formula_part(c(prepared$exposures, prepared$covariate_terms))

  mediator_bfs <- lapply(prepared$mediator_model_names, function(m) {
    brms::bf(stats::as.formula(paste(m, "~", mediator_rhs)))
  })
  mediator_joint <- Reduce(`+`, mediator_bfs)
  if (length(mediator_bfs) > 1) {
    mediator_joint <- mediator_joint + brms::set_rescor(TRUE)
  }

  list(
    mediator_joint = mediator_joint,
    outcome = stats::as.formula(paste("Y ~", outcome_rhs)),
    total = stats::as.formula(paste("Y ~", total_rhs))
  )
}

make_bayesian_dynamic_priors_for_formula <- function(data, formula, resp = NULL) {
  check_bayesian_dynamic_packages(require_brms = TRUE)

  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  y <- stats::model.response(mf)
  y_sd <- stats::sd(y, na.rm = TRUE)
  y_mean <- mean(y, na.rm = TRUE)
  if (!is.finite(y_sd) || y_sd <= 0) {
    stop("Response has zero or invalid SD for formula: ", deparse(formula), call. = FALSE)
  }

  add_resp <- function(args) {
    if (!is.null(resp)) args$resp <- resp
    args
  }

  priors <- c(
    do.call(
      brms::prior_string,
      add_resp(list(
        prior = paste0(
          "student_t(3, ",
          format_bayes_v1_prior_number(y_mean),
          ", ",
          format_bayes_v1_prior_number(10 * y_sd),
          ")"
        ),
        class = "Intercept"
      ))
    ),
    do.call(
      brms::prior_string,
      add_resp(list(
        prior = paste0("exponential(", format_bayes_v1_prior_number(1 / y_sd), ")"),
        class = "sigma"
      ))
    )
  )

  x <- stats::model.matrix(formula, data = data)
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]
  for (coef_name in colnames(x)) {
    x_sd <- stats::sd(x[, coef_name], na.rm = TRUE)
    if (!is.finite(x_sd) || x_sd <= 0) {
      next
    }
    priors <- c(
      priors,
      do.call(
        brms::prior_string,
        add_resp(list(
          prior = paste0("student_t(3, 0, ", format_bayes_v1_prior_number(10 * y_sd / x_sd), ")"),
          class = "b",
          coef = coef_name
        ))
      )
    )
  }

  priors
}

make_bayesian_dynamic_priors <- function(prepared, formulas) {
  mediator_prior_list <- list()
  for (m in prepared$mediator_model_names) {
    f <- stats::as.formula(paste(m, "~", bayesian_dynamic_formula_part(c(prepared$exposures, prepared$covariate_terms))))
    mediator_prior_list <- c(
      mediator_prior_list,
      list(make_bayesian_dynamic_priors_for_formula(
        prepared$data,
        f,
        resp = if (length(prepared$mediator_model_names) > 1) m else NULL
      ))
    )
  }
  mediator_priors <- do.call(c, mediator_prior_list)
  if (length(prepared$mediator_model_names) > 1) {
    mediator_priors <- c(mediator_priors, brms::prior_string("lkj(2)", class = "rescor"))
  }

  list(
    mediator_joint = mediator_priors,
    outcome = make_bayesian_dynamic_priors_for_formula(prepared$data, formulas$outcome),
    total = make_bayesian_dynamic_priors_for_formula(prepared$data, formulas$total)
  )
}

fit_bayesian_dynamic_models <- function(prepared,
                                        formulas = bayesian_dynamic_formulas(prepared),
                                        priors = make_bayesian_dynamic_priors(prepared, formulas),
                                        settings = bayesian_dynamic_default_settings(),
                                        progress_callback = NULL) {
  check_bayesian_dynamic_packages(require_brms = TRUE)

  notify <- function(step, detail = NULL) {
    if (is.function(progress_callback)) {
      progress_callback(step, detail)
    }
  }

  notify("mediator_model", "Fitting joint mediator model")
  mediator_fit <- fit_bayes_v1_single_brm(
    formula = formulas$mediator_joint,
    data = prepared$data,
    prior = priors$mediator_joint,
    chains = settings$chains,
    iter = settings$iter,
    warmup = settings$warmup,
    cores = settings$cores,
    seed = settings$seed + 1,
    backend = settings$backend,
    adapt_delta = settings$adapt_delta,
    max_treedepth = settings$max_treedepth,
    refresh = settings$refresh
  )

  notify("outcome_model", "Fitting outcome model")
  outcome_fit <- fit_bayes_v1_single_brm(
    formula = brms::bf(formulas$outcome),
    data = prepared$data,
    prior = priors$outcome,
    chains = settings$chains,
    iter = settings$iter,
    warmup = settings$warmup,
    cores = settings$cores,
    seed = settings$seed + 2,
    backend = settings$backend,
    adapt_delta = settings$adapt_delta,
    max_treedepth = settings$max_treedepth,
    refresh = settings$refresh
  )

  notify("total_model", "Fitting total-effect model")
  total_fit <- fit_bayes_v1_single_brm(
    formula = brms::bf(formulas$total),
    data = prepared$data,
    prior = priors$total,
    chains = settings$chains,
    iter = settings$iter,
    warmup = settings$warmup,
    cores = settings$cores,
    seed = settings$seed + 3,
    backend = settings$backend,
    adapt_delta = settings$adapt_delta,
    max_treedepth = settings$max_treedepth,
    refresh = settings$refresh
  )

  list(mediator_joint = mediator_fit, outcome = outcome_fit, total = total_fit)
}

extract_bayesian_dynamic_sigma_y <- function(draws_df) {
  draws_df[[find_bayes_v1_draw_column(draws_df, c("sigma", "sigma_Y"), "outcome sigma")]]
}

bayesian_dynamic_rescor_array <- function(draws_df, model_names, mediator_labels) {
  n_draws <- nrow(draws_df)
  k <- length(model_names)
  R <- array(0, dim = c(n_draws, k, k), dimnames = list(NULL, mediator_labels, mediator_labels))
  for (s in seq_len(n_draws)) {
    R[s, , ] <- diag(1, k)
  }
  if (k > 1) {
    for (i in seq_len(k - 1)) {
      for (j in (i + 1):k) {
        rho <- extract_bayes_v1_rescor(draws_df, model_names[[i]], model_names[[j]])
        R[, i, j] <- rho
        R[, j, i] <- rho
      }
    }
  }
  R
}

extract_bayesian_dynamic_residual_draws <- function(fits,
                                                    prepared,
                                                    n_draws = NULL,
                                                    seed = 123,
                                                    sampled_indices = NULL) {
  d_m_all <- posterior::as_draws_df(fits$mediator_joint)
  d_y_all <- posterior::as_draws_df(fits$outcome)
  d_t_all <- posterior::as_draws_df(fits$total)
  available <- c(mediator_joint = nrow(d_m_all), outcome = nrow(d_y_all), total = nrow(d_t_all))

  if (is.null(sampled_indices)) {
    set.seed(seed)
    if (is.null(n_draws)) {
      n_draws <- min(available)
    }
    idx_m <- sample.int(nrow(d_m_all), n_draws, replace = TRUE)
    idx_y <- sample.int(nrow(d_y_all), n_draws, replace = TRUE)
    idx_t <- sample.int(nrow(d_t_all), n_draws, replace = TRUE)
  } else {
    idx_m <- sampled_indices$mediator_joint
    idx_y <- sampled_indices$outcome
    idx_t <- sampled_indices$total
    n_draws <- length(idx_m)
  }

  d_m <- d_m_all[idx_m, , drop = FALSE]
  d_y <- d_y_all[idx_y, , drop = FALSE]
  d_t <- d_t_all[idx_t, , drop = FALSE]
  k <- length(prepared$mediator_model_names)
  mediator_labels <- prepared$mediators

  sig <- if (k == 1) {
    stats::setNames(as.data.frame(list(extract_bayesian_dynamic_sigma_y(d_m))), mediator_labels)
  } else {
    stats::setNames(
      as.data.frame(lapply(prepared$mediator_model_names, function(m) extract_bayes_v1_sigma(d_m, m))),
      mediator_labels
    )
  }
  beta <- stats::setNames(
    as.data.frame(lapply(prepared$mediator_model_names, function(m) extract_bayes_v1_b(d_y, m))),
    mediator_labels
  )

  alpha <- stats::setNames(vector("list", length(prepared$exposures)), prepared$exposures)
  for (exposure in prepared$exposures) {
    alpha[[exposure]] <- stats::setNames(
      as.data.frame(lapply(seq_along(prepared$mediator_model_names), function(i) {
        resp <- if (k > 1) prepared$mediator_model_names[[i]] else NULL
        extract_bayes_v1_b(d_m, exposure, resp = resp)
      })),
      mediator_labels
    )
  }

  outcome_direct_coef <- stats::setNames(
    as.data.frame(lapply(prepared$exposures, function(exposure) extract_bayes_v1_b(d_y, exposure))),
    prepared$exposures
  )
  total_effect <- stats::setNames(
    as.data.frame(lapply(prepared$exposures, function(exposure) extract_bayes_v1_b(d_t, exposure))),
    prepared$exposures
  )

  residual_draws <- list(
    n_draws = n_draws,
    mediators = mediator_labels,
    mediator_model_names = prepared$mediator_model_names,
    exposures = prepared$exposures,
    contrast_labels = prepared$exposure_labels,
    sig = sig,
    R = bayesian_dynamic_rescor_array(d_m, prepared$mediator_model_names, mediator_labels),
    sigma_y = extract_bayesian_dynamic_sigma_y(d_y),
    beta = beta,
    alpha = alpha,
    outcome_direct_coef = outcome_direct_coef,
    total_effect = total_effect,
    sampled_indices = list(mediator_joint = idx_m, outcome = idx_y, total = idx_t),
    available_draws = available
  )

  if (identical(mediator_labels, c("PC1_R", "PC2_R", "PC3")) || identical(mediator_labels, c("PC1", "PC2", "PC3"))) {
    R <- residual_draws$R
    residual_draws$r12 <- R[, 1, 2]
    residual_draws$r13 <- R[, 1, 3]
    residual_draws$r23 <- R[, 2, 3]
  }

  residual_draws
}

summarize_bayesian_dynamic_diagnostics <- function(fits,
                                                   max_treedepth = 10,
                                                   rhat_threshold = 1.01,
                                                   ess_threshold = 400) {
  summarize_bayesian_v1_diagnostics(
    fits = fits,
    max_treedepth = max_treedepth,
    rhat_threshold = rhat_threshold,
    ess_threshold = ess_threshold
  )
}

make_bayesian_dynamic_artifact <- function(analysis_obj,
                                           selected_pcs,
                                           pc_signs,
                                           settings = bayesian_dynamic_default_settings(),
                                           progress_callback = NULL) {
  notify <- function(step, detail = NULL) {
    if (is.function(progress_callback)) {
      progress_callback(step, detail)
    }
  }

  notify("prepare", "Preparing oriented PC mediator data")
  prepared <- bayesian_dynamic_prepare_data(
    analysis_obj,
    selected_pcs = selected_pcs,
    pc_signs = pc_signs,
    standardize_numeric_covariates = settings$standardize_numeric_covariates
  )
  formulas <- bayesian_dynamic_formulas(prepared)
  priors <- make_bayesian_dynamic_priors(prepared, formulas)

  fits <- fit_bayesian_dynamic_models(
    prepared = prepared,
    formulas = formulas,
    priors = priors,
    settings = settings,
    progress_callback = progress_callback
  )

  notify("posterior_processing", "Extracting posterior mediation draws")
  residual_draws <- extract_bayesian_dynamic_residual_draws(
    fits,
    prepared = prepared,
    n_draws = settings$posterior_mc_draws,
    seed = settings$seed + 100
  )
  baseline <- compute_bayesian_A_baseline_from_residual_draws(residual_draws)

  notify("diagnostics", "Computing MCMC diagnostics")
  diagnostics <- summarize_bayesian_dynamic_diagnostics(
    fits,
    max_treedepth = settings$max_treedepth
  )

  artifact <- list(
    artifact_version = "dynamic_app_artifact_v2_structural_rho",
    dataset_id = paste0("prepared_", analysis_obj$spec$data_source),
    model_id = bayesian_dynamic_model_version(),
    settings = c(
      settings,
      list(
        selected_pcs = selected_pcs,
        pc_orientation = stats::setNames(unname(pc_signs[selected_pcs]), selected_pcs),
        pc_orientation_reference = "unmodified prcomp scores and loadings",
        mediator_map = prepared$mediator_map,
        covariate_terms = prepared$covariate_terms,
        covariate_scale_info = prepared$covariate_scale_info
      )
    ),
    baseline = list(
      wide = baseline$wide,
      decomposition_summary = baseline$decomposition_summary,
      path_summary = baseline$path_summary,
      n_draws = baseline$n_draws
    ),
    sensitivity_A = list(
      parameterization = structural_rho_version(),
      residual_draws = residual_draws,
      validation = list(source = "dynamic_fit", rho0_draw_identity_expected = TRUE)
    ),
    diagnostics = diagnostics,
    created_at = Sys.time()
  )

  validate_bayesian_app_artifact(artifact)
  artifact
}

bayesian_dynamic_cache_key <- function(analysis_obj, selected_pcs, pc_signs, settings = bayesian_dynamic_default_settings()) {
  digest::digest(
    list(
      model_version = bayesian_dynamic_model_version(),
      data_source = analysis_obj$spec$data_source,
      metadata = analysis_obj$spec$analysis_data$metadata,
      feature_matrix = analysis_obj$spec$analysis_data$feature_matrix,
      exposures = analysis_obj$exposures,
      covariates = analysis_obj$covariates,
      selected_pcs = selected_pcs,
      pc_signs = unname(pc_signs[selected_pcs]),
      settings = settings[c("chains", "iter", "warmup", "seed", "adapt_delta", "max_treedepth")]
    ),
    algo = "xxhash64"
  )
}
