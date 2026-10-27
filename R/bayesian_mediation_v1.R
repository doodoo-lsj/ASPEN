# Bayesian mediation v1 helper functions.
#
# v1 keeps the current deterministic frequentist PCA preprocessing and treats
# PC1_R, PC2_R, and PC3 as fixed observed mediators. Relative to v0, only the
# mediator model changes from three separate Gaussian regressions to one
# multivariate Gaussian brms model with residual correlations enabled.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

check_bayesian_v1_packages <- function(require_brms = TRUE) {
  required <- c("posterior", "dplyr", "tidyr", "tibble", "purrr")
  if (isTRUE(require_brms)) {
    required <- c("brms", required)
  }
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Bayesian mediation v1 requires installed R packages: ",
      paste(missing, collapse = ", "),
      ". Install them before running scripts/current/run_bayesian_v1.R.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

format_bayes_v1_prior_number <- function(x) {
  formatC(x, digits = 8, format = "fg", flag = "#")
}

summarize_bayes_v1_predictor_scales <- function(data,
                                                variables = c("Y", "PC1_R", "PC2_R", "PC3", "x_coord", "y_coord")) {
  missing <- setdiff(variables, names(data))
  if (length(missing) > 0) {
    stop("Missing variables for scale summary: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  purrr::map_dfr(variables, function(v) {
    x <- data[[v]]
    tibble::tibble(
      variable = v,
      mean = mean(x, na.rm = TRUE),
      sd = stats::sd(x, na.rm = TRUE),
      min = min(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  })
}

standardize_bayes_v1_covariates <- function(data,
                                            covariates = c("x_coord", "y_coord"),
                                            standardized_names = c("x_std", "y_std")) {
  if (length(covariates) != length(standardized_names)) {
    stop("covariates and standardized_names must have the same length.", call. = FALSE)
  }
  missing <- setdiff(covariates, names(data))
  if (length(missing) > 0) {
    stop("Missing covariates for standardization: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  data_std <- data
  scale_rows <- vector("list", length(covariates))

  for (i in seq_along(covariates)) {
    v <- covariates[[i]]
    out <- standardized_names[[i]]
    mu <- mean(data[[v]], na.rm = TRUE)
    sig <- stats::sd(data[[v]], na.rm = TRUE)
    if (!is.finite(sig) || sig <= 0) {
      stop("Covariate has zero or invalid SD: ", v, call. = FALSE)
    }

    data_std[[out]] <- (data[[v]] - mu) / sig
    scale_rows[[i]] <- tibble::tibble(
      variable = v,
      standardized_variable = out,
      center = mu,
      scale = sig
    )
  }

  list(data = data_std, scale_info = dplyr::bind_rows(scale_rows))
}

prepare_bayesian_mediation_v1_data <- function(
    rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
    coord_x = "imagecol",
    coord_y = "imagerow",
    standardize_covariates = TRUE,
    standardize_mediators = FALSE) {
  needed_functions <- c(
    "load_spatial_mediation_rds",
    "prepare_analysis_metadata",
    "fit_pca_mediators",
    "build_pca_score_data",
    "build_full_analysis_data"
  )
  missing_functions <- needed_functions[!vapply(needed_functions, exists, logical(1), mode = "function")]
  if (length(missing_functions) > 0) {
    stop(
      "Missing current pipeline functions: ",
      paste(missing_functions, collapse = ", "),
      ". Source R/pca_mediation_pipeline.R before using Bayesian mediation v1.",
      call. = FALSE
    )
  }
  if (isTRUE(standardize_mediators)) {
    stop(
      "standardize_mediators = TRUE is not implemented in v1 because v1 keeps PC1_R/PC2_R/PC3 on the current mediation scale by default.",
      call. = FALSE
    )
  }

  gex <- load_spatial_mediation_rds(rds_path)
  analysis_df <- prepare_analysis_metadata(gex, coord_x = coord_x, coord_y = coord_y)
  pca <- fit_pca_mediators(gex$M_expr)
  pca_df <- build_pca_score_data(analysis_df, pca$pca_fit)
  dat_full_allpc <- build_full_analysis_data(pca_df)
  predictor_scale_summary <- summarize_bayes_v1_predictor_scales(dat_full_allpc)

  dat_bayes <- dat_full_allpc
  covariate_scale_info <- NULL
  if (isTRUE(standardize_covariates)) {
    standardized <- standardize_bayes_v1_covariates(dat_full_allpc)
    dat_bayes <- standardized$data
    covariate_scale_info <- standardized$scale_info
  }

  # brms multivariate response variable names are easier and safer to extract
  # when they avoid underscores. These aliases are exact copies of the current
  # deterministic PCA mediator scores; output tables map them back to PC1_R,
  # PC2_R, and PC3.
  dat_bayes <- dplyr::mutate(
    dat_bayes,
    PC1R = PC1_R,
    PC2R = PC2_R
  )

  pca_validation <- tibble::tibble(
    check = c("PC1_R_identical", "PC2_R_identical", "PC3_identical"),
    passed = c(
      isTRUE(all.equal(dat_bayes$PC1R, dat_bayes$PC1_R, tolerance = 0)),
      isTRUE(all.equal(dat_bayes$PC2R, dat_bayes$PC2_R, tolerance = 0)),
      isTRUE(all.equal(dat_bayes$PC3, dat_bayes$PC3, tolerance = 0))
    )
  )

  standardized_covariate_validation <- NULL
  if (isTRUE(standardize_covariates)) {
    standardized_covariate_validation <- purrr::map_dfr(c("x_std", "y_std"), function(v) {
      tibble::tibble(
        variable = v,
        mean = mean(dat_bayes[[v]], na.rm = TRUE),
        sd = stats::sd(dat_bayes[[v]], na.rm = TRUE),
        mean_close_to_0 = abs(mean(dat_bayes[[v]], na.rm = TRUE)) < 1e-12,
        sd_close_to_1 = abs(stats::sd(dat_bayes[[v]], na.rm = TRUE) - 1) < 1e-12
      )
    })
  }

  list(
    gex = gex,
    analysis_df = analysis_df,
    pca = pca,
    pca_df = pca_df,
    dat_full_allpc = dat_full_allpc,
    dat_bayes = dat_bayes,
    predictor_scale_summary = predictor_scale_summary,
    covariate_scale_info = covariate_scale_info,
    pca_validation = pca_validation,
    standardized_covariate_validation = standardized_covariate_validation
  )
}

bayesian_mediation_v1_terms <- function(standardize_covariates = TRUE) {
  if (isTRUE(standardize_covariates)) {
    c("x_std", "y_std")
  } else {
    c("x_coord", "y_coord")
  }
}

bayesian_mediation_v1_formulas <- function(standardize_covariates = TRUE) {
  coord_terms <- bayesian_mediation_v1_terms(standardize_covariates)
  coord_part <- paste(coord_terms, collapse = " + ")

  mediator_formula <- brms::bf(stats::as.formula(paste("PC1R ~ X_adjacent + X_far +", coord_part))) +
    brms::bf(stats::as.formula(paste("PC2R ~ X_adjacent + X_far +", coord_part))) +
    brms::bf(stats::as.formula(paste("PC3 ~ X_adjacent + X_far +", coord_part))) +
    brms::set_rescor(TRUE)

  list(
    mediator_joint = mediator_formula,
    outcome = stats::as.formula(paste("Y ~ X_adjacent + X_far + PC1_R + PC2_R + PC3 +", coord_part)),
    total = stats::as.formula(paste("Y ~ X_adjacent + X_far +", coord_part))
  )
}

make_bayes_v1_weak_priors_for_univariate_formula <- function(data, formula) {
  check_bayesian_v1_packages(require_brms = TRUE)

  mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  y <- stats::model.response(mf)
  y_sd <- stats::sd(y, na.rm = TRUE)
  y_mean <- mean(y, na.rm = TRUE)
  if (!is.finite(y_sd) || y_sd <= 0) {
    stop("Response has zero or invalid SD for formula: ", deparse(formula), call. = FALSE)
  }

  x <- stats::model.matrix(formula, data = data)
  x <- x[, colnames(x) != "(Intercept)", drop = FALSE]

  priors <- c(
    brms::prior_string(
      paste0(
        "student_t(3, ",
        format_bayes_v1_prior_number(y_mean),
        ", ",
        format_bayes_v1_prior_number(10 * y_sd),
        ")"
      ),
      class = "Intercept"
    ),
    brms::prior_string(
      paste0("exponential(", format_bayes_v1_prior_number(1 / y_sd), ")"),
      class = "sigma"
    )
  )

  for (coef_name in colnames(x)) {
    x_sd <- stats::sd(x[, coef_name], na.rm = TRUE)
    if (!is.finite(x_sd) || x_sd <= 0) {
      stop("Predictor has zero or invalid SD: ", coef_name, call. = FALSE)
    }
    priors <- c(
      priors,
      brms::prior_string(
        paste0("student_t(3, 0, ", format_bayes_v1_prior_number(10 * y_sd / x_sd), ")"),
        class = "b",
        coef = coef_name
      )
    )
  }

  priors
}

make_bayes_v1_weak_priors_for_response <- function(data,
                                                   response,
                                                   predictors,
                                                   brms_resp_name = response) {
  check_bayesian_v1_packages(require_brms = TRUE)

  y <- data[[response]]
  y_sd <- stats::sd(y, na.rm = TRUE)
  y_mean <- mean(y, na.rm = TRUE)
  if (!is.finite(y_sd) || y_sd <= 0) {
    stop("Response has zero or invalid SD: ", response, call. = FALSE)
  }

  priors <- c(
    brms::prior_string(
      paste0(
        "student_t(3, ",
        format_bayes_v1_prior_number(y_mean),
        ", ",
        format_bayes_v1_prior_number(10 * y_sd),
        ")"
      ),
      class = "Intercept",
      resp = brms_resp_name
    ),
    brms::prior_string(
      paste0("exponential(", format_bayes_v1_prior_number(1 / y_sd), ")"),
      class = "sigma",
      resp = brms_resp_name
    )
  )

  for (coef_name in predictors) {
    x_sd <- stats::sd(data[[coef_name]], na.rm = TRUE)
    if (!is.finite(x_sd) || x_sd <= 0) {
      stop("Predictor has zero or invalid SD: ", coef_name, call. = FALSE)
    }
    priors <- c(
      priors,
      brms::prior_string(
        paste0("student_t(3, 0, ", format_bayes_v1_prior_number(10 * y_sd / x_sd), ")"),
        class = "b",
        coef = coef_name,
        resp = brms_resp_name
      )
    )
  }

  priors
}

make_bayesian_mediation_v1_priors <- function(data,
                                              standardize_covariates = TRUE) {
  coord_terms <- bayesian_mediation_v1_terms(standardize_covariates)
  mediator_predictors <- c("X_adjacent", "X_far", coord_terms)

  mediator_priors <- c(
    make_bayes_v1_weak_priors_for_response(data, "PC1R", mediator_predictors, "PC1R"),
    make_bayes_v1_weak_priors_for_response(data, "PC2R", mediator_predictors, "PC2R"),
    make_bayes_v1_weak_priors_for_response(data, "PC3", mediator_predictors, "PC3"),
    brms::prior_string("lkj(2)", class = "rescor")
  )

  formulas <- bayesian_mediation_v1_formulas(standardize_covariates)
  list(
    mediator_joint = mediator_priors,
    outcome = make_bayes_v1_weak_priors_for_univariate_formula(data, formulas$outcome),
    total = make_bayes_v1_weak_priors_for_univariate_formula(data, formulas$total)
  )
}

fit_bayes_v1_single_brm <- function(formula,
                                    data,
                                    prior,
                                    family = stats::gaussian(),
                                    chains = 4,
                                    iter = 2000,
                                    warmup = 1000,
                                    cores = 4,
                                    seed = 123,
                                    backend = NULL,
                                    adapt_delta = NULL,
                                    max_treedepth = NULL,
                                    refresh = 100) {
  check_bayesian_v1_packages(require_brms = TRUE)

  control <- list()
  if (!is.null(adapt_delta)) {
    control$adapt_delta <- adapt_delta
  }
  if (!is.null(max_treedepth)) {
    control$max_treedepth <- max_treedepth
  }

  brm_args <- list(
    formula = formula,
    data = data,
    family = family,
    prior = prior,
    chains = chains,
    iter = iter,
    warmup = warmup,
    cores = cores,
    seed = seed,
    refresh = refresh,
    control = control
  )

  if (!is.null(backend)) {
    brm_args$backend <- backend
  }

  do.call(brms::brm, brm_args)
}

fit_bayesian_mediation_v1 <- function(data,
                                      formulas = bayesian_mediation_v1_formulas(),
                                      priors = make_bayesian_mediation_v1_priors(data),
                                      chains = 4,
                                      iter = 2000,
                                      warmup = 1000,
                                      cores = 4,
                                      seed = 123,
                                      backend = NULL,
                                      adapt_delta = NULL,
                                      max_treedepth = NULL,
                                      refresh = 100) {
  timing <- list()

  fit_one <- function(model_name, model_seed, formula, prior, family = stats::gaussian()) {
    elapsed <- system.time({
      fit <- fit_bayes_v1_single_brm(
        formula = formula,
        data = data,
        prior = prior,
        family = family,
        chains = chains,
        iter = iter,
        warmup = warmup,
        cores = cores,
        seed = model_seed,
        backend = backend,
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth,
        refresh = refresh
      )
    })
    timing[[model_name]] <<- tibble::tibble(
      component = model_name,
      elapsed_seconds = unname(elapsed[["elapsed"]])
    )
    fit
  }

  fits <- list(
    mediator_joint = fit_one(
      "joint_mediator_model",
      seed + 1,
      formulas$mediator_joint,
      priors$mediator_joint
    ),
    outcome = fit_one(
      "outcome_model",
      seed + 2,
      brms::bf(formulas$outcome),
      priors$outcome
    ),
    total = fit_one(
      "total_model",
      seed + 3,
      brms::bf(formulas$total),
      priors$total
    )
  )

  attr(fits, "timing") <- dplyr::bind_rows(timing)
  fits
}

find_bayes_v1_draw_column <- function(draws_df, candidates, label) {
  existing <- intersect(candidates, names(draws_df))
  if (length(existing) > 0) {
    return(existing[[1]])
  }
  stop(
    "Could not find posterior draw column for ",
    label,
    ". Tried: ",
    paste(candidates, collapse = ", "),
    call. = FALSE
  )
}

extract_bayes_v1_b <- function(draws_df, coef_name, resp = NULL) {
  if (is.null(resp)) {
    candidates <- c(paste0("b_", coef_name))
  } else {
    candidates <- c(
      paste0("b_", resp, "_", coef_name),
      paste0("b_", resp, ".", coef_name)
    )
  }
  draws_df[[find_bayes_v1_draw_column(draws_df, candidates, paste(resp %||% "univariate", coef_name))]]
}

extract_bayes_v1_sigma <- function(draws_df, resp) {
  candidates <- c(paste0("sigma_", resp), paste0("sigma.", resp))
  draws_df[[find_bayes_v1_draw_column(draws_df, candidates, paste("sigma", resp))]]
}

extract_bayes_v1_rescor <- function(draws_df, resp1, resp2) {
  candidates <- c(
    paste0("rescor__", resp1, "__", resp2),
    paste0("rescor__", resp2, "__", resp1),
    paste0("rescor_", resp1, "__", resp2),
    paste0("rescor_", resp2, "__", resp1),
    paste0("rescor_", resp1, "_", resp2),
    paste0("rescor_", resp2, "_", resp1)
  )
  draws_df[[find_bayes_v1_draw_column(draws_df, candidates, paste("rescor", resp1, resp2))]]
}

posterior_summary_bayes_v1 <- function(x) {
  tibble::tibble(
    mean = mean(x, na.rm = TRUE),
    median = stats::median(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    q025 = stats::quantile(x, 0.025, na.rm = TRUE, names = FALSE),
    q975 = stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE),
    Pr_gt_0 = mean(x > 0, na.rm = TRUE),
    Pr_lt_0 = mean(x < 0, na.rm = TRUE)
  )
}

summarize_bayes_v1_tidy_draws <- function(draws,
                                          group_cols = c("contrast", "mediator", "quantity")) {
  draws |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
    dplyr::summarise(
      mean = mean(value, na.rm = TRUE),
      median = stats::median(value, na.rm = TRUE),
      sd = stats::sd(value, na.rm = TRUE),
      q025 = stats::quantile(value, 0.025, na.rm = TRUE),
      q975 = stats::quantile(value, 0.975, na.rm = TRUE),
      Pr_gt_0 = mean(value > 0, na.rm = TRUE),
      Pr_lt_0 = mean(value < 0, na.rm = TRUE),
      .groups = "drop"
    )
}

compute_bayesian_mediation_v1_quantities <- function(fits,
                                                     n_draws = NULL,
                                                     seed = 123) {
  check_bayesian_v1_packages(require_brms = TRUE)
  set.seed(seed)

  draws <- list(
    mediator_joint = posterior::as_draws_df(fits$mediator_joint),
    outcome = posterior::as_draws_df(fits$outcome),
    total = posterior::as_draws_df(fits$total)
  )

  available_draws <- vapply(draws, nrow, integer(1))
  if (is.null(n_draws)) {
    n_draws <- min(available_draws)
  }
  if (n_draws <= 0) {
    stop("n_draws must be positive.", call. = FALSE)
  }

  # The mediator, outcome, and total brms models are fit independently. Thus,
  # identical row numbers are not interpreted as one joint posterior draw across
  # models. Each Monte Carlo row independently samples one draw from each fit.
  # Within the joint mediator model, alpha covariance among PC1_R/PC2_R/PC3 is
  # preserved because all mediator coefficients come from the same sampled row.
  idx <- lapply(available_draws, function(n) sample.int(n, size = n_draws, replace = TRUE))
  names(idx) <- names(available_draws)

  d_m <- draws$mediator_joint[idx$mediator_joint, , drop = FALSE]
  d_y <- draws$outcome[idx$outcome, , drop = FALSE]
  d_t <- draws$total[idx$total, , drop = FALSE]

  alpha_A_PC1_R <- extract_bayes_v1_b(d_m, "X_adjacent", "PC1R")
  alpha_F_PC1_R <- extract_bayes_v1_b(d_m, "X_far", "PC1R")
  alpha_A_PC2_R <- extract_bayes_v1_b(d_m, "X_adjacent", "PC2R")
  alpha_F_PC2_R <- extract_bayes_v1_b(d_m, "X_far", "PC2R")
  alpha_A_PC3 <- extract_bayes_v1_b(d_m, "X_adjacent", "PC3")
  alpha_F_PC3 <- extract_bayes_v1_b(d_m, "X_far", "PC3")

  beta_PC1_R <- extract_bayes_v1_b(d_y, "PC1_R")
  beta_PC2_R <- extract_bayes_v1_b(d_y, "PC2_R")
  beta_PC3 <- extract_bayes_v1_b(d_y, "PC3")

  IIE_PC1_A <- alpha_A_PC1_R * beta_PC1_R
  IIE_PC2_A <- alpha_A_PC2_R * beta_PC2_R
  IIE_PC3_A <- alpha_A_PC3 * beta_PC3
  IIE_PC1_F <- alpha_F_PC1_R * beta_PC1_R
  IIE_PC2_F <- alpha_F_PC2_R * beta_PC2_R
  IIE_PC3_F <- alpha_F_PC3 * beta_PC3

  total_IIE_A <- IIE_PC1_A + IIE_PC2_A + IIE_PC3_A
  total_IIE_F <- IIE_PC1_F + IIE_PC2_F + IIE_PC3_F

  outcome_direct_A <- extract_bayes_v1_b(d_y, "X_adjacent")
  outcome_direct_F <- extract_bayes_v1_b(d_y, "X_far")
  TE_A <- extract_bayes_v1_b(d_t, "X_adjacent")
  TE_F <- extract_bayes_v1_b(d_t, "X_far")
  direct_A <- TE_A - total_IIE_A
  direct_F <- TE_F - total_IIE_F
  PM_A <- total_IIE_A / TE_A
  PM_F <- total_IIE_F / TE_F

  wide <- tibble::tibble(
    draw = seq_len(n_draws),
    alpha_A_PC1_R = alpha_A_PC1_R,
    alpha_F_PC1_R = alpha_F_PC1_R,
    alpha_A_PC2_R = alpha_A_PC2_R,
    alpha_F_PC2_R = alpha_F_PC2_R,
    alpha_A_PC3 = alpha_A_PC3,
    alpha_F_PC3 = alpha_F_PC3,
    beta_PC1_R = beta_PC1_R,
    beta_PC2_R = beta_PC2_R,
    beta_PC3 = beta_PC3,
    IIE_PC1_A = IIE_PC1_A,
    IIE_PC2_A = IIE_PC2_A,
    IIE_PC3_A = IIE_PC3_A,
    IIE_PC1_F = IIE_PC1_F,
    IIE_PC2_F = IIE_PC2_F,
    IIE_PC3_F = IIE_PC3_F,
    total_IIE_A = total_IIE_A,
    total_IIE_F = total_IIE_F,
    direct_A = direct_A,
    direct_F = direct_F,
    outcome_direct_A = outcome_direct_A,
    outcome_direct_F = outcome_direct_F,
    TE_A = TE_A,
    TE_F = TE_F,
    PM_A = PM_A,
    PM_F = PM_F
  )

  path_long <- dplyr::bind_rows(
    tibble::tibble(draw = wide$draw, contrast = "adjacent vs inside", exposure = "X_adjacent", mediator = "PC1_R", alpha = wide$alpha_A_PC1_R, beta = wide$beta_PC1_R, value = wide$IIE_PC1_A),
    tibble::tibble(draw = wide$draw, contrast = "adjacent vs inside", exposure = "X_adjacent", mediator = "PC2_R", alpha = wide$alpha_A_PC2_R, beta = wide$beta_PC2_R, value = wide$IIE_PC2_A),
    tibble::tibble(draw = wide$draw, contrast = "adjacent vs inside", exposure = "X_adjacent", mediator = "PC3", alpha = wide$alpha_A_PC3, beta = wide$beta_PC3, value = wide$IIE_PC3_A),
    tibble::tibble(draw = wide$draw, contrast = "far vs inside", exposure = "X_far", mediator = "PC1_R", alpha = wide$alpha_F_PC1_R, beta = wide$beta_PC1_R, value = wide$IIE_PC1_F),
    tibble::tibble(draw = wide$draw, contrast = "far vs inside", exposure = "X_far", mediator = "PC2_R", alpha = wide$alpha_F_PC2_R, beta = wide$beta_PC2_R, value = wide$IIE_PC2_F),
    tibble::tibble(draw = wide$draw, contrast = "far vs inside", exposure = "X_far", mediator = "PC3", alpha = wide$alpha_F_PC3, beta = wide$beta_PC3, value = wide$IIE_PC3_F)
  ) |>
    dplyr::mutate(quantity = "IIE", .before = value)

  decomposition_long <- dplyr::bind_rows(
    tibble::tibble(draw = wide$draw, contrast = "adjacent vs inside", exposure = "X_adjacent", quantity = "TE", value = wide$TE_A),
    tibble::tibble(draw = wide$draw, contrast = "adjacent vs inside", exposure = "X_adjacent", quantity = "direct_effect", value = wide$direct_A),
    tibble::tibble(draw = wide$draw, contrast = "adjacent vs inside", exposure = "X_adjacent", quantity = "outcome_direct_coef", value = wide$outcome_direct_A),
    tibble::tibble(draw = wide$draw, contrast = "adjacent vs inside", exposure = "X_adjacent", quantity = "total_IIE", value = wide$total_IIE_A),
    tibble::tibble(draw = wide$draw, contrast = "adjacent vs inside", exposure = "X_adjacent", quantity = "PM", value = wide$PM_A),
    tibble::tibble(draw = wide$draw, contrast = "far vs inside", exposure = "X_far", quantity = "TE", value = wide$TE_F),
    tibble::tibble(draw = wide$draw, contrast = "far vs inside", exposure = "X_far", quantity = "direct_effect", value = wide$direct_F),
    tibble::tibble(draw = wide$draw, contrast = "far vs inside", exposure = "X_far", quantity = "outcome_direct_coef", value = wide$outcome_direct_F),
    tibble::tibble(draw = wide$draw, contrast = "far vs inside", exposure = "X_far", quantity = "total_IIE", value = wide$total_IIE_F),
    tibble::tibble(draw = wide$draw, contrast = "far vs inside", exposure = "X_far", quantity = "PM", value = wide$PM_F)
  )

  path_summary <- path_long |>
    dplyr::group_by(contrast, exposure, mediator) |>
    dplyr::summarise(
      alpha_mean = mean(alpha, na.rm = TRUE),
      alpha_q025 = stats::quantile(alpha, 0.025, na.rm = TRUE),
      alpha_q975 = stats::quantile(alpha, 0.975, na.rm = TRUE),
      beta_mean = mean(beta, na.rm = TRUE),
      beta_q025 = stats::quantile(beta, 0.025, na.rm = TRUE),
      beta_q975 = stats::quantile(beta, 0.975, na.rm = TRUE),
      IIE_mean = mean(value, na.rm = TRUE),
      IIE_median = stats::median(value, na.rm = TRUE),
      IIE_q025 = stats::quantile(value, 0.025, na.rm = TRUE),
      IIE_q975 = stats::quantile(value, 0.975, na.rm = TRUE),
      Pr_IIE_gt_0 = mean(value > 0, na.rm = TRUE),
      Pr_IIE_lt_0 = mean(value < 0, na.rm = TRUE),
      .groups = "drop"
    )

  decomposition_summary <- decomposition_long |>
    summarize_bayes_v1_tidy_draws(group_cols = c("contrast", "exposure", "quantity")) |>
    tidyr::pivot_wider(
      names_from = quantity,
      values_from = c(mean, median, sd, q025, q975, Pr_gt_0, Pr_lt_0),
      names_glue = "{quantity}_{.value}"
    )

  pm_diagnostics <- dplyr::bind_rows(
    tibble::tibble(contrast = "adjacent vs inside", exposure = "X_adjacent", TE = wide$TE_A, PM = wide$PM_A),
    tibble::tibble(contrast = "far vs inside", exposure = "X_far", TE = wide$TE_F, PM = wide$PM_F)
  ) |>
    dplyr::group_by(contrast, exposure) |>
    dplyr::summarise(
      min_abs_TE = min(abs(TE), na.rm = TRUE),
      n_abs_TE_lt_0_01 = sum(abs(TE) < 0.01, na.rm = TRUE),
      n_abs_TE_lt_0_001 = sum(abs(TE) < 0.001, na.rm = TRUE),
      max_abs_PM = max(abs(PM), na.rm = TRUE),
      .groups = "drop"
    )

  list(
    draws_wide = wide,
    draws_tidy = dplyr::bind_rows(
      path_long |> dplyr::select(draw, contrast, exposure, mediator, quantity, value),
      decomposition_long |> dplyr::mutate(mediator = NA_character_) |> dplyr::select(draw, contrast, exposure, mediator, quantity, value)
    ),
    path_draws = path_long,
    decomposition_draws = decomposition_long,
    path_summary = path_summary,
    decomposition_summary = decomposition_summary,
    pm_diagnostics = pm_diagnostics,
    sampled_indices = idx,
    available_draws = available_draws,
    n_draws = n_draws
  )
}

extract_bayesian_v1_app_residual_draws <- function(fits,
                                                   n_draws = NULL,
                                                   seed = 123,
                                                   sampled_indices = NULL) {
  check_bayesian_v1_packages(require_brms = TRUE)

  draws <- list(
    mediator_joint = posterior::as_draws_df(fits$mediator_joint),
    outcome = posterior::as_draws_df(fits$outcome),
    total = posterior::as_draws_df(fits$total)
  )
  available <- vapply(draws, nrow, integer(1))

  if (is.null(sampled_indices)) {
    set.seed(seed)
    if (is.null(n_draws)) {
      n_draws <- min(available)
    }
    sampled_indices <- lapply(available, function(n) sample.int(n, size = n_draws, replace = TRUE))
  } else {
    required <- c("mediator_joint", "outcome", "total")
    missing <- setdiff(required, names(sampled_indices))
    if (length(missing) > 0) {
      stop("sampled_indices is missing: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    n_draws <- length(sampled_indices$mediator_joint)
    if (length(sampled_indices$outcome) != n_draws || length(sampled_indices$total) != n_draws) {
      stop("All sampled index vectors must have the same length.", call. = FALSE)
    }
  }

  d_m <- draws$mediator_joint[sampled_indices$mediator_joint, , drop = FALSE]
  d_y <- draws$outcome[sampled_indices$outcome, , drop = FALSE]
  d_t <- draws$total[sampled_indices$total, , drop = FALSE]

  mediators <- c("PC1_R", "PC2_R", "PC3")
  exposures <- c("X_adjacent", "X_far")
  contrast_labels <- c(X_adjacent = "adjacent vs inside", X_far = "far vs inside")

  sig <- data.frame(
    PC1_R = extract_bayes_v1_sigma(d_m, "PC1R"),
    PC2_R = extract_bayes_v1_sigma(d_m, "PC2R"),
    PC3 = extract_bayes_v1_sigma(d_m, "PC3"),
    check.names = FALSE
  )
  R <- array(0, dim = c(n_draws, 3, 3), dimnames = list(NULL, mediators, mediators))
  for (s in seq_len(n_draws)) {
    R[s, , ] <- diag(1, 3)
  }
  R[, 1, 2] <- R[, 2, 1] <- extract_bayes_v1_rescor(d_m, "PC1R", "PC2R")
  R[, 1, 3] <- R[, 3, 1] <- extract_bayes_v1_rescor(d_m, "PC1R", "PC3")
  R[, 2, 3] <- R[, 3, 2] <- extract_bayes_v1_rescor(d_m, "PC2R", "PC3")

  beta <- data.frame(
    PC1_R = extract_bayes_v1_b(d_y, "PC1_R"),
    PC2_R = extract_bayes_v1_b(d_y, "PC2_R"),
    PC3 = extract_bayes_v1_b(d_y, "PC3"),
    check.names = FALSE
  )
  alpha <- list(
    X_adjacent = data.frame(
      PC1_R = extract_bayes_v1_b(d_m, "X_adjacent", "PC1R"),
      PC2_R = extract_bayes_v1_b(d_m, "X_adjacent", "PC2R"),
      PC3 = extract_bayes_v1_b(d_m, "X_adjacent", "PC3"),
      check.names = FALSE
    ),
    X_far = data.frame(
      PC1_R = extract_bayes_v1_b(d_m, "X_far", "PC1R"),
      PC2_R = extract_bayes_v1_b(d_m, "X_far", "PC2R"),
      PC3 = extract_bayes_v1_b(d_m, "X_far", "PC3"),
      check.names = FALSE
    )
  )
  outcome_direct_coef <- data.frame(
    X_adjacent = extract_bayes_v1_b(d_y, "X_adjacent"),
    X_far = extract_bayes_v1_b(d_y, "X_far"),
    check.names = FALSE
  )
  total_effect <- data.frame(
    X_adjacent = extract_bayes_v1_b(d_t, "X_adjacent"),
    X_far = extract_bayes_v1_b(d_t, "X_far"),
    check.names = FALSE
  )

  list(
    n_draws = n_draws,
    mediators = mediators,
    exposures = exposures,
    contrast_labels = contrast_labels,
    sig = sig,
    R = R,
    r12 = R[, 1, 2],
    r13 = R[, 1, 3],
    r23 = R[, 2, 3],
    sigma_y = d_y[[find_bayes_v1_draw_column(d_y, c("sigma", "sigma_Y"), "outcome sigma")]],
    beta = beta,
    alpha = alpha,
    outcome_direct_coef = outcome_direct_coef,
    total_effect = total_effect,
    sampled_indices = sampled_indices,
    available_draws = available
  )
}

summarize_mediator_residual_dependence_v1 <- function(fit) {
  check_bayesian_v1_packages(require_brms = TRUE)
  draws <- posterior::as_draws_df(fit)
  responses <- c("PC1R", "PC2R", "PC3")
  labels <- c(PC1R = "PC1_R", PC2R = "PC2_R", PC3 = "PC3")

  sigma <- stats::setNames(
    lapply(responses, function(resp) extract_bayes_v1_sigma(draws, resp)),
    responses
  )

  rescor_pairs <- tibble::tribble(
    ~resp1, ~resp2,
    "PC1R", "PC2R",
    "PC1R", "PC3",
    "PC2R", "PC3"
  )

  cor_draws <- purrr::pmap_dfr(rescor_pairs, function(resp1, resp2) {
    rho <- extract_bayes_v1_rescor(draws, resp1, resp2)
    tibble::tibble(
      draw = seq_along(rho),
      mediator_1 = labels[[resp1]],
      mediator_2 = labels[[resp2]],
      residual_correlation = rho,
      residual_covariance = rho * sigma[[resp1]] * sigma[[resp2]]
    )
  })

  cor_summary <- cor_draws |>
    dplyr::group_by(mediator_1, mediator_2) |>
    dplyr::summarise(
      correlation_mean = mean(residual_correlation, na.rm = TRUE),
      correlation_median = stats::median(residual_correlation, na.rm = TRUE),
      correlation_q025 = stats::quantile(residual_correlation, 0.025, na.rm = TRUE),
      correlation_q975 = stats::quantile(residual_correlation, 0.975, na.rm = TRUE),
      covariance_mean = mean(residual_covariance, na.rm = TRUE),
      covariance_q025 = stats::quantile(residual_covariance, 0.025, na.rm = TRUE),
      covariance_q975 = stats::quantile(residual_covariance, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  mean_cov <- diag(vapply(responses, function(resp) mean(sigma[[resp]]^2, na.rm = TRUE), numeric(1)))
  mean_cor <- diag(1, length(responses))
  dimnames(mean_cov) <- list(labels[responses], labels[responses])
  dimnames(mean_cor) <- list(labels[responses], labels[responses])

  for (i in seq_len(nrow(cor_summary))) {
    m1 <- cor_summary$mediator_1[[i]]
    m2 <- cor_summary$mediator_2[[i]]
    mean_cov[m1, m2] <- cor_summary$covariance_mean[[i]]
    mean_cov[m2, m1] <- cor_summary$covariance_mean[[i]]
    mean_cor[m1, m2] <- cor_summary$correlation_mean[[i]]
    mean_cor[m2, m1] <- cor_summary$correlation_mean[[i]]
  }

  list(
    posterior_draws = cor_draws,
    correlation_summary = cor_summary,
    covariance_mean_matrix = mean_cov,
    correlation_mean_matrix = mean_cor,
    residual_correlation_is_estimated = nrow(cor_summary) == 3
  )
}

compute_bayes_v1_frequentist_estimates <- function(data) {
  if (!exists("decompose_linear_multix", mode = "function")) {
    stop("decompose_linear_multix() is missing. Source R/pca_mediation_pipeline.R first.", call. = FALSE)
  }

  fit <- decompose_linear_multix(
    data = data,
    exposures = c("X_adjacent", "X_far"),
    outcome = "Y",
    mediators = c("PC1_R", "PC2_R", "PC3"),
    covariates = c("x_coord", "y_coord")
  )

  list(
    fit = fit,
    summary = fit$summary |>
      dplyr::mutate(
        contrast = dplyr::case_when(
          exposure == "X_adjacent" ~ "adjacent vs inside",
          exposure == "X_far" ~ "far vs inside",
          TRUE ~ exposure
        )
      ),
    path = fit$path |>
      dplyr::mutate(
        contrast = dplyr::case_when(
          exposure == "X_adjacent" ~ "adjacent vs inside",
          exposure == "X_far" ~ "far vs inside",
          TRUE ~ exposure
        )
      )
  )
}

compare_frequentist_bayesian_v1 <- function(frequentist,
                                            path_summary,
                                            decomposition_summary) {
  decomp_long <- decomposition_summary |>
    dplyr::select(
      contrast,
      bayesian_TE_mean = TE_mean,
      bayesian_direct_mean = direct_effect_mean,
      bayesian_total_IIE_mean = total_IIE_mean
    )

  decomposition <- frequentist$summary |>
    dplyr::select(
      contrast,
      frequentist_TE = TE,
      frequentist_direct = NDE,
      frequentist_NIE = NIE
    ) |>
    dplyr::left_join(decomp_long, by = "contrast")

  path <- frequentist$path |>
    dplyr::select(
      contrast,
      mediator,
      frequentist_alpha = alpha_X_to_M,
      frequentist_beta = beta_M_to_Y,
      frequentist_indirect = indirect_component
    ) |>
    dplyr::left_join(
      path_summary |>
        dplyr::select(
          contrast,
          mediator,
          bayesian_alpha_mean = alpha_mean,
          bayesian_beta_mean = beta_mean,
          bayesian_IIE_mean = IIE_mean
        ),
      by = c("contrast", "mediator")
    )

  list(decomposition = decomposition, path = path)
}

summarize_bayesian_v1_diagnostics <- function(fits,
                                              max_treedepth = 10,
                                              rhat_threshold = 1.01,
                                              ess_threshold = 400) {
  check_bayesian_v1_packages(require_brms = TRUE)

  purrr::imap_dfr(fits, function(fit, model) {
    draw_diag <- posterior::summarise_draws(
      posterior::as_draws(fit),
      "rhat",
      "ess_bulk",
      "ess_tail"
    )

    nuts <- tryCatch(brms::nuts_params(fit), error = function(e) NULL)
    divergent_transitions <- NA_integer_
    max_treedepth_hits <- NA_integer_
    if (!is.null(nuts)) {
      divergent_transitions <- sum(nuts$Parameter == "divergent__" & nuts$Value == 1, na.rm = TRUE)
      max_treedepth_hits <- sum(nuts$Parameter == "treedepth__" & nuts$Value >= max_treedepth, na.rm = TRUE)
    }

    max_Rhat <- max(draw_diag$rhat, na.rm = TRUE)
    min_bulk_ESS <- min(draw_diag$ess_bulk, na.rm = TRUE)
    min_tail_ESS <- min(draw_diag$ess_tail, na.rm = TRUE)

    tibble::tibble(
      model = model,
      max_Rhat = max_Rhat,
      min_bulk_ESS = min_bulk_ESS,
      min_tail_ESS = min_tail_ESS,
      divergent_transitions = divergent_transitions,
      max_treedepth_hits = max_treedepth_hits,
      rhat_pass = is.finite(max_Rhat) && max_Rhat <= rhat_threshold,
      bulk_ESS_pass = is.finite(min_bulk_ESS) && min_bulk_ESS >= ess_threshold,
      tail_ESS_pass = is.finite(min_tail_ESS) && min_tail_ESS >= ess_threshold,
      divergence_pass = !is.na(divergent_transitions) && divergent_transitions == 0,
      treedepth_pass = !is.na(max_treedepth_hits) && max_treedepth_hits == 0,
      overall_pass = rhat_pass && bulk_ESS_pass && tail_ESS_pass && divergence_pass && treedepth_pass
    )
  })
}
