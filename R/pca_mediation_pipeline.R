load_spatial_mediation_rds <- function(rds_path) {
  gex <- readRDS(rds_path)
  required <- c("X_cat", "Y", "meta", "M_expr")
  missing <- setdiff(required, names(gex))
  if (length(missing) > 0) {
    stop("RDS is missing required elements: ", paste(missing, collapse = ", "))
  }
  gex
}

prepare_analysis_metadata <- function(gex, coord_x = "imagecol", coord_y = "imagerow") {
  meta <- as.data.frame(gex$meta)
  missing_coords <- setdiff(c(coord_x, coord_y), names(meta))
  if (length(missing_coords) > 0) {
    stop("Metadata is missing coordinate columns: ", paste(missing_coords, collapse = ", "))
  }

  x_values <- if ("X_cat" %in% names(meta)) meta$X_cat else gex$X_cat

  dplyr::mutate(
    meta,
    X_cat = factor(x_values, levels = c("inside", "adjacent", "far")),
    Y = as.numeric(gex$Y),
    x_coord = .data[[coord_x]],
    y_coord = .data[[coord_y]]
  )
}

randomized_prcomp <- function(x,
                              rank = 6L,
                              oversample = 30L,
                              power_iterations = 9L,
                              seed = 20261001L) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  rank <- min(as.integer(rank), n - 1L, p)
  if (!is.finite(rank) || rank < 1L) {
    stop("PCA requires at least two observations and one feature.", call. = FALSE)
  }

  sketch_rank <- min(rank + as.integer(oversample), n, p)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  omega <- matrix(stats::rnorm(p * sketch_rank), nrow = p, ncol = sketch_rank)
  q <- qr.Q(qr(x %*% omega), complete = FALSE)
  for (iteration in seq_len(as.integer(power_iterations))) {
    q_feature <- qr.Q(qr(crossprod(x, q)), complete = FALSE)
    q <- qr.Q(qr(x %*% q_feature), complete = FALSE)
  }

  small_svd <- base::svd(crossprod(q, x), nu = rank, nv = rank)
  rotation <- small_svd$v[, seq_len(rank), drop = FALSE]
  scores <- (q %*% small_svd$u[, seq_len(rank), drop = FALSE]) %*%
    diag(small_svd$d[seq_len(rank)], nrow = rank)

  # Fix the otherwise arbitrary SVD signs so repeated uploads are reproducible.
  anchors <- apply(rotation, 2, function(column) column[[which.max(abs(column))]])
  signs <- ifelse(anchors < 0, -1, 1)
  rotation <- sweep(rotation, 2, signs, "*")
  scores <- sweep(scores, 2, signs, "*")
  pc_names <- paste0("PC", seq_len(rank))
  colnames(rotation) <- pc_names
  rownames(rotation) <- colnames(x)
  colnames(scores) <- pc_names
  rownames(scores) <- rownames(x)

  fit <- list(
    sdev = small_svd$d[seq_len(rank)] / sqrt(max(1, n - 1L)),
    rotation = rotation,
    center = FALSE,
    scale = FALSE,
    x = scores
  )
  class(fit) <- "prcomp"
  fit
}

fit_pca_mediators <- function(M_expr,
                              method = c("exact", "randomized"),
                              rank = 6L) {
  method <- match.arg(method)
  M_scaled <- scale(M_expr)
  pca_fit <- if (identical(method, "randomized")) {
    randomized_prcomp(M_scaled, rank = rank)
  } else {
    stats::prcomp(M_scaled, center = FALSE, scale. = FALSE)
  }
  total_variance <- sum(M_scaled^2) / max(1, nrow(M_scaled) - 1L)

  list(
    M_scaled = M_scaled,
    pca_fit = pca_fit,
    var_explained = pca_fit$sdev^2 / total_variance,
    method = method
  )
}

build_pca_score_data <- function(analysis_df, pca_fit) {
  dplyr::mutate(
    analysis_df,
    PC1_original = pca_fit$x[, 1],
    PC2_original = pca_fit$x[, 2],
    PC3_original = pca_fit$x[, 3],
    PC1_R = -PC1_original,
    PC2_R = -PC2_original,
    PC3 = PC3_original
  )
}

build_scree_data <- function(var_explained) {
  tibble::tibble(
    PC = seq_along(var_explained),
    var_explained = var_explained,
    cum_var_explained = cumsum(var_explained)
  )
}

build_loading_data <- function(pca_fit) {
  stats::setNames(
    as.data.frame(pca_fit$rotation[, 1:3]),
    c("PC1_original", "PC2_original", "PC3_original")
  ) |>
    tibble::rownames_to_column("gene") |>
    dplyr::mutate(
      PC1_R = -PC1_original,
      PC2_R = -PC2_original,
      PC3 = PC3_original
    )
}

extract_top_loading_genes <- function(loading_df,
                                      pc1_n = 150,
                                      pc2_n = 100,
                                      pc3_n = 100) {
  list(
    pc1_pos = loading_df |>
      dplyr::arrange(dplyr::desc(PC1_R)) |>
      dplyr::slice_head(n = pc1_n) |>
      dplyr::pull(gene),
    pc1_neg = loading_df |>
      dplyr::arrange(PC1_R) |>
      dplyr::slice_head(n = pc1_n) |>
      dplyr::pull(gene),
    pc2_pos = loading_df |>
      dplyr::arrange(dplyr::desc(PC2_R)) |>
      dplyr::slice_head(n = pc2_n) |>
      dplyr::pull(gene),
    pc2_neg = loading_df |>
      dplyr::arrange(PC2_R) |>
      dplyr::slice_head(n = pc2_n) |>
      dplyr::pull(gene),
    pc3_pos = loading_df |>
      dplyr::arrange(dplyr::desc(PC3)) |>
      dplyr::slice_head(n = pc3_n) |>
      dplyr::pull(gene),
    pc3_neg = loading_df |>
      dplyr::arrange(PC3) |>
      dplyr::slice_head(n = pc3_n) |>
      dplyr::pull(gene)
  )
}

build_spatial_plot_inputs <- function(pca_df, buffer_scale = 0.58, flip_y = FALSE) {
  plot_df <- dplyr::mutate(
    pca_df,
    x_plot = x_coord,
    y_plot = if (isTRUE(flip_y)) -y_coord else y_coord
  )

  xy_mat <- as.matrix(plot_df[, c("x_plot", "y_plot")])
  nn_result <- RANN::nn2(data = xy_mat, query = xy_mat, k = 2)
  nearest_distance <- stats::median(nn_result$nn.dists[, 2])

  spot_sf <- sf::st_as_sf(
    plot_df,
    coords = c("x_plot", "y_plot"),
    remove = FALSE,
    crs = NA
  )

  buffer_radius <- nearest_distance * buffer_scale
  spot_buffer_sf <- sf::st_buffer(spot_sf, dist = buffer_radius)

  category_polygon_sf <- spot_buffer_sf |>
    dplyr::group_by(X_cat) |>
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop") |>
    sf::st_make_valid()

  list(
    plot_df = plot_df,
    nearest_distance = nearest_distance,
    buffer_radius = buffer_radius,
    flip_y = isTRUE(flip_y),
    spot_sf = spot_sf,
    spot_buffer_sf = spot_buffer_sf,
    category_polygon_sf = category_polygon_sf
  )
}

plot_pc_spatial <- function(data,
                            pc_var,
                            variance_explained,
                            pc_label,
                            boundary_sf,
                            flip_y = FALSE,
                            padding_fraction = 0.06) {
  x_breaks <- pretty(data$x_plot, n = 5)
  y_breaks <- pretty(data$y_plot, n = 5)
  x_range <- range(data$x_plot, na.rm = TRUE)
  y_range <- range(data$y_plot, na.rm = TRUE)
  x_pad <- diff(x_range) * padding_fraction
  y_pad <- diff(y_range) * padding_fraction
  if (!is.finite(x_pad) || x_pad == 0) x_pad <- 1
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 1

  ggplot2::ggplot() +
    ggplot2::geom_point(
      data = data,
      ggplot2::aes(x = x_plot, y = y_plot, color = .data[[pc_var]]),
      size = 2,
      alpha = 0.88
    ) +
    ggplot2::scale_color_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      name = pc_label
    ) +
    ggplot2::scale_x_continuous(breaks = x_breaks) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = if (isTRUE(flip_y)) function(z) -z else ggplot2::waiver()
    ) +
    ggplot2::coord_equal(
      xlim = c(x_range[[1]] - x_pad, x_range[[2]] + x_pad),
      ylim = c(y_range[[1]] - y_pad, y_range[[2]] + y_pad),
      expand = FALSE
    ) +
    ggplot2::labs(
      title = paste0(
        pc_label,
        " spatial score map\n(",
        sprintf("%.2f", 100 * variance_explained),
        "% variance explained)"
      ),
      x = "x coordinate",
      y = "y coordinate"
    ) +
    ggplot2::theme_bw(base_size = 15) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid = ggplot2::element_blank()
    )
}

build_spatial_pc_plots <- function(plot_df, var_explained, boundary_sf, flip_y = FALSE) {
  list(
    p_pc1 = plot_pc_spatial(plot_df, "PC1_R", var_explained[1], "PC1_R", boundary_sf, flip_y),
    p_pc2 = plot_pc_spatial(plot_df, "PC2_R", var_explained[2], "PC2_R", boundary_sf, flip_y),
    p_pc3 = plot_pc_spatial(plot_df, "PC3", var_explained[3], "PC3", boundary_sf, flip_y)
  )
}

decompose_linear <- function(data, exposure, outcome, mediators, covariates = NULL) {
  cov_part <- if (!is.null(covariates) && length(covariates) > 0) {
    paste(covariates, collapse = " + ")
  } else {
    NULL
  }

  te_formula <- stats::as.formula(
    paste(
      outcome, "~", exposure,
      if (!is.null(cov_part)) paste("+", cov_part) else ""
    )
  )
  fit_te <- stats::lm(te_formula, data = data)
  TE <- stats::coef(fit_te)[exposure]

  alpha_tbl <- purrr::map_dfr(mediators, function(med) {
    m_formula <- stats::as.formula(
      paste(
        med, "~", exposure,
        if (!is.null(cov_part)) paste("+", cov_part) else ""
      )
    )
    fit_m <- stats::lm(m_formula, data = data)
    tibble::tibble(mediator = med, alpha_X_to_M = stats::coef(fit_m)[exposure])
  })

  y_formula <- stats::as.formula(
    paste(
      outcome, "~", exposure, "+", paste(mediators, collapse = " + "),
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

  NIE <- sum(path_tbl$indirect_component)
  outcome_direct_coef <- stats::coef(fit_y)[exposure]
  NDE <- TE - NIE
  PM <- NIE / TE

  summary_tbl <- tibble::tibble(
    mediators = paste(mediators, collapse = " + "),
    covariates = ifelse(
      is.null(covariates) || length(covariates) == 0,
      "none",
      paste(covariates, collapse = " + ")
    ),
    TE = TE,
    NDE = NDE,
    NIE = NIE,
    PM = PM,
    outcome_direct_coef = outcome_direct_coef
  )

  list(summary = summary_tbl, path = path_tbl, fit_total = fit_te, fit_outcome = fit_y)
}

build_full_analysis_data <- function(pca_df) {
  pca_df |>
    dplyr::mutate(
      X_cat = factor(X_cat, levels = c("inside", "adjacent", "far")),
      X_adjacent = as.integer(X_cat == "adjacent"),
      X_far = as.integer(X_cat == "far")
    ) |>
    dplyr::select(Y, X_cat, X_adjacent, X_far, PC1_R, PC2_R, PC3, x_coord, y_coord) |>
    tidyr::drop_na()
}

decompose_linear_multix <- function(data,
                                    exposures,
                                    outcome,
                                    mediators,
                                    covariates = NULL) {
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
      outcome_direct_coef = stats::coef(fit_y)[exposure],
      NDE = TE - NIE,
      PM = NIE / TE,
      mediators = paste(mediators, collapse = " + "),
      covariates = ifelse(
        is.null(covariates) || length(covariates) == 0,
        "none",
        paste(covariates, collapse = " + ")
      )
    ) |>
    dplyr::select(exposure, mediators, covariates, TE, NDE, NIE, PM, outcome_direct_coef)

  list(summary = summary_tbl, path = path_tbl, fit_total = fit_te, fit_outcome = fit_y)
}

build_path_contribution_table <- function(res_nocov, res_coord) {
  dplyr::bind_rows(
    dplyr::mutate(res_nocov$path, model = "PC1-PC3, no controls"),
    dplyr::mutate(res_coord$path, model = "PC1-PC3, coordinate controls")
  ) |>
    dplyr::mutate(
      contrast = dplyr::case_when(
        exposure == "X_adjacent" ~ "adjacent vs inside",
        exposure == "X_far" ~ "far vs inside",
        TRUE ~ exposure
      ),
      dplyr::across(c(alpha_X_to_M, beta_M_to_Y, indirect_component), ~ round(.x, 5))
    ) |>
    dplyr::select(contrast, model, mediator, alpha_X_to_M, beta_M_to_Y, indirect_component)
}

compute_design_diagnostics <- function(data) {
  data <- dplyr::mutate(data, X_cat = factor(X_cat, levels = c("inside", "adjacent", "far")))

  X_design <- stats::model.matrix(~ X_cat + x_coord + y_coord, data = data)
  X_scaled <- X_design
  continuous_cols <- c("x_coord", "y_coord")
  X_scaled[, continuous_cols] <- scale(X_scaled[, continuous_cols])

  fit_rank_check <- stats::lm(stats::rnorm(nrow(data)) ~ X_cat + x_coord + y_coord, data = data)
  fit_vif <- stats::lm(stats::rnorm(nrow(data)) ~ X_cat + x_coord + y_coord, data = data)

  resid_PC1_coord <- stats::resid(stats::lm(PC1_R ~ X_cat + x_coord + y_coord, data = data))
  resid_PC2_coord <- stats::resid(stats::lm(PC2_R ~ X_cat + x_coord + y_coord, data = data))
  resid_PC3_coord <- stats::resid(stats::lm(PC3 ~ X_cat + x_coord + y_coord, data = data))

  resid_PC1_x <- stats::resid(stats::lm(PC1_R ~ X_cat, data = data))
  resid_PC2_x <- stats::resid(stats::lm(PC2_R ~ X_cat, data = data))
  resid_PC3_x <- stats::resid(stats::lm(PC3 ~ X_cat, data = data))

  list(
    X_design = X_design,
    design_dim = dim(X_design),
    design_rank = qr(X_design)$rank,
    design_ncol = ncol(X_design),
    full_rank = qr(X_design)$rank == ncol(X_design),
    design_colnames = colnames(X_design),
    alias = stats::alias(fit_rank_check),
    kappa_raw = kappa(X_design),
    kappa_scaled = kappa(X_scaled),
    vif = car::vif(fit_vif),
    residual_cor_coord = stats::cor(cbind(resid_PC1_coord, resid_PC2_coord, resid_PC3_coord)),
    residual_cor_x = stats::cor(cbind(resid_PC1_x, resid_PC2_x, resid_PC3_x))
  )
}

bootstrap_decomposition <- function(data,
                                    B = 1000,
                                    seed = 123,
                                    exposures = c("X_adjacent", "X_far"),
                                    outcome = "Y",
                                    mediators = c("PC1_R", "PC2_R", "PC3"),
                                    covariates = c("x_coord", "y_coord"),
                                    progress_callback = NULL) {
  set.seed(seed)

  n <- nrow(data)
  boot_summary <- vector("list", B)
  boot_path <- vector("list", B)

  for (b in seq_len(B)) {
    idx <- sample(seq_len(n), size = n, replace = TRUE)
    dat_b <- data[idx, ]

    fit_b <- decompose_linear_multix(
      data = dat_b,
      exposures = exposures,
      outcome = outcome,
      mediators = mediators,
      covariates = covariates
    )

    boot_summary[[b]] <- fit_b$summary
    boot_path[[b]] <- fit_b$path

    if (is.function(progress_callback)) {
      progress_callback(b, B)
    }
  }

  list(
    boot_summary = boot_summary,
    boot_path = boot_path,
    boot_summary_df = dplyr::bind_rows(boot_summary, .id = "boot_id"),
    boot_path_df = dplyr::bind_rows(boot_path, .id = "boot_id")
  )
}

summarize_bootstrap_ci <- function(boot_summary_df, boot_path_df) {
  summary_ci <- boot_summary_df |>
    dplyr::group_by(exposure) |>
    dplyr::summarise(
      TE_mean = mean(TE, na.rm = TRUE),
      TE_lwr = stats::quantile(TE, 0.025, na.rm = TRUE),
      TE_upr = stats::quantile(TE, 0.975, na.rm = TRUE),
      NDE_mean = mean(NDE, na.rm = TRUE),
      NDE_lwr = stats::quantile(NDE, 0.025, na.rm = TRUE),
      NDE_upr = stats::quantile(NDE, 0.975, na.rm = TRUE),
      NIE_mean = mean(NIE, na.rm = TRUE),
      NIE_lwr = stats::quantile(NIE, 0.025, na.rm = TRUE),
      NIE_upr = stats::quantile(NIE, 0.975, na.rm = TRUE),
      PM_mean = mean(PM, na.rm = TRUE),
      PM_lwr = stats::quantile(PM, 0.025, na.rm = TRUE),
      PM_upr = stats::quantile(PM, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  path_ci <- boot_path_df |>
    dplyr::group_by(exposure, mediator) |>
    dplyr::summarise(
      alpha_mean = mean(alpha_X_to_M, na.rm = TRUE),
      alpha_lwr = stats::quantile(alpha_X_to_M, 0.025, na.rm = TRUE),
      alpha_upr = stats::quantile(alpha_X_to_M, 0.975, na.rm = TRUE),
      beta_mean = mean(beta_M_to_Y, na.rm = TRUE),
      beta_lwr = stats::quantile(beta_M_to_Y, 0.025, na.rm = TRUE),
      beta_upr = stats::quantile(beta_M_to_Y, 0.975, na.rm = TRUE),
      indirect_mean = mean(indirect_component, na.rm = TRUE),
      indirect_lwr = stats::quantile(indirect_component, 0.025, na.rm = TRUE),
      indirect_upr = stats::quantile(indirect_component, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  list(summary_ci = summary_ci, path_ci = path_ci)
}

analyze_pc3_local_cluster <- function(pca_df,
                                      M_expr,
                                      pca_fit,
                                      pc3_negative_genes = c(
                                        "CXCL9", "CXCL10", "CCL5",
                                        "GZMK", "GNLY",
                                        "ISG15", "IFITM1",
                                        "TAP1", "NLRC5", "B2M",
                                        "HLA-DRA", "CD74"
                                      )) {
  pc3_cutoff <- pca_df |>
    dplyr::filter(X_cat == "far") |>
    dplyr::summarise(cutoff = stats::quantile(PC3, 0.10, na.rm = TRUE)) |>
    dplyr::pull(cutoff)

  pc3_target_cluster_initial <- pca_df |>
    dplyr::filter(X_cat == "far", PC3 <= pc3_cutoff, x_coord >= 7600, y_coord >= 7800)

  pc3_far_cutoff_10 <- pca_df |>
    dplyr::filter(X_cat == "far") |>
    dplyr::summarise(cutoff = stats::quantile(PC3, 0.10, na.rm = TRUE)) |>
    dplyr::pull(cutoff)

  pc3_far_low <- pca_df |>
    dplyr::filter(X_cat == "far", PC3 <= pc3_far_cutoff_10)

  pc3_target_cluster <- pc3_far_low |>
    dplyr::filter(x_coord >= 7600, y_coord >= 7800)

  gene_df <- as.data.frame(M_expr) |>
    tibble::rownames_to_column("spot_id")

  pca_df2 <- pca_df |>
    dplyr::mutate(spot_id = rownames(pca_fit$x)) |>
    dplyr::left_join(gene_df, by = "spot_id")

  comparison_df <- pca_df2 |>
    dplyr::mutate(
      pc3_cluster = dplyr::if_else(
        X_cat == "far" &
          x_coord >= 7600 &
          x_coord <= 8500 &
          y_coord >= 7800 &
          y_coord <= 9300 &
          PC3 <= pc3_far_cutoff_10,
        "PC3-negative local cluster",
        "Other far spots"
      )
    ) |>
    dplyr::filter(X_cat == "far")

  list(
    pc3_cutoff = pc3_cutoff,
    pc3_target_cluster_initial = pc3_target_cluster_initial,
    pc3_far_cutoff_10 = pc3_far_cutoff_10,
    pc3_far_low = pc3_far_low,
    pc3_negative_genes = pc3_negative_genes,
    pc3_target_cluster = pc3_target_cluster,
    pca_df2 = pca_df2,
    comparison_df = comparison_df,
    gene_summary = comparison_df |>
      dplyr::group_by(pc3_cluster) |>
      dplyr::summarise(dplyr::across(dplyr::all_of(pc3_negative_genes), ~ mean(.x, na.rm = TRUE))),
    immune_summary = comparison_df |>
      dplyr::group_by(pc3_cluster) |>
      dplyr::summarise(
        n = dplyr::n(),
        mean_immune_score = mean(immune_score, na.rm = TRUE),
        median_immune_score = stats::median(immune_score, na.rm = TRUE),
        sd_immune_score = stats::sd(immune_score, na.rm = TRUE)
      ),
    cluster_counts = dplyr::count(comparison_df, pc3_cluster)
  )
}

run_pca_mediation_pipeline <- function(rds_path = "data/raw/pt16_emory_GEX_immune_FULL_v2.rds",
                                       coord_x = "imagecol",
                                       coord_y = "imagerow",
                                       bootstrap_B = 1000,
                                       bootstrap_seed = 123,
                                       make_plots = TRUE) {
  gex <- load_spatial_mediation_rds(rds_path)
  analysis_df <- prepare_analysis_metadata(gex, coord_x = coord_x, coord_y = coord_y)

  pca <- fit_pca_mediators(gex$M_expr)
  pca_df <- build_pca_score_data(analysis_df, pca$pca_fit)
  scree_df <- build_scree_data(pca$var_explained)
  loading_df <- build_loading_data(pca$pca_fit)
  top_loading_genes <- extract_top_loading_genes(loading_df)

  spatial_inputs <- build_spatial_plot_inputs(pca_df, flip_y = TRUE)
  spatial_plots <- if (isTRUE(make_plots)) {
    build_spatial_pc_plots(
      spatial_inputs$plot_df,
      pca$var_explained,
      spatial_inputs$category_polygon_sf,
      flip_y = spatial_inputs$flip_y
    )
  } else {
    NULL
  }

  dat_full_allpc <- build_full_analysis_data(pca_df)

  res_full_pc123_coord <- decompose_linear_multix(
    data = dat_full_allpc,
    exposures = c("X_adjacent", "X_far"),
    outcome = "Y",
    mediators = c("PC1_R", "PC2_R", "PC3"),
    covariates = c("x_coord", "y_coord")
  )

  res_full_pc123_nocov <- decompose_linear_multix(
    data = dat_full_allpc,
    exposures = c("X_adjacent", "X_far"),
    outcome = "Y",
    mediators = c("PC1_R", "PC2_R", "PC3"),
    covariates = NULL
  )

  full3_path_contribution_table <- build_path_contribution_table(
    res_full_pc123_nocov,
    res_full_pc123_coord
  )

  diagnostics <- compute_design_diagnostics(dat_full_allpc)

  bootstrap <- bootstrap_decomposition(
    data = dat_full_allpc,
    B = bootstrap_B,
    seed = bootstrap_seed,
    exposures = c("X_adjacent", "X_far"),
    outcome = "Y",
    mediators = c("PC1_R", "PC2_R", "PC3"),
    covariates = c("x_coord", "y_coord")
  )

  bootstrap_ci <- summarize_bootstrap_ci(
    bootstrap$boot_summary_df,
    bootstrap$boot_path_df
  )

  pc3_local_cluster <- analyze_pc3_local_cluster(
    pca_df = pca_df,
    M_expr = gex$M_expr,
    pca_fit = pca$pca_fit
  )

  list(
    gex_names = names(gex),
    gex_structure = utils::capture.output(str(gex, max.level = 1)),
    X_cat = gex$X_cat,
    Y = gex$Y,
    meta = gex$meta,
    M_expr = gex$M_expr,
    analysis_df = analysis_df,
    M_scaled = pca$M_scaled,
    pca_fit = pca$pca_fit,
    pca_df = pca_df,
    var_explained = pca$var_explained,
    scree_df = scree_df,
    loading_df = loading_df,
    top_loading_genes = top_loading_genes,
    spatial_inputs = spatial_inputs,
    spatial_plots = spatial_plots,
    dat_full_allpc = dat_full_allpc,
    x_category_counts = table(dat_full_allpc$X_cat),
    res_full_pc123_coord = res_full_pc123_coord,
    res_full_pc123_nocov = res_full_pc123_nocov,
    full3_path_contribution_table = full3_path_contribution_table,
    diagnostics = diagnostics,
    bootstrap = bootstrap,
    summary_ci = bootstrap_ci$summary_ci,
    path_ci = bootstrap_ci$path_ci,
    pc3_local_cluster = pc3_local_cluster
  )
}
