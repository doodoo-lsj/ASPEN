# Data input and setup helpers for the Shiny app.
#
# This layer keeps raw metadata, mediator feature matrices, and selected PC
# mediators conceptually separate. It performs no Bayesian fitting.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

available_excel_sheets <- function(path) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Excel upload requires the readxl package.", call. = FALSE)
  }
  readxl::excel_sheets(path)
}

read_uploaded_table <- function(path, sheet = NULL) {
  if (is.null(path) || !file.exists(path)) {
    stop("Uploaded file is not available.", call. = FALSE)
  }

  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    out <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (ext == "tsv") {
    out <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (ext == "txt") {
    lines <- readLines(path, n = 5, warn = FALSE)
    tab_count <- sum(grepl("\t", lines, fixed = TRUE))
    comma_count <- sum(grepl(",", lines, fixed = TRUE))
    if (tab_count > comma_count && tab_count > 0) {
      out <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    } else if (comma_count > tab_count && comma_count > 0) {
      out <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    } else {
      out <- utils::read.table(path, header = TRUE, sep = "", check.names = FALSE, stringsAsFactors = FALSE)
    }
  } else if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Excel upload requires the readxl package. Please upload CSV/TSV/TXT or install readxl.", call. = FALSE)
    }
    sheet <- sheet %||% readxl::excel_sheets(path)[[1]]
    out <- as.data.frame(readxl::read_excel(path, sheet = sheet), check.names = FALSE)
  } else {
    stop("Unsupported file format. Please upload CSV, TSV, TXT, XLSX, or XLS.", call. = FALSE)
  }

  validate_basic_table(out, "uploaded table")
  out
}

validate_basic_table <- function(data, label) {
  if (!is.data.frame(data)) {
    stop(label, " is not a tabular data frame.", call. = FALSE)
  }
  if (nrow(data) == 0 || ncol(data) == 0) {
    stop(label, " is empty.", call. = FALSE)
  }
  if (anyDuplicated(names(data)) > 0) {
    stop(label, " has duplicated column names.", call. = FALSE)
  }
  invisible(TRUE)
}

numeric_column_names <- function(data) {
  names(data)[vapply(data, is.numeric, logical(1))]
}

table_preview <- function(data, n_rows = 6, n_cols = 8) {
  data[seq_len(min(nrow(data), n_rows)), seq_len(min(ncol(data), n_cols)), drop = FALSE]
}

column_type_table <- function(data) {
  tibble::tibble(
    Column = names(data),
    Type = vapply(data, function(x) paste(class(x), collapse = "/"), character(1))
  )
}

example_data_inputs <- function(rds_path) {
  gex <- load_spatial_mediation_rds(rds_path)
  obs_id <- rownames(gex$M_expr)
  if (is.null(obs_id)) {
    obs_id <- names(gex$Y)
  }
  if (is.null(obs_id)) {
    obs_id <- as.character(seq_along(gex$Y))
  }

  main_data <- as.data.frame(gex$meta, stringsAsFactors = FALSE)
  if (!"barcode" %in% names(main_data)) {
    main_data$barcode <- obs_id
  }
  main_data$Y <- as.numeric(gex$Y)
  if (!"X_cat" %in% names(main_data)) {
    main_data$X_cat <- gex$X_cat
  }

  feature_matrix <- as.data.frame(gex$M_expr, check.names = FALSE)
  feature_matrix <- tibble::add_column(feature_matrix, barcode = obs_id, .before = 1)

  list(
    data_source = "example",
    main_data = main_data,
    feature_data = feature_matrix,
    defaults = list(
      observation_id = "barcode",
      exposure = "X_cat",
      outcome = "Y",
      spatial_x = "imagecol",
      spatial_y = "imagerow",
      covariates = character(0),
      reference = "inside",
      feature_observation_id = "barcode",
      excluded_feature_columns = character(0)
    )
  )
}

safe_dummy_name <- function(level) {
  paste0("X_", make.names(as.character(level)))
}

make_exposure_dummy_names <- function(levels, reference) {
  contrast_levels <- setdiff(levels, reference)
  stats::setNames(vapply(contrast_levels, safe_dummy_name, character(1)), contrast_levels)
}

check_selected_columns_exist <- function(data, columns, label) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    stop(label, " selected columns are missing: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

check_unique_roles <- function(roles) {
  required <- unlist(roles[c("exposure", "outcome", "spatial_x", "spatial_y")], use.names = FALSE)
  conflict <- required[duplicated(required)]
  if (length(conflict) > 0) {
    stop("Variable roles conflict. A column cannot be assigned to incompatible required roles: ", paste(unique(conflict), collapse = ", "), call. = FALSE)
  }
  covariates <- roles$covariates %||% character(0)
  cov_conflict <- intersect(covariates, required)
  if (length(cov_conflict) > 0) {
    stop("Additional covariates cannot reuse columns assigned to exposure, outcome, or spatial coordinates: ", paste(unique(cov_conflict), collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

validate_ids <- function(values, label) {
  if (any(is.na(values) | values == "")) {
    stop(label, " contains missing observation IDs.", call. = FALSE)
  }
  if (anyDuplicated(values) > 0) {
    stop(label, " contains duplicated observation IDs.", call. = FALSE)
  }
  invisible(TRUE)
}

build_analysis_spec <- function(data_source,
                                main_data,
                                feature_data,
                                roles,
                                feature_id = NULL,
                                excluded_feature_columns = character(0)) {
  validate_basic_table(main_data, "Main dataset")
  validate_basic_table(feature_data, "Mediator gene / feature matrix")

  required_role_names <- c("exposure", "outcome", "spatial_x", "spatial_y", "reference")
  missing_roles <- setdiff(required_role_names, names(roles))
  if (length(missing_roles) > 0) {
    stop("Missing required variable-role selections: ", paste(missing_roles, collapse = ", "), call. = FALSE)
  }

  role_columns <- unlist(roles[c("observation_id", "exposure", "outcome", "spatial_x", "spatial_y", "covariates")], use.names = FALSE)
  role_columns <- role_columns[nzchar(role_columns)]
  check_selected_columns_exist(main_data, role_columns, "Main dataset")
  check_unique_roles(roles)

  if (!roles$exposure %in% names(main_data)) {
    stop("Exposure column is missing.", call. = FALSE)
  }
  exposure_values <- main_data[[roles$exposure]]
  if (any(is.na(exposure_values))) {
    stop("Exposure contains missing values.", call. = FALSE)
  }
  exposure_levels <- unique(as.character(exposure_values))
  if (length(exposure_levels) < 2) {
    stop("Exposure must contain at least two observed categories.", call. = FALSE)
  }
  if (!roles$reference %in% exposure_levels) {
    stop("Reference category is not an observed exposure level.", call. = FALSE)
  }
  exposure_levels <- c(roles$reference, setdiff(exposure_levels, roles$reference))

  for (nm in c("outcome", "spatial_x", "spatial_y")) {
    col <- roles[[nm]]
    if (!is.numeric(main_data[[col]])) {
      stop(c(outcome = "Outcome", spatial_x = "Spatial X coordinate", spatial_y = "Spatial Y coordinate")[[nm]],
           " must be numeric.", call. = FALSE)
    }
    if (any(!is.finite(main_data[[col]]))) {
      stop(c(outcome = "Outcome", spatial_x = "Spatial X coordinate", spatial_y = "Spatial Y coordinate")[[nm]],
           " contains missing or non-finite values.", call. = FALSE)
    }
  }

  main_id_col <- roles$observation_id %||% ""
  feature_id <- feature_id %||% ""
  if (nzchar(main_id_col) && nzchar(feature_id)) {
    check_selected_columns_exist(feature_data, feature_id, "Mediator gene / feature matrix")
    main_ids <- as.character(main_data[[main_id_col]])
    feature_ids <- as.character(feature_data[[feature_id]])
    validate_ids(main_ids, "Main dataset observation ID")
    validate_ids(feature_ids, "Feature matrix observation ID")
    missing_in_feature <- setdiff(main_ids, feature_ids)
    missing_in_main <- setdiff(feature_ids, main_ids)
    if (length(missing_in_feature) > 0 || length(missing_in_main) > 0) {
      stop(
        "Observation IDs do not match. Missing in feature matrix: ", length(missing_in_feature),
        "; missing in main dataset: ", length(missing_in_main), ".",
        call. = FALSE
      )
    }
    feature_data <- feature_data[match(main_ids, feature_ids), , drop = FALSE]
    observation_ids <- main_ids
    match_warning <- NULL
  } else if (!nzchar(main_id_col) && !nzchar(feature_id)) {
    if (nrow(main_data) != nrow(feature_data)) {
      stop("No observation IDs were supplied and row counts differ, so observations cannot be matched.", call. = FALSE)
    }
    observation_ids <- as.character(seq_len(nrow(main_data)))
    match_warning <- "No observation ID was provided. The app is assuming that the rows of the main dataset and mediator feature matrix are in exactly the same order."
  } else {
    stop("Select observation ID columns in both datasets, or leave both as None for positional matching.", call. = FALSE)
  }

  excluded_feature_columns <- unique(c(excluded_feature_columns, feature_id))
  feature_columns <- setdiff(names(feature_data), excluded_feature_columns)
  numeric_features <- feature_columns[vapply(feature_data[feature_columns], is.numeric, logical(1))]
  nonnumeric_features <- setdiff(feature_columns, numeric_features)
  if (length(numeric_features) == 0) {
    stop("No usable numeric mediator feature columns are available for PCA.", call. = FALSE)
  }

  feature_matrix <- as.matrix(feature_data[, numeric_features, drop = FALSE])
  storage.mode(feature_matrix) <- "double"
  if (any(!is.finite(feature_matrix))) {
    stop("Mediator feature matrix contains missing or non-finite values that prevent PCA.", call. = FALSE)
  }
  feature_sds <- apply(feature_matrix, 2, stats::sd)
  zero_variance <- names(feature_sds)[!is.finite(feature_sds) | feature_sds == 0]
  if (length(zero_variance) > 0) {
    stop("Mediator feature matrix contains zero-variance features: ", paste(utils::head(zero_variance, 10), collapse = ", "), call. = FALSE)
  }
  if (nrow(feature_matrix) < 3 || ncol(feature_matrix) < 3) {
    stop("Insufficient observations or features for the current PC1-PC3 PCA mediation workflow.", call. = FALSE)
  }
  rownames(feature_matrix) <- observation_ids

  dummy_names <- make_exposure_dummy_names(exposure_levels, roles$reference)
  contrast_labels <- stats::setNames(
    paste(names(dummy_names), "vs", roles$reference),
    unname(dummy_names)
  )

  metadata <- tibble::tibble(
    observation_id = observation_ids,
    X_cat = factor(as.character(main_data[[roles$exposure]]), levels = exposure_levels),
    Y = as.numeric(main_data[[roles$outcome]]),
    x_coord = as.numeric(main_data[[roles$spatial_x]]),
    y_coord = as.numeric(main_data[[roles$spatial_y]])
  )
  for (i in seq_along(dummy_names)) {
    metadata[[dummy_names[[i]]]] <- as.integer(metadata$X_cat == names(dummy_names)[[i]])
  }

  covariate_columns <- roles$covariates %||% character(0)
  covariate_map <- character(0)
  if (length(covariate_columns) > 0) {
    for (cov in covariate_columns) {
      canonical <- make.names(paste0("cov_", cov), unique = TRUE)
      metadata[[canonical]] <- main_data[[cov]]
      covariate_map[[canonical]] <- cov
    }
  }

  list(
    data_source = data_source,
    observation_id = if (nzchar(main_id_col)) main_id_col else NULL,
    exposure = list(
      variable = roles$exposure,
      reference = roles$reference,
      levels = exposure_levels,
      dummy_names = unname(dummy_names),
      contrast_labels = contrast_labels
    ),
    outcome = roles$outcome,
    coordinates = list(x = roles$spatial_x, y = roles$spatial_y),
    covariates = covariate_map,
    mediator_features = list(
      observation_id = if (nzchar(feature_id)) feature_id else NULL,
      feature_names = colnames(feature_matrix),
      excluded_columns = excluded_feature_columns,
      nonnumeric_excluded = nonnumeric_features,
      n_features = ncol(feature_matrix)
    ),
    matching = list(
      main_n = nrow(main_data),
      feature_n = nrow(feature_data),
      matched_n = length(observation_ids),
      warning = match_warning
    ),
    analysis_data = list(
      metadata = metadata,
      feature_matrix = feature_matrix,
      exposure_dummies = unname(dummy_names),
      covariates = c("x_coord", "y_coord", names(covariate_map))
    ),
    original_labels = list(
      exposure = roles$exposure,
      outcome = roles$outcome,
      spatial_x = roles$spatial_x,
      spatial_y = roles$spatial_y,
      covariates = covariate_map
    )
  )
}

validate_example_bayesian_config <- function(spec) {
  identical(spec$data_source, "example") &&
    identical(spec$exposure$variable, "X_cat") &&
    identical(spec$outcome, "Y") &&
    identical(spec$coordinates$x, "imagecol") &&
    identical(spec$coordinates$y, "imagerow") &&
    identical(spec$exposure$reference, "inside") &&
    identical(unname(spec$exposure$dummy_names), c("X_adjacent", "X_far")) &&
    length(spec$covariates) == 0
}

validation_status_table <- function(spec) {
  warning_text <- spec$matching$warning
  tibble::tibble(
    Check = c(
      "Main dataset loaded",
      "Exposure is categorical",
      "Outcome is numeric",
      "Coordinates are numeric",
      "Observations matched",
      "Numeric mediator features available",
      "Reference category",
      "Ready for PCA"
    ),
    Status = c(
      paste0("OK: ", spec$matching$main_n, " observations"),
      paste0("OK: ", length(spec$exposure$levels), " levels"),
      "OK",
      "OK",
      paste0("OK: ", spec$matching$matched_n, " matched", if (!is.null(warning_text)) " (positional matching warning)" else ""),
      paste0("OK: ", spec$mediator_features$n_features, " features"),
      paste0("OK: ", spec$exposure$reference),
      "OK"
    )
  )
}

build_full_analysis_data_from_spec <- function(pca_df, spec, mediators = c("PC1_R", "PC2_R", "PC3")) {
  required <- c("Y", "X_cat", spec$analysis_data$exposure_dummies, mediators, "x_coord", "y_coord", names(spec$covariates))
  missing <- setdiff(required, names(pca_df))
  if (length(missing) > 0) {
    stop("Prepared analysis data is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }

  pca_df |>
    dplyr::select(dplyr::all_of(required)) |>
    tidyr::drop_na()
}

compute_analysis_from_spec <- function(spec, mediation_mediators = c("PC1_R", "PC2_R", "PC3")) {
  # Keep the validated example axes exact; uploaded data only needs the six axes exposed by the app.
  pca_method <- if (identical(spec$data_source, "example")) "exact" else "randomized"
  pca <- fit_pca_mediators(
    spec$analysis_data$feature_matrix,
    method = pca_method,
    rank = 6L
  )
  pca_df <- build_pca_score_data(spec$analysis_data$metadata, pca$pca_fit)
  scree_df <- build_scree_data(pca$var_explained)
  loading_df <- build_loading_data(pca$pca_fit)
  # Preserve historical column aliases while using unmodified PCA axes in the app.
  # The publication pipeline keeps its separately specified reversed orientation.
  pca_df$PC1_R <- pca_df$PC1_original
  pca_df$PC2_R <- pca_df$PC2_original
  loading_df$PC1_R <- loading_df$PC1_original
  loading_df$PC2_R <- loading_df$PC2_original
  flip_y <- identical(spec$data_source, "example")
  spatial_inputs <- build_spatial_plot_inputs(pca_df, flip_y = flip_y)
  spatial_plots <- build_spatial_pc_plots(
    spatial_inputs$plot_df,
    pca$var_explained,
    spatial_inputs$category_polygon_sf,
    flip_y = spatial_inputs$flip_y
  )

  dat_full_allpc <- build_full_analysis_data_from_spec(
    pca_df,
    spec,
    mediators = mediation_mediators
  )

  mediation_coord <- decompose_linear_multix(
    data = dat_full_allpc,
    exposures = spec$analysis_data$exposure_dummies,
    outcome = "Y",
    mediators = mediation_mediators,
    covariates = spec$analysis_data$covariates
  )

  list(
    spec = spec,
    pca = pca,
    pca_df = pca_df,
    scree_df = scree_df,
    loading_df = loading_df,
    top_loading_genes = extract_top_loading_genes(loading_df),
    spatial_inputs = spatial_inputs,
    spatial_plots = spatial_plots,
    dat_full_allpc = dat_full_allpc,
    mediation_coord = mediation_coord,
    exposures = spec$analysis_data$exposure_dummies,
    mediators = mediation_mediators,
    covariates = spec$analysis_data$covariates,
    contrast_labels = spec$exposure$contrast_labels,
    bayesian_artifact_allowed = validate_example_bayesian_config(spec)
  )
}
