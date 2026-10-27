#!/usr/bin/env Rscript
# Match runApp's private app environment; sourcing app.R globally masks scope bugs.
app <- shiny::shinyAppDir("app")
server <- app$serverFuncSource()
stopifnot(!identical(environment(server), globalenv()),
          !exists("mediator_label", envir = globalenv(), inherits = FALSE))
example_analysis <- get("example_analysis", envir = environment(server))
selection_signs <- get("pc_signs_for_selection", envir = environment(server))
export_test_dir <- Sys.getenv("EXPORT_TEST_DIR", tempfile("export-tests-"))
dir.create(export_test_dir, recursive = TRUE, showWarnings = FALSE)
close_export <- function(x, y) stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = 1e-10)))
shiny::testServer(server, {
  stopifnot(!mediation_result_available())
  analysis_state(example_analysis)
  session$setInputs(analysis_framework = "bayesian", sensitivity_scenario = "common", sensitivity_active_rho = 0.15,
                   mediation_contrast_filter = "adjacent vs inside", mediation_mediator_filter = "PC1",
                   sensitivity_contrast_filter = "far vs inside", bootstrap_B = 20)
  confirmed_pca_spec(current_pca_spec())
  bundle <- result_export_bundle()
  stopifnot(length(unique(bundle$effects$contrast)) == 2L, length(unique(bundle$paths$mediator)) == 2L,
            length(unique(bundle$curves$rho)) == 13L, all(abs(bundle$selected$rho - 0.15) < 1e-12),
            all(c("common", "one_at_a_time") %in% bundle$curves$scenario))
  baseline <- get_bayesian_mediation_summary(bayesian_current_artifact())
  close_export(bundle$effects$estimate[bundle$effects$quantity == "Total indirect effect"], baseline$total_IIE_est)
  workbook <- file.path(export_test_dir, export_filename("xlsx"))
  downloaded <- output$download_result_tables
  stopifnot(basename(downloaded) == paste0("Bayesian_", Sys.Date(), ".xlsx"),
            grepl("Download Tables", output$result_export_controls$html, fixed = TRUE))
  file.copy(downloaded, workbook, overwrite = TRUE)
  stopifnot(length(openxlsx::getSheetNames(workbook)) == 7L)
  exported_headers <- unlist(lapply(openxlsx::getSheetNames(workbook), function(sheet) {
    names(openxlsx::read.xlsx(workbook, sheet = sheet, rows = 1))
  }), use.names = FALSE)
  stopifnot(!any(c("Posterior.Pr.effect...0.", "Feasible.fraction") %in% exported_headers))
  readback <- openxlsx::read.xlsx(workbook, sheet = "Overall Effects")
  stopifnot(is.numeric(readback$Estimate))
  close_export(readback$Estimate, bundle$effects$estimate)
  context <- result_context_figures()
  expected_context <- c("Spatial_Exposure.png", "Spatial_Outcome.png", "PCA_Explained_Variance.png",
                        "PC1_Spatial.png", "PC1_Loadings.png", "PC2_Spatial.png", "PC2_Loadings.png")
  stopifnot(setequal(names(context), expected_context),
            identical(context$Spatial_Exposure.png$data, overview_exposure_plot()$data),
            identical(context$PCA_Explained_Variance.png$data, explained_variance_plot()$data))
  write_result_figures(bundle, file.path(export_test_dir, "Bayesian_common.zip"), context)
  session$setInputs(sensitivity_scenario = "one_at_a_time", sensitivity_one_mediator = "PC2")
  one <- result_export_bundle()
  stopifnot(all(one$selected$varied_mediator == "PC2"))
  downloaded_figures <- output$download_result_figures
  stopifnot(basename(downloaded_figures) == paste0("Bayesian_", Sys.Date(), ".zip"))
  file.copy(downloaded_figures, file.path(export_test_dir, "Bayesian_one.zip"), overwrite = TRUE)
  common_files <- utils::unzip(file.path(export_test_dir, "Bayesian_common.zip"), list = TRUE)$Name
  one_files <- utils::unzip(file.path(export_test_dir, "Bayesian_one.zip"), list = TRUE)$Name
  stopifnot(length(common_files) == 10L, length(one_files) == 11L,
            all(expected_context %in% common_files), all(expected_context %in% one_files))
  saveRDS(one, file.path(export_test_dir, "bayesian_bundle.rds"))

  session$setInputs(analysis_framework = "frequentist")
  freq <- list(key = frequentist_current_cache_key(), result = make_frequentist_dynamic_fit(example_analysis, confirmed_selected_pc_ids(), confirmed_pc_sign_vector()))
  frequentist_fit_state(freq)
  ensure_frequentist_bootstrap(freq)
  before <- frequentist_bootstrap_cache()
  f <- result_export_bundle()
  stopifnot(identical(before, frequentist_bootstrap_cache()), all(is.finite(f$effects$lower)), all(is.finite(f$paths$lower)),
            all(is.finite(f$curves$lower)), all(f$selected$varied_mediator == "PC2"))
  close_export(f$effects$estimate[f$effects$quantity == "Total indirect effect"], freq$result$fit$summary$NIE)
  close_export(f$effects$lower[f$effects$quantity == "Total indirect effect"], current_bootstrap_summary_ci()$NIE_lwr)
  ci <- current_bootstrap_sensitivity_ci()$one_summary_ci
  for (i in seq_len(nrow(ci))) {
    row <- subset(f$curves, scenario == "one_at_a_time" & contrast == ci$contrast[i] &
                    varied_mediator == ci$varied_mediator[i] & abs(rho - ci$rho[i]) < 1e-10 & quantity == "Total indirect effect")
    stopifnot(nrow(row) == 1L)
    close_export(row$lower, ci$adjusted_indirect_lwr[i])
  }
  downloaded <- output$download_result_tables
  downloaded_figures <- output$download_result_figures
  stopifnot(basename(downloaded) == paste0("Frequentist_", Sys.Date(), ".xlsx"),
            basename(downloaded_figures) == paste0("Frequentist_", Sys.Date(), ".zip"))
  file.copy(downloaded, file.path(export_test_dir, basename(downloaded)), overwrite = TRUE)
  file.copy(downloaded_figures, file.path(export_test_dir, "Frequentist_one.zip"), overwrite = TRUE)
  frequentist_bootstrap_cache(list())
  without_ci <- result_export_bundle()
  stopifnot(all(is.na(without_ci$effects$lower)), identical(f$effects$estimate, without_ci$effects$estimate))
  # Spatial/feature exports use the confirmed set, including nonconsecutive PCs.
  confirmed_pca_spec(list(selected_pcs = c("PC1", "PC3"),
                          pc_signs = selection_signs(c("PC1", "PC3"), analysis())))
  nonconsecutive <- result_context_figures()
  stopifnot("PC3_Loadings.png" %in% names(nonconsecutive), !"PC2_Loadings.png" %in% names(nonconsecutive))
  confirmed_pca_spec(list(selected_pcs = paste0("PC", 1:6),
                          pc_signs = selection_signs(paste0("PC", 1:6), analysis())))
  six <- result_context_figures()
  stopifnot(length(six) == 15L, "PC6_Spatial.png" %in% names(six))
  confirmed_pca_spec(list(selected_pcs = paste0("PC", 1:3),
                          pc_signs = selection_signs(paste0("PC", 1:3), analysis())))
  three <- result_context_figures()
  close_export(three$PC1_Spatial.png$layers[[1]]$data$PC1, oriented_spatial_df()$PC1 * -1)
  session$setInputs(loading_n = 50)
  tall <- result_context_figures()$PC1_Loadings.png
  stopifnot(nrow(tall$data) == 100L, attr(tall, "export_height") > 10)
  confirmed_pca_spec(NULL)
  session$flushReact()
  stopifnot(!mediation_result_available(), is.null(output$result_export_controls))
})
test_curve <- data.frame(contrast = "test", scenario = "common", varied_mediator = "All mediators", mediator = "All mediators",
                         quantity = "Total indirect effect", rho = c(-0.1, 0, 0.1), estimate = c(-1, NA, 1))
stopifnot(is.na(export_zero_crossings(test_curve)$rho))
test_curve$estimate <- c(-1, 0, 1)
stopifnot(identical(export_zero_crossings(test_curve)$rho, 0))
cat("PASS: all-contrast tables, numeric workbook round trip, both frameworks/scenarios, cached CI reuse, spatial/PCA/loading ZIPs and arbitrary confirmed PCs.\n", export_test_dir, "\n")
