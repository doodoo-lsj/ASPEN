#!/usr/bin/env Rscript
# Exercise the nested Shiny source path and the result adapters used by the UI.
options(warn=2)
app <- source("app/app.R")$value
stopifnot(inherits(app, "shiny.appobj"))
for (entry in bayesian_validated_artifacts) {
  a <- entry$artifact
  stopifnot(identical(a$sensitivity_A$parameterization, structural_rho_version()))
  for (mode in c("common", "one_at_a_time")) {
    med <- if (mode == "one_at_a_time") bayes_A_bundle_mediators(a$sensitivity_A$residual_draws)[1] else NULL
    selected <- bayesian_sensitivity_selected_cache(a, paste(entry$pc_ids, collapse="_"), .15, mode, med)
    stopifnot(is.list(selected))
    curve <- bayesian_sensitivity_curve_cache(a, paste(entry$pc_ids, collapse="_"), mode, med)
    stopifnot(nrow(curve$summary) > 0,
              all(c("total_IIE_est", "direct_est", "TE_est") %in% names(curve$summary)),
              all(is.finite(curve$summary$total_IIE_est)),
              max(abs(curve$summary$TE_est - curve$summary$direct_est -
                        curve$summary$total_IIE_est)) < 1e-12)
  }
  diagnostics <- bayesian_app_mcmc_diagnostics(a)
  stopifnot(nrow(diagnostics) == 3L, "overall_pass" %in% names(diagnostics))
}
shiny::testServer(server, {
  analysis_state(example_analysis)
  confirmed_pca_spec(list(selected_pcs=c("PC1","PC2","PC3"),
                          pc_signs=pc_signs_for_selection(c("PC1","PC2","PC3"), example_analysis)))
  result <- bayesian_current_result()
  stopifnot(identical(result$source, "validated_artifact"),
            grepl(structural_rho_version(), result$key, fixed=TRUE),
            identical(bayesian_current_artifact()$baseline, bayesian_app_artifact$baseline))
})
cat("PASS: Shiny source, validated artifact adapters, common/OAT caches, diagnostics and server selection.\n")
