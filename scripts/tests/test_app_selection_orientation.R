#!/usr/bin/env Rscript
options(warn = 2)
source("app/app.R")
close <- function(x, y, tol = 1e-9) {
  stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol)))
}
pcs <- paste0("PC", 1:6)
uploaded <- example_analysis
uploaded$spec$data_source <- "upload"
uploaded$bayesian_artifact_allowed <- FALSE
modified_example <- example_analysis
modified_example$spec$exposure$reference <- "far"
modified_example$bayesian_artifact_allowed <- FALSE
modified_features <- example_analysis
modified_features$spec$analysis_data$feature_matrix <-
  modified_features$spec$analysis_data$feature_matrix[, -1, drop = FALSE]
original <- setNames(rep(1, 6), pcs)
for (mask in 1:63) {
  selected <- pcs[as.logical(intToBits(mask)[1:6])]
  expected <- setNames(rep(1, 6), pcs)
  if (setequal(selected, pcs[1:3])) expected[1:2] <- -1
  stopifnot(identical(pc_signs_for_selection(selected, example_analysis), expected),
            identical(pc_signs_for_selection(rev(selected), example_analysis), expected),
            identical(pc_signs_for_selection(selected), original),
            identical(pc_signs_for_selection(selected, modified_features), original))
  for (other in list(uploaded, modified_example)) {
    other_signs <- pc_sign_vector(list(pca_mediation_pcs = selected), other)
    stopifnot(identical(other_signs, original))
    other_prepared <- bayesian_dynamic_prepare_data(other, selected, other_signs)
    for (j in seq_along(selected)) {
      close(other_prepared$data[[paste0("M", j)]], other$pca$pca_fit$x[, selected[j]])
    }
  }
  signs <- pc_sign_vector(list(pca_mediation_pcs = selected), example_analysis)
  stopifnot(identical(signs, expected))
  prepared <- bayesian_dynamic_prepare_data(example_analysis, selected, signs)
  scores <- orient_spatial_scores(example_analysis$spatial_inputs$plot_df,
                                  signs, example_analysis$pca$pca_fit)
  loadings <- orient_loading_data(example_analysis$loading_df, signs,
                                  example_analysis$pca$pca_fit)
  for (j in seq_along(selected)) {
    pc <- selected[j]
    close(prepared$data[[paste0("M", j)]],
          example_analysis$pca$pca_fit$x[, pc] * signs[[pc]])
    close(scores[[pc]], prepared$data[[paste0("M", j)]])
    close(loadings[[pc]], example_analysis$pca$pca_fit$rotation[, pc] * signs[[pc]])
  }
}
reversed <- pc_signs_for_selection(pcs[1:3], example_analysis)
f3 <- make_frequentist_example_fit(example_analysis, reversed)
f0 <- make_frequentist_dynamic_fit(example_analysis, pcs[1:3], original)
close(f3$fit$summary$NIE, f0$fit$summary$NIE)
close(f3$fit$summary$NDE, f0$fit$summary$NDE)
close(f3$fit$summary$TE, f0$fit$summary$TE)
close(f3$fit$path$alpha_X_to_M,
      f0$fit$path$alpha_X_to_M * reversed[f0$fit$path$mediator])
close(f3$fit$path$beta_M_to_Y,
      f0$fit$path$beta_M_to_Y * reversed[f0$fit$path$mediator])
stopifnot(!identical(frequentist_dynamic_cache_key(example_analysis, pcs[1:3], reversed),
                     frequentist_dynamic_cache_key(example_analysis, pcs[1:3], original)))
shiny::testServer(server, {
  analysis_state(example_analysis)
  for (selected in list(pcs[1:2], pcs[1:3], pcs[1:4], pcs[c(1, 3, 4)], pcs[1])) {
    session$setInputs(pca_mediation_pcs = selected)
    spec <- current_pca_spec()
    stopifnot(identical(spec$pc_signs, pc_signs_for_selection(selected, example_analysis)))
    confirmed_pca_spec(spec)
    entry <- bayesian_validated_artifact_for_selection(selected)
    if (is.null(entry)) {
      stopifnot(!bayesian_artifact_available(), !frequentist_example_default_available())
    } else {
      stopifnot(bayesian_artifact_available(),
                identical(bayesian_current_artifact(), entry$artifact))
      shown <- bayesian_mediation_path()
      expected <- get_bayesian_mediation_path(entry$artifact)
      close(shown$alpha_est, expected$alpha_est)
      close(shown$beta_est, expected$beta_est)
      close(shown$alpha_lwr, expected$alpha_lwr)
      close(shown$beta_upr, expected$beta_upr)
      if (setequal(selected, pcs[1:3])) {
        stopifnot(frequentist_example_default_available())
        close(frequentist_current_result()$result$prepared$data$M1,
              -example_analysis$pca$pca_fit$x[, 1])
        close(oriented_spatial_df()$PC1, -example_analysis$pca$pca_fit$x[, 1])
        close(oriented_loading_df()$PC2, -example_analysis$pca$pca_fit$rotation[, 2])
      }
    }
  }
  # Never reuse the publication artifact for an incompatible confirmed orientation.
  confirmed_pca_spec(list(selected_pcs = pcs[1:3], pc_signs = original))
  stopifnot(!bayesian_artifact_available(), !frequentist_example_default_available())
  session$setInputs(pca_mediation_pcs = pcs[1:3])
  for (other in list(uploaded, modified_example)) {
    analysis_state(other)
    spec <- current_pca_spec()
    stopifnot(identical(spec$pc_signs, original))
    confirmed_pca_spec(spec)
    stopifnot(!bayesian_artifact_available(), !frequentist_example_default_available())
    close(oriented_spatial_df()$PC1, other$pca$pca_fit$x[, 1])
    close(oriented_loading_df()$PC2, other$pca$pca_fit$rotation[, 2])
  }
  analysis_state(example_analysis)
  stopifnot(identical(current_pca_spec()$pc_signs, reversed))
})
cat("PASS: all 63 PC subsets in Emory, uploaded and modified-example contexts; reversal only for the Emory publication three-PC result; model inputs, maps/loadings, UI paths and cache guards.\n")
