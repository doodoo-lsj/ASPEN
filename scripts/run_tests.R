stopifnot(file.exists("app/app.R"), file.exists("scripts/tests/test_bundle_startup.R"))
tests <- c("test_bundle_startup.R", "test_structural_rho.R",
           "test_app_structural_rho.R", "test_app_selection_orientation.R",
           "test_app_result_export.R")
output <- "test-output"
dir.create(output, showWarnings = FALSE)
results <- lapply(tests, function(test) {
  cat("Running", test, "\n")
  start <- Sys.time()
  log <- file.path(output, paste0(test, ".log"))
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    shQuote(file.path("scripts/tests", test)), stdout = log, stderr = log)
  cat(if (status == 0L) "PASS" else "FAIL", test, "\n")
  data.frame(test = test, exit_status = status,
             elapsed_seconds = as.numeric(difftime(Sys.time(), start, units = "secs")))
})
results <- do.call(rbind, results)
write.csv(results, file.path(output, "test_results.csv"), row.names = FALSE)
if (any(results$exit_status != 0L)) stop("Tests failed; see test-output logs.")
cat("All app packaging tests passed. Full Bayesian MCMC was not rerun.\n")
