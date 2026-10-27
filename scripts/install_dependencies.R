packages <- c("shiny", "DT", "dplyr", "tidyr", "tibble", "ggplot2", "purrr",
              "digest", "readxl", "openxlsx", "zip", "posterior", "sf", "RANN",
              "testthat")
if ("--bayesian" %in% commandArgs(trailingOnly = TRUE)) {
  packages <- c(packages, "brms", "rstan")
}
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Packages still unavailable: ", paste(missing, collapse = ", "))
cat("Required packages are installed. New Bayesian fits also require a Stan toolchain.\n")
