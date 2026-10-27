manifest <- read.csv("docs/bundle_checksums.csv", stringsAsFactors = FALSE)
stopifnot(all(file.exists(manifest$path)))
observed <- unname(tools::md5sum(manifest$path))
if (anyNA(observed) || !identical(observed, manifest$md5)) {
  stop("Bundle files differ from the packaged checksums. Review changes before updating the manifest.")
}
cat("PASS:", nrow(manifest), "app source/data/test files match the packaged checksums.\n")
