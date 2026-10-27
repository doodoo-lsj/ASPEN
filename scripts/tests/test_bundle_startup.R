root <- normalizePath(".")
app <- shiny::shinyAppDir("app")
server <- app$serverFuncSource()
env <- environment(server)
stopifnot(identical(get("project_root", env), root))
analysis <- get("example_analysis", env)
stopifnot(nrow(analysis$pca$pca_fit$x) == 1507L,
          nrow(analysis$pca$pca_fit$rotation) == 1985L,
          identical(get("selected_pc_ids", env)(list(), analysis), c("PC1", "PC2")))
artifacts <- get("bayesian_validated_artifacts", env)
stopifnot(length(artifacts) == 2L,
          identical(artifacts$PC1_PC2$pc_signs, c(PC1 = 1, PC2 = 1)),
          identical(artifacts$PC1_PC2_PC3$pc_signs, c(PC1 = -1, PC2 = -1, PC3 = 1)))
stopifnot(file.exists("R/app_result_export.R"))
cat("PASS: standalone root, both example artifacts, dimensions and default PC selection.\n")
