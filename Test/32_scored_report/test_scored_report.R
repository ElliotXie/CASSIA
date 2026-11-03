# Run from the repository root with: Rscript Test/32_scored_report/test_scored_report.R
stopifnot(requireNamespace("reticulate", quietly = TRUE))
stopifnot(requireNamespace("jsonlite", quietly = TRUE))

Sys.setenv(RETICULATE_PYTHON = Sys.which("python"))
python_dir <- normalizePath("CASSIA_R/inst/python", winslash = "/", mustWork = TRUE)
cassia <- reticulate::import_from_path("CASSIA", path = python_dir)

test_dir <- tempfile("cassia_scored_report_")
dir.create(test_dir)
on.exit(unlink(test_dir, recursive = TRUE), add = TRUE)

input_file <- file.path(test_dir, "results_summary.csv")
output_file <- file.path(test_dir, "scored.csv")
write.csv(data.frame(
  "Cluster ID" = "cluster 1",
  "Predicted General Cell Type" = "T cell",
  "Marker List" = "CD3D",
  Tissue = "Large Intestine",
  Species = "Human",
  Score = 92,
  Scoring_Reasoning = "Score reason from R input",
  check.names = FALSE
), input_file, row.names = FALSE)
jsonlite::write_json(list("cluster 1" = list(
  annotations = list("Annotation from R companion JSON"),
  validations = list("VALIDATION PASSED: R companion JSON")
)), file.path(test_dir, "results_conversations.json"), auto_unbox = TRUE)

cassia$runCASSIA_score_batch(input_file, output_file, model = "test-model")
report <- paste(readLines(file.path(test_dir, "scored_report.html"), warn = FALSE),
                collapse = "\n")
stopifnot(grepl("Annotation from R companion JSON", report, fixed = TRUE))
stopifnot(grepl("Score reason from R input", report, fixed = TRUE))
