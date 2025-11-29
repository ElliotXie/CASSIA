devtools::install_github("ElliotXie/CASSIA/CASSIA_R", ref = "dev_test_branch")

library(CASSIA)

markers <- loadExampleMarkers(processed = FALSE)

outputdir="C:/Users/ellio/OneDrive - UW-Madison/CASSIA_enjoy/CASSIA/Test/16_manual_package_Test/results"

runCASSIA_pipeline(
  output_dir = file.path(outputdir, "test_output"),
  output_file_name = "test_output",
  tissue = "large intestine",
  species = "human",
  marker = markers,
  max_workers = 3,
  annotation_model = "google/gemini-2.5-flash",
  annotation_provider = "openrouter",
  score_model = "google/gemini-2.5-flash",
  score_provider = "openrouter",
  score_threshold = 98,
)
