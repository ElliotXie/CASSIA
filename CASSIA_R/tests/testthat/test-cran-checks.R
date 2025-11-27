# Basic CRAN submission tests
# These tests verify package integrity without requiring API calls or Python

test_that("CASSIA package loads", {
  expect_true(requireNamespace("CASSIA", quietly = TRUE))
})
test_that("Core annotation functions are exported", {
  expect_true(exists("runCASSIA", where = asNamespace("CASSIA")))
  expect_true(exists("runCASSIA_batch", where = asNamespace("CASSIA")))
  expect_true(exists("runCASSIA_pipeline", where = asNamespace("CASSIA")))
})

test_that("API key functions are exported", {
  expect_true(exists("setLLMApiKey", where = asNamespace("CASSIA")))
  expect_true(exists("setAnthropicApiKey", where = asNamespace("CASSIA")))
  expect_true(exists("set_openai_api_key", where = asNamespace("CASSIA")))
})

test_that("Setup functions are exported", {
  expect_true(exists("setup_cassia_env", where = asNamespace("CASSIA")))
  expect_true(exists("set_python_env", where = asNamespace("CASSIA")))
})

test_that("Data loading functions are exported", {
  expect_true(exists("loadExampleMarkers", where = asNamespace("CASSIA")))
  expect_true(exists("loadBuiltinMarkers", where = asNamespace("CASSIA")))
  expect_true(exists("listAvailableMarkers", where = asNamespace("CASSIA")))
})

test_that("Model functions are exported", {
  expect_true(exists("list_models", where = asNamespace("CASSIA")))
  expect_true(exists("get_recommended_model", where = asNamespace("CASSIA")))
  expect_true(exists("get_model_info", where = asNamespace("CASSIA")))
})

test_that("Example markers data loads correctly", {
  skip_on_cran()
  markers <- loadExampleMarkers()
  expect_true(is.data.frame(markers) || is.list(markers))
})

test_that("Available markers can be listed", {
  skip_on_cran()
  skip_if_not(reticulate::py_available(), "Python not available")
  available <- listAvailableMarkers()
  expect_true(is.character(available))
  expect_true(length(available) > 0)
})
