# Test script to verify backward compatibility
# This script tests that existing code patterns continue to work

# Test 1: setup_cassia_env() with no parameters (should work with new fallback logic)
message("Test 1: setup_cassia_env() with no parameters")
tryCatch({
  # This should now try virtualenv first, then fall back to conda
  CASSIA::setup_cassia_env()
  message("✓ setup_cassia_env() with no parameters works")
}, error = function(e) {
  message("✗ setup_cassia_env() with no parameters failed: ", e$message)
})

# Test 2: setup_cassia_env() with conda_env parameter (existing pattern)
message("\nTest 2: setup_cassia_env() with conda_env parameter")
tryCatch({
  # This should work exactly as before
  CASSIA::setup_cassia_env(conda_env = "cassia_env")
  message("✓ setup_cassia_env() with conda_env parameter works")
}, error = function(e) {
  message("✗ setup_cassia_env() with conda_env parameter failed: ", e$message)
})

# Test 3: set_python_env() with existing environment
message("\nTest 3: set_python_env() with environment name")
tryCatch({
  # This should work with both virtualenv and conda environments
  result <- CASSIA::set_python_env("cassia_env")
  if (result) {
    message("✓ set_python_env() works")
  } else {
    message("✗ set_python_env() returned FALSE")
  }
}, error = function(e) {
  message("✗ set_python_env() failed: ", e$message)
})

# Test 4: check_python_env() should work regardless of environment type
message("\nTest 4: check_python_env()")
tryCatch({
  result <- CASSIA::check_python_env()
  if (result) {
    message("✓ check_python_env() works")
  } else {
    message("✗ check_python_env() returned FALSE")
  }
}, error = function(e) {
  message("✗ check_python_env() failed: ", e$message)
})

# Test 5: New method parameter functionality
message("\nTest 5: New method parameter functionality")
tryCatch({
  # Test explicit virtualenv method
  CASSIA::setup_cassia_env(conda_env = "cassia_env_venv", method = "virtualenv")
  message("✓ setup_cassia_env() with method='virtualenv' works")
}, error = function(e) {
  message("✗ setup_cassia_env() with method='virtualenv' failed: ", e$message)
})

tryCatch({
  # Test explicit conda method
  CASSIA::setup_cassia_env(conda_env = "cassia_env_conda", method = "conda")
  message("✓ setup_cassia_env() with method='conda' works")
}, error = function(e) {
  message("✗ setup_cassia_env() with method='conda' failed: ", e$message)
})

# Test 6: Options compatibility
message("\nTest 6: Options compatibility")
# Check that old options still work
old_env <- getOption("CASSIA.conda_env", default = NULL)
new_env <- getOption("CASSIA.env_name", default = NULL)
method <- getOption("CASSIA.env_method", default = NULL)

message("CASSIA.conda_env option: ", ifelse(is.null(old_env), "NULL", old_env))
message("CASSIA.env_name option: ", ifelse(is.null(new_env), "NULL", new_env))
message("CASSIA.env_method option: ", ifelse(is.null(method), "NULL", method))

message("\n--- Backward Compatibility Test Complete ---")