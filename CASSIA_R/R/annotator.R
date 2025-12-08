#' @import reticulate

# Use reticulate to import the CASSIA Python package
# All functions are exported at package level via __init__.py
py_cassia <- NULL

# =============================================================================
# Internal Diagnostic Functions for Setup Error Handling
# =============================================================================

# Get platform information for platform-specific error messages
.get_platform_info <- function() {
  list(
    os = .Platform$OS.type,
    system = Sys.info()["sysname"],
    is_windows = .Platform$OS.type == "windows",
    is_mac = Sys.info()["sysname"] == "Darwin",
    is_linux = Sys.info()["sysname"] == "Linux"
  )
}

# Check if Python is available on the system
.check_python_available <- function() {
  tryCatch({
    # First try py_discover_config (doesn't initialize Python)
    config <- reticulate::py_discover_config(required_module = NULL, use_environment = NULL)
    if (!is.null(config) && !is.null(config$python)) {
      return(list(available = TRUE, path = config$python, version = config$version))
    }

    # If py_discover_config fails, try to find Python in PATH
    # without initializing reticulate
    python_path <- Sys.which("python")
    if (python_path != "" && file.exists(python_path)) {
      # Get version by running python --version
      version_output <- tryCatch({
        system2(python_path, "--version", stdout = TRUE, stderr = TRUE)
      }, error = function(e) NULL)

      if (!is.null(version_output) && length(version_output) > 0) {
        # Parse version from output like "Python 3.10.11"
        version_match <- regmatches(version_output[1], regexpr("\\d+\\.\\d+", version_output[1]))
        if (length(version_match) > 0) {
          return(list(available = TRUE, path = python_path, version = version_match))
        }
      }
    }

    # Also try python3 on Unix-like systems
    python3_path <- Sys.which("python3")
    if (python3_path != "" && file.exists(python3_path)) {
      version_output <- tryCatch({
        system2(python3_path, "--version", stdout = TRUE, stderr = TRUE)
      }, error = function(e) NULL)

      if (!is.null(version_output) && length(version_output) > 0) {
        version_match <- regmatches(version_output[1], regexpr("\\d+\\.\\d+", version_output[1]))
        if (length(version_match) > 0) {
          return(list(available = TRUE, path = python3_path, version = version_match))
        }
      }
    }

    return(list(available = FALSE, path = NULL, version = NULL))
  }, error = function(e) {
    return(list(available = FALSE, path = NULL, version = NULL, error = e$message))
  })
}

# Check if Python version meets requirements (>= 3.9)
.check_python_version <- function(required_version = "3.9") {
  py_info <- .check_python_available()
  if (!py_info$available) {
    return(list(ok = FALSE, reason = "Python not found", current = NULL, required = required_version))
  }

  tryCatch({
    # Convert version to string if it's a package_version object
    version_str <- as.character(py_info$version)
    # Extract major.minor from version string (e.g., "3.10.1" -> c(3, 10))
    version_parts <- as.numeric(strsplit(version_str, "\\.")[[1]][1:2])
    required_parts <- as.numeric(strsplit(required_version, "\\.")[[1]][1:2])

    version_ok <- (version_parts[1] > required_parts[1]) ||
                  (version_parts[1] == required_parts[1] && version_parts[2] >= required_parts[2])

    return(list(ok = version_ok, current = version_str, required = required_version))
  }, error = function(e) {
    return(list(ok = FALSE, reason = e$message, current = as.character(py_info$version), required = required_version))
  })
}

# Check if environment managers (virtualenv/conda) are available
.check_env_managers_available <- function() {
  virtualenv_ok <- tryCatch({
    reticulate::virtualenv_list()
    TRUE
  }, error = function(e) FALSE)

  conda_ok <- tryCatch({
    conda_bin <- reticulate::conda_binary()
    !is.null(conda_bin) && file.exists(conda_bin)
  }, error = function(e) FALSE)

  list(
    virtualenv = virtualenv_ok,
    conda = conda_ok,
    any_available = virtualenv_ok || conda_ok
  )
}

# Error message templates with platform-specific fixes
.error_messages <- list(
  python_not_found = list(
    what = "Python is not installed or not found in your system PATH.",
    why = "CASSIA requires Python 3.9+ to run its cell annotation engine.",
    fix_windows = c(
      "1. Download Python from https://www.python.org/downloads/",
      "2. During installation, CHECK the box 'Add Python to PATH'",
      "3. Restart R/RStudio after installation",
      "4. Run: library(CASSIA)"
    ),
    fix_mac = c(
      "Option 1 - Using Homebrew (recommended):",
      "  Open Terminal and run: brew install python@3.10",
      "",
      "Option 2 - Download installer:",
      "  https://www.python.org/downloads/",
      "",
      "After installing, restart R/RStudio and run: library(CASSIA)"
    ),
    fix_linux = c(
      "Ubuntu/Debian:",
      "  sudo apt update && sudo apt install python3 python3-pip python3-venv",
      "",
      "CentOS/RHEL:",
      "  sudo yum install python3 python3-pip",
      "",
      "After installing, restart R/RStudio and run: library(CASSIA)"
    )
  ),

  python_version = list(
    what = "Your Python version is not compatible with CASSIA.",
    why = "CASSIA requires Python 3.9 or higher.",
    fix_windows = c(
      "1. Download Python 3.10+ from https://www.python.org/downloads/",
      "2. During installation, CHECK 'Add Python to PATH'",
      "3. You may need to uninstall the old Python version first",
      "4. Restart R/RStudio and run: library(CASSIA)"
    ),
    fix_mac = c(
      "1. Install Python 3.10+ using Homebrew:",
      "   brew install python@3.10",
      "",
      "2. Or download from https://www.python.org/downloads/",
      "3. Restart R/RStudio and run: library(CASSIA)"
    ),
    fix_linux = c(
      "Ubuntu/Debian:",
      "  sudo apt update && sudo apt install python3.10 python3.10-venv",
      "",
      "Or use pyenv to install a newer Python version.",
      "After installing, restart R/RStudio and run: library(CASSIA)"
    )
  ),

  no_env_manager = list(
    what = "Cannot create Python environment (no virtualenv or conda found).",
    why = "CASSIA needs an isolated Python environment to manage its dependencies.",
    fix_windows = c(
      "Option 1 - Install virtualenv (recommended):",
      "  Open Command Prompt and run: pip install virtualenv",
      "",
      "Option 2 - Install Miniconda:",
      "  Download from https://docs.conda.io/en/latest/miniconda.html",
      "",
      "After installing, restart R/RStudio and run: library(CASSIA)"
    ),
    fix_mac = c(
      "Option 1 - Install virtualenv:",
      "  pip3 install virtualenv",
      "",
      "Option 2 - Install Miniconda:",
      "  Download from https://docs.conda.io/en/latest/miniconda.html",
      "",
      "After installing, restart R/RStudio and run: library(CASSIA)"
    ),
    fix_linux = c(
      "Ubuntu/Debian:",
      "  sudo apt install python3-venv",
      "  # or: pip3 install virtualenv",
      "",
      "After installing, restart R/RStudio and run: library(CASSIA)"
    )
  ),

  env_setup_failed = list(
    what = "Failed to create the CASSIA Python environment.",
    why = "This may be due to permission issues or network problems.",
    fix_windows = c(
      "Try these solutions:",
      "",
      "1. Run R/RStudio as Administrator (one time only):",
      "   Right-click R/RStudio -> 'Run as administrator'",
      "   Then run: library(CASSIA)",
      "",
      "2. Check your internet connection (packages need to be downloaded)",
      "",
      "3. Try manually creating the environment:",
      "   Open Command Prompt and run:",
      "   python -m venv %USERPROFILE%\\.virtualenvs\\cassia_env"
    ),
    fix_mac = c(
      "Try these solutions:",
      "",
      "1. Check your internet connection",
      "",
      "2. Try manually creating the environment:",
      "   python3 -m venv ~/.virtualenvs/cassia_env",
      "",
      "3. If using conda:",
      "   conda create -n cassia_env python=3.10"
    ),
    fix_linux = c(
      "Try these solutions:",
      "",
      "1. Ensure python3-venv is installed:",
      "   sudo apt install python3-venv",
      "",
      "2. Check your internet connection",
      "",
      "3. Try manually creating the environment:",
      "   python3 -m venv ~/.virtualenvs/cassia_env"
    )
  ),

  package_install_failed = list(
    what = "Failed to install required Python packages.",
    why = "CASSIA needs packages like openai, pandas, and numpy to function.",
    fix_windows = c(
      "1. Check your internet connection",
      "",
      "2. If behind a corporate firewall, configure proxy in R:",
      "   Sys.setenv(http_proxy = 'http://your-proxy:port')",
      "   Sys.setenv(https_proxy = 'http://your-proxy:port')",
      "",
      "3. Try installing packages manually:",
      "   Open Command Prompt and run:",
      "   pip install openai pandas numpy anthropic requests matplotlib seaborn"
    ),
    fix_mac = c(
      "1. Check your internet connection",
      "",
      "2. Try installing packages manually:",
      "   pip3 install openai pandas numpy anthropic requests matplotlib seaborn"
    ),
    fix_linux = c(
      "1. Check your internet connection",
      "",
      "2. Try installing packages manually:",
      "   pip3 install openai pandas numpy anthropic requests matplotlib seaborn"
    )
  )
)

# Display formatted error message with platform-specific instructions
.show_setup_error <- function(error_type, details = NULL) {
  platform <- .get_platform_info()
  msg_template <- .error_messages[[error_type]]

  if (is.null(msg_template)) {
    warning("Unknown setup error occurred. Please check your Python installation.")
    return(invisible(NULL))
  }

  # Build the problem description
  problem_text <- msg_template$what
  if (!is.null(details)) {
    if (!is.null(details$current)) {
      problem_text <- paste0(problem_text, " (Found: Python ", details$current, ")")
    }
    if (!is.null(details$error)) {
      problem_text <- paste0(problem_text, "\nError: ", details$error)
    }
  }

  # Build formatted output
  lines <- c(
    "",
    "============================================================",
    "  CASSIA Setup Error",
    "============================================================",
    "",
    paste0("PROBLEM: ", problem_text),
    "",
    paste0("WHY THIS MATTERS: ", msg_template$why),
    "",
    "HOW TO FIX:",
    ""
  )

  # Add platform-specific instructions
  if (platform$is_windows) {
    lines <- c(lines, msg_template$fix_windows)
  } else if (platform$is_mac) {
    lines <- c(lines, msg_template$fix_mac)
  } else {
    lines <- c(lines, msg_template$fix_linux)
  }

  lines <- c(lines, "", "============================================================", "")

  # Print as a single message
  message(paste(lines, collapse = "\n"))
  invisible(NULL)
}

# =============================================================================
# Package Load Function
# =============================================================================

# Package-level variable to store load status for .onAttach messages
.cassia_load_status <- new.env(parent = emptyenv())

# Helper function to check if running during R CMD check
.is_cran_check <- function() {
  # Multiple ways to detect R CMD check environment
  if (!identical(Sys.getenv("NOT_CRAN"), "true")) {
    # Check for common R CMD check indicators
    if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) return(TRUE)
    if (nzchar(Sys.getenv("R_TESTS"))) return(TRUE)
    if (identical(Sys.getenv("R_CMD"), "true")) return(TRUE)
  }
  FALSE
}

.onLoad <- function(libname, pkgname) {

  # Initialize load status
  .cassia_load_status$message <- NULL
  .cassia_load_status$success <- FALSE
  .cassia_load_status$python_available <- FALSE

  # =========================================================================
  # Skip Python setup during R CMD check (CRAN policy)
  # =========================================================================
  if (.is_cran_check()) {
    .cassia_load_status$message <- "cran_check"
    return(invisible(NULL))
  }

  # =========================================================================
  # Step 1: Pre-flight checks - Verify Python is available and compatible
  # =========================================================================

  # Check if Python is available
  py_check <- .check_python_available()
  if (!py_check$available) {
    .cassia_load_status$message <- "python_not_found"
    return(invisible(NULL))
  }

  # Check Python version (must be >= 3.9)
  version_check <- .check_python_version("3.9")
  if (!version_check$ok) {
    .cassia_load_status$message <- "python_version"
    return(invisible(NULL))
  }

  # =========================================================================
  # Step 2: Set up or activate Python environment
  # =========================================================================

  # Get the environment name from the package configuration
  env_name <- getOption("CASSIA.env_name", default = NULL)
  if (is.null(env_name)) {
    env_name <- getOption("CASSIA.conda_env", default = "cassia_env")
  }

  # Set up the Python environment
  tryCatch({
    # Check if environment exists in either virtualenv or conda
    env_exists <- .check_env_exists(env_name)

    if (!env_exists$virtualenv && !env_exists$conda) {
      # Check if environment managers are available before trying to create
      env_managers <- .check_env_managers_available()
      if (!env_managers$any_available) {
        .cassia_load_status$message <- "no_env_manager"
        return(invisible(NULL))
      }

      # Environment doesn't exist, create it using the new setup function
      .cassia_load_status$message <- "env_setup"
      setup_cassia_env(conda_env = env_name)
    } else {
      # Environment exists, activate it using the appropriate method
      if (env_exists$virtualenv) {
        reticulate::use_virtualenv(env_name, required = TRUE)
        options(CASSIA.env_name = env_name)
        options(CASSIA.env_method = "virtualenv")
      } else if (env_exists$conda) {
        reticulate::use_condaenv(env_name, required = TRUE)
        options(CASSIA.env_name = env_name)
        options(CASSIA.env_method = "conda")
        # Maintain backward compatibility
        options(CASSIA.conda_env = env_name)
      }
    }

    # =========================================================================
    # Step 3: Import CASSIA Python package
    # =========================================================================
    # All functions are exported at the package level via CASSIA/__init__.py
    py_cassia <<- reticulate::import_from_path("CASSIA", path = system.file("python", package = "CASSIA"))

    # Success!
    .cassia_load_status$success <- TRUE
    .cassia_load_status$python_available <- TRUE

  }, error = function(e) {
    # If setup fails, try to run setup_cassia_env() automatically
    .cassia_load_status$message <- "retry_setup"
    tryCatch({
      setup_cassia_env(conda_env = env_name)

      # Try to import CASSIA Python package again after successful setup
      py_cassia <<- reticulate::import_from_path("CASSIA", path = system.file("python", package = "CASSIA"))

      # Success after retry!
      .cassia_load_status$success <- TRUE
      .cassia_load_status$python_available <- TRUE

    }, error = function(e2) {
      # Store error but don't show it - package should still load
      .cassia_load_status$message <- "setup_failed"
      .cassia_load_status$error <- e2$message
    })
  })
}

.onAttach <- function(libname, pkgname) {
  # Skip messages during R CMD check
  if (.is_cran_check()) {
    return(invisible(NULL))
  }

  # Display startup messages based on load status
  if (.cassia_load_status$success) {
    packageStartupMessage("CASSIA loaded successfully! Happy annotating!")
  } else if (!is.null(.cassia_load_status$message)) {
    msg <- switch(.cassia_load_status$message,
      "python_not_found" = "Note: Python not found. Run setup_cassia_env() to configure.",
      "python_version" = "Note: Python >= 3.9 required. Run setup_cassia_env() to configure.",
      "no_env_manager" = "Note: No virtualenv or conda found. Please install one first.",
      "env_setup" = "Setting up CASSIA Python environment...",
      "setup_failed" = "Note: Python setup incomplete. Run setup_cassia_env() to configure.",
      NULL
    )
    if (!is.null(msg)) {
      packageStartupMessage(msg)
    }
  }
}



# Helper function to try virtualenv setup
.try_virtualenv_setup <- function(env_name, python_version, pip_packages) {
  tryCatch({
    # Check if virtualenv already exists
    existing_envs <- tryCatch(reticulate::virtualenv_list(), error = function(e) character(0))
    
    if (!env_name %in% existing_envs) {
      # Create virtualenv
      reticulate::virtualenv_create(envname = env_name, python_version = python_version)
    }
    
    # Use the virtualenv
    reticulate::use_virtualenv(env_name, required = TRUE)
    
    # Install packages
    reticulate::virtualenv_install(envname = env_name, packages = pip_packages)
    
    # Set options to track the environment and method
    options(CASSIA.env_name = env_name)
    options(CASSIA.env_method = "virtualenv")
    
    message("Successfully set up CASSIA environment using virtualenv: ", env_name)
    return(TRUE)
    
  }, error = function(e) {
    warning("Virtualenv setup failed: ", e$message)
    return(FALSE)
  })
}

# Helper function to try conda setup (current logic)
.try_conda_setup <- function(conda_env, python_version, pip_packages) {
  tryCatch({
    # Check if conda environment exists
    existing_envs <- tryCatch(reticulate::conda_list()$name, error = function(e) character(0))
    
    if (!conda_env %in% existing_envs) {
      # Create conda environment
      reticulate::conda_create(envname = conda_env, python_version = python_version)
    }
    
    # Use the conda environment
    reticulate::use_condaenv(conda_env, required = TRUE)
    
    # Install packages
    reticulate::py_install(pip_packages, pip = TRUE)
    
    # Set options to track the environment and method
    options(CASSIA.env_name = conda_env)
    options(CASSIA.env_method = "conda")
    
    message("Successfully set up CASSIA environment using conda: ", conda_env)
    return(TRUE)
    
  }, error = function(e) {
    warning("Conda setup failed: ", e$message)
    return(FALSE)
  })
}

# Helper function to get environment name
.get_env_name <- function(env_name) {
  if (is.null(env_name)) {
    # Check for existing environment name from options
    existing_env <- getOption("CASSIA.env_name", default = NULL)
    if (!is.null(existing_env)) {
      return(existing_env)
    }
    # Fall back to legacy conda env option
    return(getOption("CASSIA.conda_env", default = "cassia_env"))
  }
  return(env_name)
}

# Helper function to check if environment exists
.check_env_exists <- function(env_name) {
  # Check virtualenv
  virtualenv_exists <- tryCatch({
    existing_virtualenvs <- reticulate::virtualenv_list()
    env_name %in% existing_virtualenvs
  }, error = function(e) FALSE)
  
  # Check conda
  conda_exists <- tryCatch({
    existing_condas <- reticulate::conda_list()$name
    env_name %in% existing_condas
  }, error = function(e) FALSE)
  
  return(list(virtualenv = virtualenv_exists, conda = conda_exists))
}

#' Set up CASSIA Python Environment
#'
#' This function sets up the required Python environment for CASSIA.
#' It can be used to create a new environment or update an existing one.
#' By default, it tries virtualenv first (simpler, more reliable), then falls back to conda.
#' 
#' @param conda_env The name of the environment to use. If NULL, uses the default from package configuration.
#' @param python_version The Python version to use. Default is "3.10".
#' @param pip_packages A character vector of pip packages to install.
#' @param method The method to use for environment setup. Options: "auto" (try virtualenv first, then conda), "virtualenv", "conda". Default is "auto".
#'
#' @return Invisible NULL. Called for side effects.
#' @export
setup_cassia_env <- function(conda_env = NULL, python_version = "3.10",
                           pip_packages = c("openai", "pandas", "numpy",
                                          "requests", "anthropic", "matplotlib", "seaborn"),
                           method = "auto") {
  
  # Get environment name
  env_name <- .get_env_name(conda_env)
  
  # Validate method parameter
  if (!method %in% c("auto", "virtualenv", "conda")) {
    stop("Method must be one of: 'auto', 'virtualenv', 'conda'")
  }
  
  # Check if environment already exists
  env_exists <- .check_env_exists(env_name)
  
  success <- FALSE
  
  if (method == "auto") {
    # Try virtualenv first, then conda
    if (env_exists$virtualenv) {
      message("Using existing virtualenv: ", env_name)
      success <- .try_virtualenv_setup(env_name, python_version, pip_packages)
    } else if (env_exists$conda) {
      message("Using existing conda environment: ", env_name)
      success <- .try_conda_setup(env_name, python_version, pip_packages)
    } else {
      # Neither exists, try to create virtualenv first
      message("Creating new environment. Trying virtualenv first...")
      success <- .try_virtualenv_setup(env_name, python_version, pip_packages)
      
      if (!success) {
        message("Virtualenv failed. Trying conda as fallback...")
        success <- .try_conda_setup(env_name, python_version, pip_packages)
      }
    }
  } else if (method == "virtualenv") {
    success <- .try_virtualenv_setup(env_name, python_version, pip_packages)
  } else if (method == "conda") {
    success <- .try_conda_setup(env_name, python_version, pip_packages)
  }
  
  if (!success) {
    stop("Failed to set up CASSIA Python environment with method: ", method)
  }
  
  # Maintain backward compatibility - keep the old conda_env option
  if (getOption("CASSIA.env_method", default = "conda") == "conda") {
    options(CASSIA.conda_env = env_name)
  }
  
  invisible(NULL)
}


#' Set API Key for LLM Provider
#'
#' @param api_key Character string containing the API key
#' @param provider Character string specifying the provider ('openai', 'anthropic', or 'openrouter')
#' @param persist Logical indicating whether to save the key to .Renviron (default: FALSE)
#' @return Invisible NULL. Called for side effects.
#' @export
setLLMApiKey <- function(api_key, provider = "anthropic", persist = FALSE) {
  # Set environment variable based on provider
  if (provider == "openai") {
    Sys.setenv(OPENAI_API_KEY = api_key)
    env_var_name <- "OPENAI_API_KEY"
  } else if (provider == "anthropic") {
    Sys.setenv(ANTHROPIC_API_KEY = api_key)
    env_var_name <- "ANTHROPIC_API_KEY"
  } else if (provider == "openrouter") {
    Sys.setenv(OPENROUTER_API_KEY = api_key)
    env_var_name <- "OPENROUTER_API_KEY"
  } else if (startsWith(provider, "http")) {
    # Handle custom HTTP endpoints
    Sys.setenv(CUSTOMIZED_API_KEY = api_key)
    env_var_name <- "CUSTOMIZED_API_KEY"
  } else {
    stop("Unsupported provider. Use 'openai', 'anthropic', 'openrouter', or an HTTP URL for custom endpoints")
  }
  
  # Set API key in Python CASSIA if available
  if (!is.null(py_cassia)) {
    py_cassia$set_api_key(api_key, provider)
  }
  
  # Persist to .Renviron if requested
  if (persist) {
    renviron_path <- path.expand("~/.Renviron")
    
    # Read existing content
    if (file.exists(renviron_path)) {
      lines <- readLines(renviron_path)
      # Remove existing entry if present
      lines <- lines[!grepl(paste0("^", env_var_name, "="), lines)]
    } else {
      lines <- character()
    }
    
    # Add new entry
    lines <- c(lines, paste0(env_var_name, "='", api_key, "'"))
    
    # Write back to .Renviron
    writeLines(lines, renviron_path)
    
    # Reload .Renviron
    readRenviron("~/.Renviron")
    
    message("API key has been saved to .Renviron and will be loaded in future sessions")
  }
  
  invisible(NULL)
}


#' Set OpenAI API Key
#'
#' @param api_key Character string containing the OpenAI API key
#' @param persist Logical indicating whether to save the key to .Renviron (default: FALSE)
#' @return Invisible NULL. Called for side effects.
#' @export
set_openai_api_key <- function(api_key, persist = FALSE) {
  setLLMApiKey(api_key, provider = "openai", persist = persist)
}


#' Set Anthropic API Key
#'
#' @param api_key Character string containing the Anthropic API key
#' @param persist Logical indicating whether to save the key to .Renviron (default: FALSE)
#' @return Invisible NULL. Called for side effects.
#' @export
setAnthropicApiKey <- function(api_key, persist = FALSE) {
  setLLMApiKey(api_key, provider = "anthropic", persist = persist)
}


#' Set OpenRouter API Key
#'
#' @param api_key Character string containing the OpenRouter API key
#' @param persist Logical indicating whether to save the key to .Renviron (default: FALSE)
#' @return Invisible NULL. Called for side effects.
#' @export
setOpenRouterApiKey <- function(api_key, persist = FALSE) {
  setLLMApiKey(api_key, provider = "openrouter", persist = persist)
}






#' Run Cell Type Analysis
#'
#' @param model Character string specifying the model to use.
#' @param temperature Numeric value for temperature parameter.
#' @param marker_list List of marker genes. Can be a character vector, comma-separated string, or
#'   a data frame with a 'gene' column.
#' @param tissue Character string specifying the tissue type.
#' @param species Character string specifying the species.
#' @param additional_info Additional information as a character string.
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter', default='openai')
#'
#' @param validator_involvement Validator involvement level: "v0" for high involvement (stronger validation), "v1" for moderate involvement (default: "v1")
#' @param use_reference Logical. Whether to use reference-based annotation for complex cases (default: FALSE)
#' @param reasoning Reasoning effort level: "high", "medium", or "low". Default: NULL (no extended reasoning)
#' @return A list containing three elements: structured_output, conversation_history, and reference_info.
#' @export
runCASSIA <- function(model = "google/gemini-2.5-flash-preview", temperature, marker_list, tissue, species, additional_info = NULL, provider = "openrouter", validator_involvement = "v1", use_reference = FALSE, reasoning = NULL) {
  # Convert marker_list to character vector if it's a data frame
  if (is.data.frame(marker_list)) {
    # Try common column names for gene markers
    if ("gene" %in% names(marker_list)) {
      marker_list <- as.character(marker_list$gene)
    } else if ("Gene" %in% names(marker_list)) {
      marker_list <- as.character(marker_list$Gene)
    } else if ("marker" %in% names(marker_list)) {
      marker_list <- as.character(marker_list$marker)
    } else if ("Marker" %in% names(marker_list)) {
      marker_list <- as.character(marker_list$Marker)
    } else {
      # Use first column as markers
      marker_list <- as.character(marker_list[[1]])
    }
  }

  tryCatch({
    result <- py_cassia$runCASSIA(
      model = model,
      temperature = temperature,
      marker_list = marker_list,
      tissue = tissue,
      species = species,
      additional_info = additional_info,
      provider = provider,
      validator_involvement = validator_involvement,
      use_reference = use_reference,
      reasoning = reasoning
    )

    # Convert structured_output (result[[1]])
    structured_output <- as.list(result[[1]])

    # Convert conversation_history (result[[2]])
    conversation_history <- lapply(result[[2]], function(entry) {
      if (is.list(entry) && length(entry) == 2) {
        list(agent = as.character(entry[[1]]), message = as.character(entry[[2]]))
      } else {
        warning("Unexpected entry structure in conversation_history")
        list(agent = "Unknown", message = "Conversion failed")
      }
    })

    # Convert reference_info (result[[3]])
    reference_info <- as.list(result[[3]])

    # Return flattened result with main fields accessible directly
    # Also include nested structured_output for backward compatibility
    return(list(
      main_cell_type = structured_output$main_cell_type,
      sub_cell_types = structured_output$sub_cell_types,
      possible_mixed_cell_types = structured_output$possible_mixed_cell_types,
      num_markers = structured_output$num_markers,
      iterations = structured_output$iterations,
      structured_output = structured_output,
      conversation_history = conversation_history,
      reference_info = reference_info
    ))
  }, error = function(e) {
    error_msg <- paste("Error in runCASSIA:", e$message, "\n",
                       "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}

#' Run Cell Type Analysis Multiple Times
#'
#' @param n Number of times to run the analysis.
#' @param tissue Character string specifying the tissue type.
#' @param species Character string specifying the species.
#' @param additional_info Additional information as a character string.
#' @param temperature Numeric value for temperature parameter. Default: 0.3 (updated to match Python).
#' @param marker_list List of marker genes.
#' @param model Character string specifying the model to use.
#' @param max_workers Maximum number of workers for parallel processing.
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter')
#' @param validator_involvement Validator involvement level: "v0" for high involvement (stronger validation), "v1" for moderate involvement (default: "v1")
#' @param use_reference Logical. Whether to use reference-based annotation for complex cases (default: FALSE)
#' @param reasoning Reasoning effort level: "high", "medium", or "low". Default: NULL (no extended reasoning)
#'
#' @return A list containing results from multiple runs, each with analysis_result, conversation_history, and reference_info.
#' @export
runCASSIA_n_times <- function(n, tissue, species, additional_info, temperature = 0.3, marker_list,
                           model = "google/gemini-2.5-flash-preview", max_workers = 10, provider = "openrouter", validator_involvement = "v1", use_reference = FALSE, reasoning = NULL) {
  tryCatch({
    result <- py_cassia$runCASSIA_n_times(
      n = as.integer(n),
      tissue = tissue,
      species = species,
      additional_info = additional_info,
      temperature = as.numeric(temperature),
      marker_list = marker_list,
      model = model,
      max_workers = as.integer(max_workers),
      provider = provider,
      validator_involvement = validator_involvement,
      use_reference = use_reference,
      reasoning = reasoning
    )

    # Convert the result to an R list
    converted_result <- lapply(seq_len(n), function(i) {
      run_result <- result[[as.character(i-1)]]  # Python uses 0-based indexing
      list(
        analysis_result = list(
          main_cell_type = run_result[[1]]$main_cell_type,
          sub_cell_types = run_result[[1]]$sub_cell_types,
          possible_mixed_cell_types = run_result[[1]]$possible_mixed_cell_types,
          num_markers = run_result[[1]]$num_markers,
          iterations = run_result[[1]]$iterations
        ),
        conversation_history = lapply(run_result[[2]], function(entry) {
          list(agent = entry[[1]], message = entry[[2]])
        }),
        reference_info = as.list(run_result[[3]])
      )
    })

    names(converted_result) <- as.character(seq_len(n) - 1)  # Match Python's 0-based indexing
    return(converted_result)
  }, error = function(e) {
    error_msg <- paste("Error in runCASSIA_n_times:", e$message, "\n",
                       "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}




#' Process Cell Type Analysis Results
#'
#' @param tissue Character string specifying the tissue type.
#' @param species Character string specifying the species.
#' @param additional_info Additional information as a character string.
#' @param temperature Numeric value for temperature parameter.
#' @param marker_list List of marker genes.
#' @param model Character string specifying the model to use.
#' @param max_workers Maximum number of workers for parallel processing.
#' @param n Number of times to run the analysis.
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter')
#' @param validator_involvement Validator involvement level: "v0" for high involvement (stronger validation), "v1" for moderate involvement (default: "v1")
#' @param use_reference Logical. Whether to use reference-based annotation for complex cases (default: FALSE)
#' @param generate_report Logical. Whether to generate an HTML report (default: TRUE)
#' @param report_output_path Character string. Path to save the HTML report (default: 'uq_report.html')
#'
#' @return A list containing processed results including variance analysis.
#' @export
runCASSIA_n_times_similarity_score <- function(tissue, species, additional_info, temperature, marker_list, model = "google/gemini-2.5-flash-preview", max_workers, n, provider = "openrouter", validator_involvement = "v1", use_reference = FALSE, generate_report = TRUE, report_output_path = NULL) {
  tryCatch({
    # Call the Python function with the new parameter structure
    processed_results <- py_cassia$runCASSIA_n_times_similarity_score(
      tissue = tissue,
      species = species,
      additional_info = additional_info,
      temperature = as.numeric(temperature),
      marker_list = marker_list,
      model = model,
      max_workers = as.integer(max_workers),
      n = as.integer(n),
      provider = provider,
      validator_involvement = validator_involvement,
      use_reference = use_reference,
      generate_report = generate_report,
      report_output_path = report_output_path
    )
    
    # Convert the processed results back to R
    r_results <- list(
      unified_results = as.character(processed_results$unified_results),
      consensus_types = list(
        general = as.character(processed_results$consensus_types[[1]]),
        sub = as.character(processed_results$consensus_types[[2]])
      ),
      general_celltype_llm = as.character(processed_results$general_celltype_llm),
      sub_celltype_llm = as.character(processed_results$sub_celltype_llm),
      Possible_mixed_celltypes_llm = as.list(processed_results$Possible_mixed_celltypes_llm),
      similarity_score = as.numeric(processed_results$similarity_score),
      original_results = lapply(processed_results$original_results, function(x) as.character(x)),
      llm_response = as.character(processed_results$llm_response)
    )
    
    return(r_results)
  }, error = function(e) {
    error_msg <- paste("Error in runCASSIA_n_times_similarity_score:", e$message, "\n",
                       "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}







#' Run Cell Type Analysis Batch
#'
#' @param marker A data frame or path to the CSV file containing marker data.
#' @param output_name Name of the output JSON file.
#' @param model Character string specifying the model to use.
#' @param temperature Numeric value for temperature parameter.
#' @param tissue Character string specifying the tissue type.
#' @param species Character string specifying the species.
#' @param additional_info Additional information as a character string.
#' @param celltype_column Name of the column containing cell types.
#' @param gene_column_name Name of the column containing gene names.
#' @param max_workers Maximum number of workers for parallel processing.
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter')
#' @param n_genes Number of top genes to use (default: 50)
#' @param max_retries Maximum number of retries for failed analyses (default: 1)
#' @param validator_involvement Validator involvement level: "v0" for high involvement (stronger validation), "v1" for moderate involvement (default: "v1")
#' @param ranking_method Method used to rank marker genes: "avg_log2FC" (default), "p_val_adj", "pct_diff", or "Score".
#' @param ascending Logical value indicating sort direction. If NULL (default), uses method-specific default.
#' @param reasoning Reasoning effort level: "high", "medium", or "low". Default: NULL (no extended reasoning)
#'
#' @return None. This function creates output files and prints execution time.
#' @export
runCASSIA_batch <- function(marker, output_name = "cell_type_analysis_results.json",
                          model = "google/gemini-2.5-flash-preview", temperature = 0, tissue = "lung",
                          species = "human", additional_info = NULL,
                          celltype_column = NULL, gene_column_name = NULL,
                          max_workers = 10, provider = "openrouter", n_genes = 50,
                          max_retries = 1, validator_involvement = "v1",
                          ranking_method = "avg_log2FC", ascending = NULL, reasoning = NULL) {
  execution_time <- system.time({
    # Convert R dataframe to Python if df_input is a dataframe
if (is.data.frame(marker)) {
  # Determine the cluster column to check
  cluster_col <- if (!is.null(celltype_column)) celltype_column else "cluster"

  # If the cluster column exists and is a factor, convert it to character
  if (cluster_col %in% names(marker) && is.factor(marker[[cluster_col]])) {
    marker[[cluster_col]] <- as.character(marker[[cluster_col]])
  }
  pd <- reticulate::import("pandas")
  # Direct conversion with convert=TRUE
  marker <- reticulate::r_to_py(marker, convert = TRUE)
} else if (!is.character(marker)) {
  stop("marker must be either a data frame or a character vector")
}
    
    py_cassia$runCASSIA_batch(
      marker = marker,  # Changed parameter name to match Python function
      output_name = output_name,
      model = model,
      temperature = temperature,
      tissue = tissue,
      species = species,
      additional_info = additional_info,
      celltype_column = celltype_column,
      gene_column_name = gene_column_name,
      max_workers = as.integer(max_workers),
      provider = provider,
      n_genes = as.integer(n_genes),
      max_retries = as.integer(max_retries),
      validator_involvement = validator_involvement,
      ranking_method = ranking_method,
      ascending = ascending,
      reasoning = reasoning
    )
  })
  
  print(paste("Execution time for runCASSIA_batch:"))
  print(execution_time)
}





#' Run Batch Analysis Multiple Times
#'
#' @param n Number of times to run the batch analysis.
#' @param marker Path to the CSV file containing marker data.
#' @param output_name Prefix for output JSON files.
#' @param model Character string specifying the model to use.
#' @param temperature Numeric value for temperature parameter.
#' @param tissue Character string specifying the tissue type.
#' @param species Character string specifying the species.
#' @param additional_info Additional information as a character string.
#' @param celltype_column Name of the column containing cell types.
#' @param gene_column_name Name of the column containing gene names.
#' @param max_workers Maximum number of workers for parallel processing.
#' @param batch_max_workers Maximum number of workers for batch processing.
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter')
#' @param max_retries Maximum number of retries for failed analyses (default: 1)
#' @param validator_involvement Validator involvement level: "v0" for high involvement (stronger validation), "v1" for moderate involvement (default: "v1")
#'
#' @return None. This function creates output files and prints execution time.
#' @export
runCASSIA_batch_n_times <- function(n, marker, output_name = "cell_type_analysis_results", 
                                  model = "google/gemini-2.5-flash-preview", temperature = 0, tissue = "lung", 
                                  species = "human", additional_info = NULL, 
                                  celltype_column = NULL, gene_column_name = NULL, 
                                  max_workers = 10, batch_max_workers = 5, 
                                  provider = "openrouter", max_retries = 1, validator_involvement = "v1") {

  if (is.data.frame(marker)) {
    # Determine the cluster column to check
    cluster_col <- if (!is.null(celltype_column)) celltype_column else "cluster"

    # If the cluster column exists and is a factor, convert it to character
    if (cluster_col %in% names(marker) && is.factor(marker[[cluster_col]])) {
      marker[[cluster_col]] <- as.character(marker[[cluster_col]])
    }
    pd <- reticulate::import("pandas")
    # Direct conversion with convert=TRUE
    marker <- reticulate::r_to_py(marker, convert = TRUE)
  } else if (!is.character(marker)) {
    stop("marker must be either a data frame or a character vector")
  }

  execution_time <- system.time({
    tryCatch({
      py_cassia$runCASSIA_batch_n_times(
        as.integer(n), marker, output_name, model, temperature, tissue, 
        species, additional_info, celltype_column, gene_column_name, 
        as.integer(max_workers), as.integer(batch_max_workers), provider,
        as.integer(max_retries), validator_involvement
      )
    }, error = function(e) {
      stop(paste("Error in runCASSIA_batch_n_times:", e$message))
    })
  })
  
  print(paste("Execution time for runCASSIA_batch_n_times:"))
  print(execution_time)
}




#' Process and Save Batch Results
#'
#' @param marker Path to the marker file.
#' @param file_pattern Pattern to match result files.
#' @param output_name Name of the output CSV file.
#' @param celltype_column Name of the column containing cell types.
#' @param max_workers Maximum number of workers for parallel processing.
#' @param model Model to use for processing (default: "gpt-4o")
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter')
#' @param main_weight Weight for the main cell type.
#' @param sub_weight Weight for the sub cell type.
#' @param generate_report Logical. Whether to generate an HTML report (default: TRUE)
#' @param report_output_path Character string. Path to save the HTML report (default: 'uq_batch_report.html')
#'
#' @return None. This function processes and saves results to a CSV file and prints execution time.
#' @export
runCASSIA_similarity_score_batch <- function(marker, file_pattern, output_name,
                                               celltype_column = NULL, max_workers = 10, model = "google/gemini-2.5-flash-preview", provider = "openrouter", main_weight=0.5, sub_weight=0.5, generate_report = TRUE, report_output_path = NULL) {


  if (is.data.frame(marker)) {
    # Determine the cluster column to check
    cluster_col <- if (!is.null(celltype_column)) celltype_column else "cluster"

    # If the cluster column exists and is a factor, convert it to character
    if (cluster_col %in% names(marker) && is.factor(marker[[cluster_col]])) {
      marker[[cluster_col]] <- as.character(marker[[cluster_col]])
    }
    pd <- reticulate::import("pandas")
    # Direct conversion with convert=TRUE
    marker <- reticulate::r_to_py(marker, convert = TRUE)
  } else if (!is.character(marker)) {
    stop("marker must be either a data frame or a character vector")
  }

  execution_time <- system.time({
    tryCatch({
      py_cassia$runCASSIA_similarity_score_batch(
        marker = marker,
        file_pattern = file_pattern,
        output_name = output_name,
        celltype_column = celltype_column,
        max_workers = as.integer(max_workers),
        model = model,
        provider = provider,
        main_weight = main_weight,
        sub_weight = sub_weight,
        generate_report = generate_report,
        report_output_path = report_output_path
      )
    }, error = function(e) {
      stop(paste("Error in runCASSIA_similarity_score_batch:", e$message))
    })
  })

  print(paste("Execution time for runCASSIA_batch_get_similarity_score:"))
  print(execution_time)
}







#' Generate Cell Type Analysis Report
#'
#' @param full_result_path Path to the full results CSV file
#' @param marker Path to the marker genes CSV file or data frame with marker data
#' @param cluster_name Name of the cluster to analyze
#' @param major_cluster_info General information about the dataset (e.g., "Human PBMC")
#' @param output_name Name of the output HTML file
#' @param num_iterations Number of iterations for marker analysis (default=5)
#' @param model Model to use for analysis (default="google/gemini-2.5-flash-preview")
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter')
#' @param temperature Sampling temperature (0-1)
#' @param conversation_history_mode Mode for extracting conversation history ("full", "final", or "none")
#' @param search_strategy Search strategy - "breadth" (test multiple hypotheses) or "depth" (one hypothesis at a time) (default: "breadth")
#' @param report_style Style of report ("per_iteration" or "total_summary") (default: "per_iteration")
#' @param validator_involvement Validator involvement level: "v0" for high involvement (stronger validation), "v1" for moderate involvement (default: "v1")
#' @param ... Additional arguments passed to the Python function
#'
#' @return None. This function generates output files.
#' @export
runCASSIA_annotationboost <- function(full_result_path, 
                                     marker, 
                                     cluster_name, 
                                     major_cluster_info, 
                                     output_name, 
                                     num_iterations = 5,
                                     model = "google/gemini-2.5-flash-preview",
                                     provider = "openrouter",
                                     temperature = 0,
                                     conversation_history_mode = "final",
                                     search_strategy = "breadth",
                                     report_style = "per_iteration",
                                     validator_involvement = "v1",
                                     ...) {

  if (is.data.frame(marker)) {
    # Factors can cause issues with reticulate, convert to character
    if ("cluster" %in% names(marker) && is.factor(marker[["cluster"]])) {
      marker[["cluster"]] <- as.character(marker[["cluster"]])
    }
    pd <- reticulate::import("pandas")
    # Direct conversion with convert=TRUE
    marker <- reticulate::r_to_py(marker, convert = TRUE)
  } else if (!is.character(marker)) {
    stop("marker must be either a data frame or a character vector")
  }
                                                      
  tryCatch({
    py_cassia$runCASSIA_annotationboost(
      full_result_path = full_result_path,
      marker = marker,
      cluster_name = cluster_name,
      major_cluster_info = major_cluster_info,
      output_name = output_name,
      num_iterations = as.integer(num_iterations),
      model = model,
      provider = provider,
      temperature = as.numeric(temperature),
      conversation_history_mode = conversation_history_mode,
      search_strategy = search_strategy,
      report_style = report_style
    )

    invisible(NULL)

  }, error = function(e) {
    error_msg <- paste("Error in runCASSIA_annotationboost:", e$message)
    stop(error_msg)
  })
}




#' Generate Cell Type Analysis Report with Additional Task
#'
#' @param full_result_path Path to the full results CSV file
#' @param marker Path to the marker genes CSV file or data frame with marker data
#' @param cluster_name Name of the cluster to analyze
#' @param major_cluster_info General information about the dataset (e.g., "Human PBMC")
#' @param output_name Name of the output HTML file
#' @param num_iterations Number of iterations for marker analysis (default=5)
#' @param model Model to use for analysis (default="google/gemini-2.5-flash-preview")
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter')
#' @param additional_task Additional task to perform during analysis
#' @param temperature Sampling temperature (0-1)
#' @param conversation_history_mode Mode for extracting conversation history ("full", "final", or "none")
#' @param search_strategy Search strategy - "breadth" (test multiple hypotheses) or "depth" (one hypothesis at a time) (default: "breadth")
#' @param report_style Style of report ("per_iteration" or "total_summary") (default: "per_iteration")
#' @param validator_involvement Validator involvement level: "v0" for high involvement (stronger validation), "v1" for moderate involvement (default: "v1")
#' @param ... Additional parameters for future compatibility
#' 
#' @return None. This function generates output files.
#' @export
runCASSIA_annotationboost_additional_task <- function(full_result_path, 
                                                     marker, 
                                                     cluster_name, 
                                                     major_cluster_info, 
                                                     output_name, 
                                                     num_iterations = 5,
                                                     model = "google/gemini-2.5-flash-preview",
                                                     provider = "openrouter",
                                                     additional_task = "check if this is a cancer cluster",
                                                     temperature = 0,
                                                     conversation_history_mode = "final",
                                                     search_strategy = "breadth",
                                                     report_style = "per_iteration",
                                                     validator_involvement = "v1",
                                                     ...) {

  if (is.data.frame(marker)) {
    # Factors can cause issues with reticulate, convert to character
    if ("cluster" %in% names(marker) && is.factor(marker[["cluster"]])) {
      marker[["cluster"]] <- as.character(marker[["cluster"]])
    }
    pd <- reticulate::import("pandas")
    # Direct conversion with convert=TRUE
    marker <- reticulate::r_to_py(marker, convert = TRUE)
  } else if (!is.character(marker)) {
    stop("marker must be either a data frame or a character vector")
  }
                                                      
  tryCatch({
    py_cassia$runCASSIA_annotationboost_additional_task(
      full_result_path = full_result_path,
      marker = marker,
      cluster_name = cluster_name,
      major_cluster_info = major_cluster_info,
      output_name = output_name,
      num_iterations = as.integer(num_iterations),
      model = model,
      provider = provider,
      additional_task = additional_task,
      temperature = as.numeric(temperature),
      conversation_history_mode = conversation_history_mode,
      search_strategy = search_strategy,
      report_style = report_style,
      validator_involvement = validator_involvement,
      ...
    )
    
    invisible(NULL)
    
  }, error = function(e) {
    error_msg <- paste("Error in runCASSIA_annotationboost_additional_task:", e$message)
    stop(error_msg)
  })
}


#' Run Scoring with Progress Updates
#'
#' @param input_file Path to input CSV file
#' @param output_file Path to output CSV file (optional)
#' @param max_workers Maximum number of parallel workers
#' @param model Model to use
#' @param provider AI provider to use ('openai', 'anthropic', or 'openrouter')
#' @param max_retries Maximum number of retries for failed analyses (default: 1)
#'
#' @return None
#' @export
runCASSIA_score_batch <- function(input_file, 
                                    output_file = NULL, 
                                    max_workers = 4, 
                                    model = "deepseek/deepseek-chat-v3-0324",
                                    provider = "openrouter",
                                    max_retries = 1) {
  tryCatch({
    results <- py_cassia$runCASSIA_score_batch(
      input_file = input_file,
      output_file = output_file,
      max_workers = as.integer(max_workers),
      model = model,
      provider = provider,
      max_retries = as.integer(max_retries)
    )
    
  })
}

#' Generate HTML Reports from Scored Results
#'
#' This function processes a CSV file containing scored annotation results and generates 
#' individual HTML reports for each row, along with an index page.
#'
#' @param csv_path Character string. Path to the CSV file containing scored results.
#' @param output_name Character string. Base name for the output index file (default: "CASSIA_reports_summary").
#'
#' @details 
#' The function generates:
#' 1. Individual HTML reports for each annotation result
#' 2. An index page linking to all generated reports
#' 
#' The CSV file should contain columns:
#' - Conversation History
#' - Scoring_Reasoning
#' - Score
#'
#' @return None. Files are written to the current directory.
#' @export
#'
#' @examples
#' \dontrun{
#' runCASSIA_generate_score_report("path/to/scored_results.csv")
#' }
runCASSIA_generate_score_report <- function(csv_path, output_name = "CASSIA_reports_summary") {
  tryCatch({
    py_cassia$runCASSIA_generate_score_report(
      csv_path = csv_path,
      index_name = output_name
    )
    message("Reports generated successfully. Check the current directory for the generated HTML files.")
  }, error = function(e) {
    error_msg <- paste("Error in generating reports:", e$message, "\n",
                      "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}


#' Run Complete Cell Analysis Pipeline
#'
#' @param output_file_name Base name for output files
#' @param tissue Tissue type being analyzed
#' @param species Species being analyzed
#' @param marker Marker data (data frame or file path)
#' @param max_workers Maximum number of concurrent workers (default: 4)
#' @param annotation_model Model to use for initial annotation (default: "google/gemini-2.5-flash-preview")
#' @param annotation_provider Provider for initial annotation (default: "openrouter")
#' @param score_model Model to use for scoring (default: "deepseek/deepseek-chat-v3-0324")
#' @param score_provider Provider for scoring (default: "openrouter")
#' @param annotationboost_model Model to use for boosting low-scoring annotations (default: "google/gemini-2.5-flash-preview")
#' @param annotationboost_provider Provider for boosting low-scoring annotations (default: "openrouter")
#' @param score_threshold Threshold for identifying low-scoring clusters (default: 75)
#' @param additional_info Additional information for analysis (default: NULL)
#' @param max_retries Maximum number of retries for failed analyses (default: 1)
#' @param do_merge_annotations Whether to run the merging annotations step (default: TRUE)
#' @param merge_model Model to use for merging annotations (default: "deepseek/deepseek-chat-v3-0324")
#' @param merge_provider Provider to use for merging annotations (default: "openrouter")
#' @param conversation_history_mode Mode for extracting conversation history ("full", "final", or "none") (default: "final")
#' @param search_strategy Search strategy for annotation boost - "breadth" (test multiple hypotheses) or "depth" (one hypothesis at a time) (default: "breadth")
#' @param report_style Style of report for annotation boost ("per_iteration" or "total_summary") (default: "per_iteration")
#' @param ranking_method Method to rank genes ('avg_log2FC', 'p_val_adj', 'pct_diff', 'Score') (default: "avg_log2FC")
#' @param ascending Sort direction (NULL uses default for each method) (default: NULL)
#' @param validator_involvement Validator involvement level: "v0" for high involvement (stronger validation), "v1" for moderate involvement (default: "v1")
#' @param output_dir Directory where the output folder will be created. If NULL, uses current working directory. (default: NULL)
#'
#' @return None. Creates output files and generates reports.
#' @export
runCASSIA_pipeline <- function(
    output_file_name,
    tissue,
    species,
    marker,
    max_workers = 4,
    annotation_model = "google/gemini-2.5-flash",
    annotation_provider = "openrouter",
    score_model = "google/gemini-2.5-flash",
    score_provider = "openrouter",
    annotationboost_model = "google/gemini-2.5-flash",
    annotationboost_provider = "openrouter",
    score_threshold = 75,
    additional_info = NULL,
    max_retries = 1,
    do_merge_annotations = TRUE,
    merge_model = "google/gemini-2.5-flash",
    merge_provider = "openrouter",
    conversation_history_mode = "final",
    search_strategy = "breadth",
    report_style = "per_iteration",
    ranking_method = "avg_log2FC",
    ascending = NULL,
    validator_involvement = "v1",
    output_dir = NULL
) {
  # Convert R dataframe to Python if marker is a dataframe
  if (is.data.frame(marker)) {
    # Determine the cluster column to check
    cluster_col <- "cluster"
    if ("cluster" %in% names(marker)) {
      cluster_col <- "cluster"
    } else if (length(names(marker)) > 0) {
      cluster_col <- names(marker)[1]  # Use first column if no cluster column
    }

    # If the cluster column exists and is a factor, convert it to character
    if (cluster_col %in% names(marker) && is.factor(marker[[cluster_col]])) {
      marker[[cluster_col]] <- as.character(marker[[cluster_col]])
    }
    pd <- reticulate::import("pandas")
    # Direct conversion with convert=TRUE
    marker <- reticulate::r_to_py(marker, convert = TRUE)
  } else if (!is.character(marker)) {
    stop("marker must be either a data frame or a character vector")
  }

  tryCatch({
    py_cassia$runCASSIA_pipeline(
      output_file_name = output_file_name,
      tissue = tissue,
      species = species,
      marker = marker,
      max_workers = as.integer(max_workers),
      annotation_model = annotation_model,
      annotation_provider = annotation_provider,
      score_model = score_model,
      score_provider = score_provider,
      annotationboost_model = annotationboost_model,
      annotationboost_provider = annotationboost_provider,
      score_threshold = as.numeric(score_threshold),
      additional_info = if (is.null(additional_info)) "None" else additional_info,
      max_retries = as.integer(max_retries),
      merge_annotations = do_merge_annotations,
      merge_model = merge_model,
      merge_provider = merge_provider,
      conversation_history_mode = conversation_history_mode,
      ranking_method = ranking_method,
      ascending = ascending,
      report_style = report_style,
      validator_involvement = validator_involvement,
      output_dir = output_dir
    )
    
    invisible(NULL)
    
  }, error = function(e) {
    error_msg <- paste("Error in runCASSIA_pipeline:", e$message)
    if (!is.null(reticulate::py_last_error())) {
      error_msg <- paste(error_msg, "\nPython traceback:", reticulate::py_last_error())
    }
    stop(error_msg)
  })
}



#' Compare Cell Types Using Different Models
#'
#' @param tissue Character string specifying the tissue type
#' @param celltypes Character vector of 2-4 cell types to compare
#' @param marker Character vector or string of marker genes
#' @param species Character string specifying the species (default: "human")
#' @param model_list Character vector of model names to use (optional)
#' @param output_file Character string specifying output file path (optional)
#'
#' @return A list containing responses from different models
#' @export
compareCelltypes <- function(tissue, celltypes, marker, species = "human", model_list = NULL, output_file = NULL) {
  tryCatch({
    # Input validation
    if (length(celltypes) < 2 || length(celltypes) > 5) {
      stop("Please provide 2-5 cell types to compare")
    }
    
    # Convert marker_set to string if it's a vector
    if (is.vector(marker)) {
      marker <- paste(marker, collapse = ", ")
    }
    
    # Call the Python function
    responses <- py_cassia$compareCelltypes(
      tissue = tissue,
      celltypes = celltypes,
      marker_set = marker,
      species = species,
      model_list = model_list,
      output_file = output_file
    )

  }, error = function(e) {
    error_msg <- paste("Error in compareCelltypes:", e$message, "\n",
                      "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}


#' Symphony Compare - Advanced Multi-Model Cell Type Comparison with Consensus Building
#'
#' Orchestrate multiple AI models to compare cell types with automatic consensus building.
#' This function conducts a comprehensive cell type comparison using multiple AI models in parallel,
#' automatically triggering discussion rounds when models disagree on the best matching cell type.
#' Think of it as a virtual panel of expert biologists debating and reaching consensus.
#'
#' @param tissue Character string specifying the tissue type (e.g., "blood", "brain", "liver")
#' @param celltypes Character vector of 2-4 cell types to compare
#' @param marker_set Character vector or string of gene markers to analyze
#' @param species Character string specifying the species (default: "human")
#' @param model_preset Character string specifying preset model configuration. Options:
#'   \itemize{
#'     \item "symphony": High-performance ensemble (Claude, GPT-4, Gemini Pro)
#'     \item "quartet": Balanced 4-model ensemble
#'     \item "budget": Cost-effective models
#'     \item "custom": Use custom_models list
#'   }
#' @param custom_models Character vector of custom models to use (when model_preset="custom")
#' @param output_dir Character string specifying directory to save results (default: current directory)
#' @param output_basename Character string for base name of output files (auto-generated if NULL)
#' @param enable_discussion Logical indicating whether to enable automatic discussion rounds when no consensus (default: TRUE)
#' @param max_discussion_rounds Integer specifying maximum discussion rounds to perform (default: 2)
#' @param consensus_threshold Numeric value (0-1) specifying fraction of models that must agree for consensus (default: 0.8)
#' @param generate_report Logical indicating whether to generate interactive HTML report (default: TRUE)
#' @param api_key Character string for OpenRouter API key (uses environment variable if NULL)
#' @param verbose Logical indicating whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#'   \itemize{
#'     \item results: List of all model responses and scores
#'     \item consensus: The consensus cell type (if reached)
#'     \item confidence: Confidence level of the consensus (0-1)
#'     \item csv_file: Path to the generated CSV file
#'     \item html_file: Path to the generated HTML report (if enabled)
#'     \item summary: Summary statistics of the comparison
#'     \item dataframe: R data frame with structured results
#'   }
#'
#' @examples
#' \dontrun{
#' # Basic usage - let Symphony Compare handle everything
#' results <- symphonyCompare(
#'   tissue = "peripheral blood",
#'   celltypes = c("T cell", "B cell", "NK cell", "Monocyte"),
#'   marker_set = c("CD3", "CD4", "CD8", "CD19", "CD20", "CD16", "CD56", "CD14"),
#'   species = "human"
#' )
#' 
#' # Access the results
#' cat("Consensus:", results$consensus, "\n")
#' cat("Confidence:", sprintf("%.1f%%", results$confidence * 100), "\n")
#' 
#' # Advanced usage with custom settings
#' results <- symphonyCompare(
#'   tissue = "brain",
#'   celltypes = c("Neuron", "Astrocyte", "Microglia", "Oligodendrocyte"),
#'   marker_set = c("RBFOX3", "GFAP", "IBA1", "OLIG2", "MAP2", "S100B", "CD11B", "MBP"),
#'   species = "mouse",
#'   model_preset = "quartet",  # Use 4 models instead of 3
#'   enable_discussion = TRUE,  # Enable automatic discussion rounds
#'   max_discussion_rounds = 3,  # Allow up to 3 discussion rounds
#'   consensus_threshold = 0.75,  # 75% of models must agree
#'   output_dir = "./symphony_results",
#'   verbose = TRUE
#' )
#' }
#'
#' @export
symphonyCompare <- function(tissue, celltypes, marker_set, species = "human", 
                          model_preset = "symphony", custom_models = NULL,
                          output_dir = NULL, output_basename = NULL,
                          enable_discussion = TRUE, max_discussion_rounds = 2L,
                          consensus_threshold = 0.8, generate_report = TRUE,
                          api_key = NULL, verbose = TRUE) {
  tryCatch({
    # Input validation
    if (length(celltypes) < 2 || length(celltypes) > 4) {
      stop("Please provide 2-4 cell types to compare")
    }
    
    # Convert marker_set to string if it's a vector
    if (is.vector(marker_set)) {
      marker_set <- paste(marker_set, collapse = ", ")
    }
    
    # Convert R types to Python-compatible types
    max_discussion_rounds <- as.integer(max_discussion_rounds)
    consensus_threshold <- as.numeric(consensus_threshold)
    enable_discussion <- as.logical(enable_discussion)
    generate_report <- as.logical(generate_report)
    verbose <- as.logical(verbose)
    
    # Call the Python function
    results <- py_cassia$symphonyCompare(
      tissue = tissue,
      celltypes = celltypes,
      marker_set = marker_set,
      species = species,
      model_preset = model_preset,
      custom_models = custom_models,
      output_dir = output_dir,
      output_basename = output_basename,
      enable_discussion = enable_discussion,
      max_discussion_rounds = max_discussion_rounds,
      consensus_threshold = consensus_threshold,
      generate_report = generate_report,
      api_key = api_key,
      verbose = verbose
    )
    
    # Convert Python dataframe to R dataframe if it exists
    if (!is.null(results$dataframe)) {
      results$dataframe <- reticulate::py_to_r(results$dataframe)
    }
    
    # Return the results
    return(results)
    
  }, error = function(e) {
    error_msg <- paste("Error in symphonyCompare:", e$message, "\n",
                      "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}






#' Process Subclusters
#'
#' @param marker Marker data (data frame or file path)
#' @param major_cluster_info Description of the major cluster type
#' @param output_name Base name for the output file (will add .csv if not present)
#' @param model Model name for Claude API (default: "claude-3-5-sonnet-20241022")
#' @param temperature Temperature parameter for API calls (default: 0)
#' @param provider AI provider to use (default: "anthropic")
#' @param n_genes Number of top genes to use (default: 50)
#'
#' @return None. This function processes subclusters and saves results to a CSV file.
#' @export
runCASSIA_subclusters <- function(marker, major_cluster_info, output_name, 
                               model = "google/gemini-2.5-flash-preview", temperature = 0, 
                               provider = "openrouter", n_genes = 50L) {
  py_cassia$runCASSIA_subclusters(
    marker = marker,
    major_cluster_info = major_cluster_info,
    output_name = output_name,
    model = model,
    temperature = as.numeric(temperature),
    provider = provider,
    n_genes = as.integer(n_genes)
  )
}

#' Run Analysis Multiple Times for Subclusters
#'
#' @param n Number of times to run the analysis
#' @param marker Marker data (data frame or file path)
#' @param major_cluster_info Description of the major cluster type
#' @param base_output_name Base name for output CSV files
#' @param model Model name for Claude API (default: "claude-3-5-sonnet-20241022")
#' @param temperature Temperature parameter for API calls (default: 0)
#' @param provider AI provider to use (default: "anthropic")
#' @param max_workers Maximum number of workers for parallel processing (default: 5)
#' @param n_genes Number of top genes to use (default: 50)
#'
#' @return None. This function runs the analysis multiple times and saves results to CSV files.
#' @export
runCASSIA_n_subcluster <- function(n, marker, major_cluster_info, base_output_name, 
                                               model = "google/gemini-2.5-flash-preview", temperature = 0, 
                                               provider = "openrouter", max_workers = 5,n_genes=50L) {
  py_cassia$runCASSIA_n_subcluster(
    n = as.integer(n),
    marker = marker,
    major_cluster_info = major_cluster_info,
    base_output_name = base_output_name,
    model = model,
    temperature = temperature,
    provider = provider,
    max_workers = as.integer(max_workers),
    n_genes = as.integer(n_genes)
  )
}

#' Generate Subclustering HTML Report
#' 
#' Generate a beautiful HTML report for subclustering batch results, showing annotation, 
#' reasoning, and top markers in an interactive format.
#' 
#' @param csv_path Character string specifying the path to the subclustering CSV results file
#' @param html_report_path Character string specifying the output path for the HTML report. 
#'   If NULL, will use the same name as csv_path but with .html extension
#' @param model_name Character string specifying the model name to display in the report. 
#'   If NULL, will extract from the CSV filename
#' 
#' @return None. This function generates an HTML report file.
#' 
#' @examples
#' \dontrun{
#' # Generate a subclustering report from CSV results
#' generate_subclustering_report(
#'   csv_path = "subclustering_results.csv",
#'   html_report_path = "subclustering_report.html",
#'   model_name = "Gemini-2.5-Flash"
#' )
#' 
#' # Auto-generate HTML filename and model name
#' generate_subclustering_report("subclustering_results.csv")
#' }
#' 
#' @export
generate_subclustering_report <- function(csv_path, html_report_path = NULL, model_name = NULL) {
  tryCatch({
    py_cassia$generate_subclustering_report(
      csv_path = csv_path,
      html_report_path = html_report_path,
      model_name = model_name
    )
  }, error = function(e) {
    error_msg <- paste("Error in generate_subclustering_report:", e$message, "\n",
                      "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}

#' Generate HTML Evaluation Report
#' 
#' Generate a comprehensive HTML report for evaluation results with visualizations,
#' metrics, and sample analyses.
#' 
#' @param result_df Data frame containing the evaluation results
#' @param gold_col Character string specifying the column name for ground truth/gold standard values
#' @param pred_col Character string specifying the column name for predicted values
#' @param score_col Character string specifying the column name for evaluation scores (default: "score")
#' @param reasoning_col Character string specifying the column name for reasoning/explanation text (default: "reasoning")
#' @param metrics Named list of pre-calculated metrics (optional, will be calculated if NULL)
#' @param html_report_path Character string specifying the output path for HTML report (default: "report.html")
#' @param model_name Character string specifying the model name to display in the report (optional)
#' 
#' @return None. This function generates an HTML report file.
#' 
#' @examples
#' \dontrun{
#' # Generate evaluation report
#' generate_html_report(
#'   result_df = evaluation_results,
#'   gold_col = "Cluster ID",
#'   pred_col = "Predicted General Cell Type",
#'   score_col = "evaluation_score",
#'   reasoning_col = "explanation",
#'   html_report_path = "evaluation_report.html",
#'   model_name = "Claude-3.5-Sonnet"
#' )
#' }
#' 
#' @export
generate_html_report <- function(result_df, gold_col, pred_col, score_col = "score", 
                                reasoning_col = "reasoning", metrics = NULL, 
                                html_report_path = "report.html", model_name = NULL) {
  tryCatch({
    py_cassia$generate_html_report(
      result_df = result_df,
      gold_col = gold_col,
      pred_col = pred_col,
      score_col = score_col,
      reasoning_col = reasoning_col,
      metrics = metrics,
      html_report_path = html_report_path,
      model_name = model_name
    )
  }, error = function(e) {
    error_msg <- paste("Error in generate_html_report:", e$message, "\n",
                      "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}

#' Calculate Evaluation Metrics
#' 
#' Calculate comprehensive metrics from evaluation results including mean scores,
#' distributions, and performance ratios.
#' 
#' @param eval_df Data frame containing evaluation results
#' @param score_col Character string specifying the column name for evaluation scores (default: "score")
#' 
#' @return A named list containing calculated metrics:
#'   \itemize{
#'     \item mean_score: Average score
#'     \item median_score: Median score  
#'     \item min_score: Minimum score
#'     \item max_score: Maximum score
#'     \item std_score: Standard deviation of scores
#'     \item count: Number of evaluations
#'     \item For 0-5 scale: perfect_ratio, very_good_ratio, good_ratio, partial_ratio, poor_ratio, nonsensical_ratio
#'   }
#' 
#' @examples
#' \dontrun{
#' # Calculate metrics from evaluation results
#' metrics <- calculate_evaluation_metrics(
#'   eval_df = evaluation_results,
#'   score_col = "evaluation_score"
#' )
#' 
#' # Print key metrics
#' cat("Mean Score:", metrics$mean_score, "\n")
#' cat("Perfect Predictions:", sprintf("%.1f%%", metrics$perfect_ratio * 100), "\n")
#' }
#' 
#' @export
calculate_evaluation_metrics <- function(eval_df, score_col = "score") {
  tryCatch({
    result <- py_cassia$calculate_evaluation_metrics(
      eval_df = eval_df,
      score_col = score_col
    )
    # Convert Python dict to R list
    return(reticulate::py_to_r(result))
  }, error = function(e) {
    error_msg <- paste("Error in calculate_evaluation_metrics:", e$message, "\n",
                      "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}

# =============================================================================
# Merging Annotation Functions
# =============================================================================

#' Merge Cell Cluster Annotations
#'
#' Agent function that reads a CSV file with cell cluster annotations and merges/groups them
#' using an LLM to suggest biologically meaningful groupings.
#'
#' @param csv_path Path to the CSV file containing cluster annotations
#' @param output_path Path to save the results (if NULL, returns DataFrame without saving)
#' @param provider LLM provider to use ("openai", "anthropic", or "openrouter") (default: "openai")
#' @param model Specific model to use (if NULL, uses default for provider)
#' @param api_key API key for the provider (if NULL, gets from environment)
#' @param additional_context Optional domain-specific context to help with annotation
#' @param batch_size Number of clusters to process in each LLM call (default: 20)
#' @param detail_level Level of detail for groupings: "broad", "detailed", or "very_detailed" (default: "broad")
#'
#' @return A data frame with original annotations and suggested cell groupings
#'
#' @examples
#' \dontrun{
#' # Basic usage - merge annotations with broad groupings
#' result <- merge_annotations("annotation_results.csv", "merged_results.csv")
#'
#' # Use detailed groupings
#' result <- merge_annotations(
#'   csv_path = "annotation_results.csv",
#'   output_path = "merged_detailed.csv",
#'   detail_level = "detailed",
#'   provider = "anthropic"
#' )
#' }
#'
#' @export
merge_annotations <- function(csv_path,
                             output_path = NULL,
                             provider = "openai",
                             model = NULL,
                             api_key = NULL,
                             additional_context = NULL,
                             batch_size = 20L,
                             detail_level = "broad") {
  tryCatch({
    result <- py_cassia$merge_annotations(
      csv_path = csv_path,
      output_path = output_path,
      provider = provider,
      model = model,
      api_key = api_key,
      additional_context = additional_context,
      batch_size = as.integer(batch_size),
      detail_level = detail_level
    )
    # Convert Python DataFrame to R data frame
    return(reticulate::py_to_r(result))
  }, error = function(e) {
    error_msg <- paste("Error in merge_annotations:", e$message, "\n",
                      "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}


#' Merge All Annotation Levels
#'
#' Process all three detail levels (broad, detailed, very_detailed) and return a combined DataFrame
#' with all three grouping columns.
#'
#' @param csv_path Path to the CSV file containing cluster annotations
#' @param output_path Path to save the results (if NULL, returns DataFrame without saving)
#' @param provider LLM provider to use ("openai", "anthropic", or "openrouter") (default: "openai")
#' @param model Specific model to use (if NULL, uses default for provider)
#' @param api_key API key for the provider (if NULL, gets from environment)
#' @param additional_context Optional domain-specific context to help with annotation
#' @param batch_size Number of clusters to process in each LLM call (default: 20)
#'
#' @return A data frame with original annotations and all three grouping levels
#'
#' @examples
#' \dontrun{
#' # Process all grouping levels at once
#' result <- merge_annotations_all(
#'   csv_path = "annotation_results.csv",
#'   output_path = "merged_all_levels.csv",
#'   provider = "openrouter"
#' )
#' }
#'
#' @export
merge_annotations_all <- function(csv_path,
                                  output_path = NULL,
                                  provider = "openai",
                                  model = NULL,
                                  api_key = NULL,
                                  additional_context = NULL,
                                  batch_size = 20L) {
  tryCatch({
    result <- py_cassia$merge_annotations_all(
      csv_path = csv_path,
      output_path = output_path,
      provider = provider,
      model = model,
      api_key = api_key,
      additional_context = additional_context,
      batch_size = as.integer(batch_size)
    )
    # Convert Python DataFrame to R data frame
    return(reticulate::py_to_r(result))
  }, error = function(e) {
    error_msg <- paste("Error in merge_annotations_all:", e$message, "\n",
                      "Python traceback:", reticulate::py_last_error())
    stop(error_msg)
  })
}


# =============================================================================
# Logging Configuration Functions
# =============================================================================

#' Set CASSIA Log Level
#'
#' Control the verbosity of CASSIA error and warning messages.
#' By default, CASSIA shows informational messages, warnings, and errors.
#' Use this function to adjust the level of detail in log output.
#'
#' @param level Character string specifying the log level. One of:
#'   \itemize{
#'     \item "DEBUG" - Show all messages including detailed debug info
#'     \item "INFO" - Show informational messages, warnings, and errors (default)
#'     \item "WARNING" - Show only warnings and errors
#'     \item "ERROR" - Show only errors
#'     \item "QUIET" - Suppress almost all messages
#'   }
#'
#' @return Invisible NULL. Called for side effects.
#'
#' @examples
#' \dontrun{
#' # Enable verbose debugging
#' set_log_level("DEBUG")
#'
#' # Show only errors (quiet mode)
#' set_log_level("ERROR")
#'
#' # Suppress all messages
#' set_log_level("QUIET")
#'
#' # Reset to default
#' set_log_level("INFO")
#' }
#'
#' @export
set_log_level <- function(level = "INFO") {
  if (is.null(py_cassia)) {
    stop("CASSIA Python environment not initialized. Please restart R or run library(CASSIA) again.")
  }

  valid_levels <- c("DEBUG", "INFO", "WARNING", "ERROR", "QUIET")
  level <- toupper(level)

  if (!level %in% valid_levels) {
    stop(paste("Invalid log level. Must be one of:", paste(valid_levels, collapse = ", ")))
  }

  py_cassia$set_log_level(level)
  invisible(NULL)
}

#' Get CASSIA Logger
#'
#' Get or create a CASSIA logger instance for custom logging.
#' This is primarily for advanced users who want to add custom log messages
#' that integrate with CASSIA's logging system.
#'
#' @param name Character string specifying the logger name.
#'   Will be prefixed with "CASSIA." automatically.
#'   Default is "CASSIA" for the root logger.
#'
#' @return A Python logging.Logger object
#'
#' @examples
#' \dontrun{
#' # Get the default CASSIA logger
#' logger <- get_logger()
#'
#' # Get a named logger for custom module
#' my_logger <- get_logger("my_analysis")
#' }
#'
#' @export
get_logger <- function(name = "CASSIA") {
  if (is.null(py_cassia)) {
    stop("CASSIA Python environment not initialized. Please restart R or run library(CASSIA) again.")
  }

  py_cassia$get_logger(name)
}

#' Show a Warning to the User
#'
#' Display a warning message that is always visible, regardless of log level.
#' Uses Python's warnings module to ensure visibility.
#'
#' @param message Character string containing the warning message to display.
#'
#' @return Invisible NULL. Called for side effects.
#'
#' @examples
#' \dontrun{
#' warn_user("This operation may take a long time with large datasets.")
#' }
#'
#' @export
warn_user <- function(message) {
  if (is.null(py_cassia)) {
    warning(paste("[CASSIA]", message))
    return(invisible(NULL))
  }

  py_cassia$warn_user(message)
  invisible(NULL)
}