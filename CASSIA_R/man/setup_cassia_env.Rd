% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotator.R
\name{setup_cassia_env}
\alias{setup_cassia_env}
\title{Set up CASSIA Python Environment}
\usage{
setup_cassia_env(
  conda_env = NULL,
  python_version = "3.10",
  pip_packages = c("openai", "pandas", "numpy", "scikit-learn", "requests", "anthropic",
    "charset-normalizer"),
  method = "auto"
)
}
\arguments{
\item{conda_env}{The name of the environment to use. If NULL, uses the default from package configuration.}

\item{python_version}{The Python version to use. Default is "3.10".}

\item{pip_packages}{A character vector of pip packages to install.}

\item{method}{The method to use for environment setup. Options: "auto" (try virtualenv first, then conda), "virtualenv", "conda". Default is "auto".}
}
\value{
Invisible NULL. Called for side effects.
}
\description{
This function sets up the required Python environment for CASSIA.
It can be used to create a new environment or update an existing one.
By default, it tries virtualenv first (simpler, more reliable), then falls back to conda.
}
