% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotator.R
\name{runCASSIA_score_batch}
\alias{runCASSIA_score_batch}
\title{Run Scoring with Progress Updates}
\usage{
runCASSIA_score_batch(
  input_file,
  output_file = NULL,
  max_workers = 4,
  model = "deepseek/deepseek-chat-v3-0324",
  provider = "openrouter",
  max_retries = 1
)
}
\arguments{
\item{input_file}{Path to input CSV file}

\item{output_file}{Path to output CSV file (optional)}

\item{max_workers}{Maximum number of parallel workers}

\item{model}{Model to use}

\item{provider}{AI provider to use ('openai', 'anthropic', or 'openrouter')}

\item{max_retries}{Maximum number of retries for failed analyses (default: 1)}
}
\value{
None
}
\description{
Run Scoring with Progress Updates
}
