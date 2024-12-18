% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotator.R
\name{runCASSIA_subclusters}
\alias{runCASSIA_subclusters}
\title{Process Subclusters}
\usage{
runCASSIA_subclusters(
  marker,
  major_cluster_info,
  output_name,
  model = "claude-3-5-sonnet-20241022",
  temperature = 0,
  provider = "anthropic",
  n_genes = 50L
)
}
\arguments{
\item{major_cluster_info}{Description of the major cluster type}

\item{model}{Model name for Claude API (default: "claude-3-5-sonnet-20241022")}

\item{temperature}{Temperature parameter for API calls (default: 0)}

\item{provider}{AI provider to use (default: "anthropic")}

\item{csv_file_path}{Path to input CSV file containing marker data}

\item{output_path}{Path where results CSV should be saved}
}
\value{
None. This function processes subclusters and saves results to a CSV file.
}
\description{
Process Subclusters
}
