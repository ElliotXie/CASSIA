% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotator.R
\name{symphonyCompare}
\alias{symphonyCompare}
\title{Symphony Compare - Advanced Multi-Model Cell Type Comparison with Consensus Building}
\usage{
symphonyCompare(
  tissue,
  celltypes,
  marker_set,
  species = "human",
  model_preset = "symphony",
  custom_models = NULL,
  output_dir = NULL,
  output_basename = NULL,
  enable_discussion = TRUE,
  max_discussion_rounds = 2L,
  consensus_threshold = 0.8,
  generate_report = TRUE,
  api_key = NULL,
  verbose = TRUE
)
}
\arguments{
\item{tissue}{Character string specifying the tissue type (e.g., "blood", "brain", "liver")}

\item{celltypes}{Character vector of 2-4 cell types to compare}

\item{marker_set}{Character vector or string of gene markers to analyze}

\item{species}{Character string specifying the species (default: "human")}

\item{model_preset}{Character string specifying preset model configuration. Options:
\itemize{
\item "symphony": High-performance ensemble (Claude, GPT-4, Gemini Pro)
\item "quartet": Balanced 4-model ensemble
\item "budget": Cost-effective models
\item "custom": Use custom_models list
}}

\item{custom_models}{Character vector of custom models to use (when model_preset="custom")}

\item{output_dir}{Character string specifying directory to save results (default: current directory)}

\item{output_basename}{Character string for base name of output files (auto-generated if NULL)}

\item{enable_discussion}{Logical indicating whether to enable automatic discussion rounds when no consensus (default: TRUE)}

\item{max_discussion_rounds}{Integer specifying maximum discussion rounds to perform (default: 2)}

\item{consensus_threshold}{Numeric value (0-1) specifying fraction of models that must agree for consensus (default: 0.8)}

\item{generate_report}{Logical indicating whether to generate interactive HTML report (default: TRUE)}

\item{api_key}{Character string for OpenRouter API key (uses environment variable if NULL)}

\item{verbose}{Logical indicating whether to print progress messages (default: TRUE)}
}
\value{
A list containing:
\itemize{
\item results: List of all model responses and scores
\item consensus: The consensus cell type (if reached)
\item confidence: Confidence level of the consensus (0-1)
\item csv_file: Path to the generated CSV file
\item html_file: Path to the generated HTML report (if enabled)
\item summary: Summary statistics of the comparison
\item dataframe: R data frame with structured results
}
}
\description{
Orchestrate multiple AI models to compare cell types with automatic consensus building.
This function conducts a comprehensive cell type comparison using multiple AI models in parallel,
automatically triggering discussion rounds when models disagree on the best matching cell type.
Think of it as a virtual panel of expert biologists debating and reaching consensus.
}
\examples{
\dontrun{
# Basic usage - let Symphony Compare handle everything
results <- symphonyCompare(
  tissue = "peripheral blood",
  celltypes = c("T cell", "B cell", "NK cell", "Monocyte"),
  marker_set = c("CD3", "CD4", "CD8", "CD19", "CD20", "CD16", "CD56", "CD14"),
  species = "human"
)

# Access the results
cat("Consensus:", results$consensus, "\n")
cat("Confidence:", sprintf("\%.1f\%\%", results$confidence * 100), "\n")

# Advanced usage with custom settings
results <- symphonyCompare(
  tissue = "brain",
  celltypes = c("Neuron", "Astrocyte", "Microglia", "Oligodendrocyte"),
  marker_set = c("RBFOX3", "GFAP", "IBA1", "OLIG2", "MAP2", "S100B", "CD11B", "MBP"),
  species = "mouse",
  model_preset = "quartet",  # Use 4 models instead of 3
  enable_discussion = TRUE,  # Enable automatic discussion rounds
  max_discussion_rounds = 3,  # Allow up to 3 discussion rounds
  consensus_threshold = 0.75,  # 75\% of models must agree
  output_dir = "./symphony_results",
  verbose = TRUE
)
}

}
