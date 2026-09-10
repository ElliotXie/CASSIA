#!/usr/bin/env Rscript

# Build a balanced, label-blinded Seurat cohort directly from CELLxGENE Census.
# Truth is written separately and must only be supplied to score_integrated.py
# after the agent run has completed.

parse_args <- function(args) {
  defaults <- list(
    census_version = "2025-11-08",
    organism = "Homo sapiens",
    max_per_label = 150L,
    seed = 1729L
  )
  index <- 1L
  while (index <= length(args)) {
    key <- args[[index]]
    if (key == "--help") {
      defaults$help <- TRUE
      index <- index + 1L
      next
    }
    if (!startsWith(key, "--") || index == length(args)) {
      stop("Arguments must be --name value pairs", call. = FALSE)
    }
    name <- gsub("-", "_", substring(key, 3L))
    defaults[[name]] <- args[[index + 1L]]
    index <- index + 2L
  }
  defaults$max_per_label <- as.integer(defaults$max_per_label)
  defaults$seed <- as.integer(defaults$seed)
  defaults
}

usage <- function() {
  cat(paste(
    "Usage:",
    "  Rscript prepare_census_cohort.R --dataset-id UUID --labels labels.csv --out-dir DIR [options]",
    "",
    "Required:",
    "  --dataset-id UUID       CELLxGENE dataset_id",
    "  --labels CSV            source_fine_label,truth_broad_label allowlist",
    "  --out-dir DIR           Writes blind/input.rds and evaluator/truth.csv",
    "",
    "Options:",
    "  --census-version TAG    Default: 2025-11-08",
    "  --organism NAME         Default: Homo sapiens",
    "  --max-per-label N       Default: 150",
    "  --seed N                Default: 1729",
    sep = "\n"
  ))
}

require_value <- function(options, name) {
  value <- options[[name]]
  if (is.null(value) || !nzchar(value)) {
    stop(sprintf("Missing required --%s", gsub("_", "-", name)), call. = FALSE)
  }
  value
}

options <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(options$help)) {
  usage()
  quit(status = 0L)
}

dataset_id <- require_value(options, "dataset_id")
labels_path <- normalizePath(require_value(options, "labels"), mustWork = TRUE)
out_dir <- normalizePath(require_value(options, "out_dir"), mustWork = FALSE)
if (is.na(options$max_per_label) || options$max_per_label < 1L) {
  stop("--max-per-label must be at least 1", call. = FALSE)
}
if (is.na(options$seed)) stop("--seed must be an integer", call. = FALSE)

suppressPackageStartupMessages(library(cellxgene.census))
suppressPackageStartupMessages(library(Seurat))
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite is required", call. = FALSE)
}

labels <- read.csv(labels_path, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("source_fine_label", "truth_broad_label")
missing_columns <- setdiff(required_columns, colnames(labels))
if (length(missing_columns)) {
  stop(sprintf("Label allowlist is missing: %s", paste(missing_columns, collapse = ", ")), call. = FALSE)
}
if (anyDuplicated(labels$source_fine_label) || any(!nzchar(labels$source_fine_label))) {
  stop("source_fine_label values must be unique and non-empty", call. = FALSE)
}

dir.create(file.path(out_dir, "blind"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "evaluator"), recursive = TRUE, showWarnings = FALSE)

census <- open_soma(census_version = options$census_version)
on.exit(census$close(), add = TRUE)
organism_key <- if (options$organism == "Homo sapiens") {
  "homo_sapiens"
} else if (options$organism == "Mus musculus") {
  "mus_musculus"
} else {
  stop("--organism must be 'Homo sapiens' or 'Mus musculus'", call. = FALSE)
}

filter <- sprintf("dataset_id == '%s' && is_primary_data == TRUE", dataset_id)
obs <- census$get("census_data")$get(organism_key)$get("obs")$read(
  value_filter = filter,
  column_names = c(
    "soma_joinid", "observation_joinid", "cell_type",
    "cell_type_ontology_term_id", "donor_id", "assay", "disease"
  )
)$concat()
obs <- as.data.frame(obs)
obs <- obs[obs$cell_type %in% labels$source_fine_label, , drop = FALSE]
if (!nrow(obs)) stop("No cells matched the label allowlist", call. = FALSE)

missing_labels <- setdiff(labels$source_fine_label, unique(obs$cell_type))
if (length(missing_labels)) {
  stop(sprintf("Allowlisted labels absent from Census slice: %s", paste(missing_labels, collapse = "; ")), call. = FALSE)
}

set.seed(options$seed)
selected_rows <- unlist(lapply(split(seq_len(nrow(obs)), obs$cell_type), function(indices) {
  sample(indices, min(length(indices), options$max_per_label), replace = FALSE)
}), use.names = FALSE)
selected <- obs[selected_rows, , drop = FALSE]
selected <- selected[order(selected$cell_type, selected$soma_joinid), , drop = FALSE]

seurat_obj <- get_seurat(
  census,
  organism = options$organism,
  X_layers = c(counts = "raw"),
  obs_coords = selected$soma_joinid,
  obs_column_names = c(
    "observation_joinid", "cell_type", "cell_type_ontology_term_id",
    "donor_id", "assay", "disease"
  ),
  obsm_layers = FALSE,
  var_index = "feature_name"
)

metadata <- seurat_obj[[]]
truth <- data.frame(
  cell_id = rownames(metadata),
  truth_fine_label = as.character(metadata$cell_type),
  truth_broad_label = labels$truth_broad_label[match(metadata$cell_type, labels$source_fine_label)],
  truth_ontology_term_id = as.character(metadata$cell_type_ontology_term_id),
  stringsAsFactors = FALSE
)
if (anyNA(truth$truth_broad_label)) stop("Broad-label mapping failed", call. = FALSE)

truth_path <- file.path(out_dir, "evaluator", "truth.csv")
write.csv(truth, truth_path, row.names = FALSE, quote = TRUE)
Sys.chmod(truth_path, mode = "0600")

allowed_metadata <- intersect(c("orig.ident", "nCount_RNA", "nFeature_RNA", "donor_id", "assay", "disease"), colnames(metadata))
seurat_obj@meta.data <- metadata[, allowed_metadata, drop = FALSE]
if (any(c("cell_type", "cell_type_ontology_term_id") %in% colnames(seurat_obj[[]]))) {
  stop("Truth columns remained in blind object", call. = FALSE)
}

blind_path <- file.path(out_dir, "blind", "input.rds")
saveRDS(seurat_obj, blind_path, compress = "gzip")

counts <- as.data.frame(table(truth$truth_fine_label), stringsAsFactors = FALSE)
colnames(counts) <- c("truth_fine_label", "n_cells")
manifest <- list(
  schema_version = "cassia.census-heldout-cohort.v1",
  dataset_id = dataset_id,
  census_version = options$census_version,
  organism = options$organism,
  seed = options$seed,
  max_per_label = options$max_per_label,
  n_cells = nrow(truth),
  n_labels = nrow(counts),
  label_counts = counts,
  blind_rds = normalizePath(blind_path),
  blind_rds_md5 = unname(tools::md5sum(blind_path)),
  truth_csv = normalizePath(truth_path),
  truth_csv_md5 = unname(tools::md5sum(truth_path)),
  removed_truth_columns = c("cell_type", "cell_type_ontology_term_id"),
  allowed_metadata_columns = allowed_metadata
)
manifest_path <- file.path(out_dir, "cohort_manifest.json")
jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE)

cat(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE), "\n")
