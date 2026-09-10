#!/usr/bin/env Rscript
# Default CASSIA agent entrypoint: execute exactly one versioned transaction.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("transaction_worker.R requires <workdir> <request.json>")
}

command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)
if (length(file_arg) != 1L) stop("could not locate transaction_worker.R")
script_path <- sub("^--file=", "", file_arg[[1]])
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))

workdir <- args[[1]]
run_mode <- "once"
port <- NA_integer_
request_path <- args[[2]]
source(file.path(script_dir, "runtime.R"), local = FALSE)
