#!/usr/bin/env Rscript
# Optional CASSIA accelerator: keep the shared R runtime behind localhost TCP.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("daemon.R requires <workdir> <port>")
}

command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)
if (length(file_arg) != 1L) stop("could not locate daemon.R")
script_path <- sub("^--file=", "", file_arg[[1]])
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))

workdir <- args[[1]]
run_mode <- "daemon"
port <- as.integer(args[[2]])
request_path <- NULL
source(file.path(script_dir, "runtime.R"), local = FALSE)
