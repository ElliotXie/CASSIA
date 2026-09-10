# Shared CASSIA agent R runtime, sourced by transport-specific entrypoints.
#
# The entrypoint must define workdir, run_mode, port, and request_path before
# sourcing this file. transaction_worker.R selects one-shot mode; daemon.R
# selects the optional TCP accelerator.

# Optional daemon lifecycle:
#   1. Bind a localhost port (OS-picked), write port + pid into
#      <workdir>/session.json so clients can find us.
#   2. Loop: accept one TCP connection, read one JSON line, dispatch,
#      write one JSON line back, close connection.
#   3. Exit cleanly on shutdown op.
#
# Protocol per request:  {"id": str, "op": str, "args": {...}}
# Protocol per response: {"id": str, "ok": bool, "elapsed_s": num,
#                         "result": ... | "error": str}

suppressPackageStartupMessages({
  library(jsonlite)
})

required_runtime_values <- c(
  "workdir", "run_mode", "port", "request_path", "script_dir"
)
missing_runtime_values <- required_runtime_values[
  !vapply(
    required_runtime_values,
    function(name) exists(name, envir = .GlobalEnv, inherits = FALSE),
    logical(1)
  )
]
if (length(missing_runtime_values) > 0L) {
  stop(sprintf(
    "runtime.R must be sourced by a CASSIA entrypoint; missing: %s",
    paste(missing_runtime_values, collapse = ", ")
  ))
}
if (!(run_mode %in% c("once", "daemon"))) {
  stop(sprintf("unsupported CASSIA R runtime mode: %s", run_mode))
}
workdir <- normalizePath(workdir, mustWork = TRUE)
if (identical(run_mode, "daemon") && (is.na(port) || port <= 0L)) {
  stop(sprintf("bad port: %s", port))
}
if (identical(run_mode, "once") &&
    (is.null(request_path) || !file.exists(request_path))) {
  stop("one-shot mode requires an existing request JSON file")
}
source(file.path(script_dir, "policy.R"), local = FALSE)

# In daemon mode stderr is the log and stdout is unused. In one-shot mode the
# final stdout JSON object is the transaction response.
sink(stderr(), type = "message")
log_path     <- file.path(workdir, "daemon.log")
audit_path   <- file.path(workdir, "audit.jsonl")
ready_path   <- file.path(workdir, "daemon.ready")

log_msg <- function(...) {
  line <- paste0("[", format(Sys.time(), "%H:%M:%OS3"), "] ",
                 paste(..., sep = " "), "\n")
  cat(line, file = stderr())
  try(cat(line, file = log_path, append = TRUE), silent = TRUE)
}

audit_append <- function(entry) {
  try(cat(toJSON(entry, auto_unbox = TRUE, null = "null", na = "null"),
          "\n", file = audit_path, append = TRUE, sep = ""),
      silent = TRUE)
}

state <- new.env(parent = emptyenv())
state$seurat <- NULL
state$rds_path <- NULL
state$phase <- "INIT"
state$step <- 0L
# clusters: list of lists with $uuid, $seurat_id, $n_cells, $alive,
#           $origin, $parent_uuids, $born_at_step, $died_at_step, $died_reason,
#           $label (NULL or list with $name $confidence $markers $reason
#                   $set_at_step $skip), $skip (TRUE if explicitly skipped).
state$clusters <- list()
state$libs_loaded <- FALSE
state$autozyme_active <- FALSE
# Top-N markers cache, keyed by cluster uuid → character vector of gene names.
# Populated on every `markers` op so qa() can verify cited markers cheaply.
state$marker_cache <- new.env(parent = emptyenv())
state$marker_cache_top_n <- 50L
# Counter-based UUID generator (Seurat's FindSubCluster calls set.seed(0)
# internally, which collides sample.int-based uuids inside the same op).
# Do modulo on the numeric to avoid 32-bit integer overflow → NA → "  NA"
# in the hex prefix.
state$uuid_counter <- 0L
state$session_prefix <- sprintf("%04x",
  as.integer((as.numeric(Sys.time()) * 1e3) %% 65536))
state$checkpoint_schema_version <- "cassia.agent.checkpoint.v1"
state$policy <- cassia_default_policy()

# ── helpers ──────────────────────────────────────────────────────────────

uuid8 <- function() {
  # 8-char hex slug: session_prefix (4) + monotonic counter (4). Guaranteed
  # unique within a session and independent of Seurat's internal RNG resets.
  state$uuid_counter <- state$uuid_counter + 1L
  paste0(state$session_prefix,
         sprintf("%04x", bitwAnd(state$uuid_counter, 0xFFFFL)))
}

ensure_libs <- function(use_autozyme = FALSE) {
  if (!state$libs_loaded) {
    suppressPackageStartupMessages(library(Seurat))
    state$libs_loaded <- TRUE
  }
  if (use_autozyme && !state$autozyme_active) {
    suppressPackageStartupMessages(library(autozyme))
    tryCatch(autozyme::inject_all(), error = function(e) {
      log_msg("autozyme inject_all failed:", conditionMessage(e))
    })
    state$autozyme_active <- TRUE
  }
}

require_seurat <- function() {
  if (is.null(state$seurat)) stop("no Seurat object loaded; call init first")
}

marker_cache_as_list <- function() {
  keys <- ls(envir = state$marker_cache, all.names = TRUE)
  out <- setNames(vector("list", length(keys)), keys)
  for (key in keys) out[[key]] <- state$marker_cache[[key]]
  out
}

restore_marker_cache <- function(values) {
  cache <- new.env(parent = emptyenv())
  if (!is.null(values) && length(values) > 0L) {
    for (key in names(values)) cache[[key]] <- values[[key]]
  }
  cache
}

checkpoint_payload <- function() {
  list(
    schema_version = state$checkpoint_schema_version,
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    seurat = state$seurat,
    rds_path = state$rds_path,
    phase = state$phase,
    step = state$step,
    clusters = state$clusters,
    marker_cache = marker_cache_as_list(),
    marker_cache_top_n = state$marker_cache_top_n,
    uuid_counter = state$uuid_counter,
    session_prefix = state$session_prefix,
    autozyme_active = state$autozyme_active,
    policy = state$policy
  )
}

atomic_save_rds <- function(value, path) {
  path <- normalizePath(path, mustWork = FALSE)
  parent <- dirname(path)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = parent)
  on.exit(try(unlink(tmp), silent = TRUE), add = TRUE)
  saveRDS(value, tmp)
  if (file.exists(path) && unlink(path) != 0L) {
    stop(sprintf("could not replace checkpoint: %s", path))
  }
  if (!file.rename(tmp, path)) {
    stop(sprintf("could not atomically install checkpoint: %s", path))
  }
  path
}

atomic_write_csv <- function(value, path) {
  path <- normalizePath(path, mustWork = FALSE)
  parent <- dirname(path)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = parent)
  on.exit(try(unlink(tmp), silent = TRUE), add = TRUE)
  write.csv(value, tmp, row.names = FALSE, quote = TRUE)
  if (file.exists(path) && unlink(path) != 0L) {
    stop(sprintf("could not replace CSV: %s", path))
  }
  if (!file.rename(tmp, path)) stop(sprintf("could not install CSV: %s", path))
  path
}

atomic_write_json <- function(value, path) {
  path <- normalizePath(path, mustWork = FALSE)
  parent <- dirname(path)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE)
  tmp <- tempfile(pattern = paste0(basename(path), "."), tmpdir = parent)
  on.exit(try(unlink(tmp), silent = TRUE), add = TRUE)
  writeLines(toJSON(value, auto_unbox = TRUE, null = "null", na = "null",
                    pretty = TRUE, force = TRUE), tmp)
  if (file.exists(path) && unlink(path) != 0L) {
    stop(sprintf("could not replace JSON: %s", path))
  }
  if (!file.rename(tmp, path)) stop(sprintf("could not install JSON: %s", path))
  path
}

cluster_by_uuid <- function(uuid) {
  for (c in state$clusters) {
    if (identical(c$uuid, uuid)) return(c)
  }
  NULL
}

alive_clusters <- function() {
  Filter(function(c) isTRUE(c$alive), state$clusters)
}

cluster_uuid_for_seurat_id <- function(sid) {
  for (c in state$clusters) {
    if (isTRUE(c$alive) && identical(as.character(c$seurat_id), as.character(sid))) {
      return(c$uuid)
    }
  }
  NA_character_
}

register_initial_clusters <- function() {
  obj <- state$seurat
  ids <- levels(Idents(obj))
  sizes <- as.list(table(Idents(obj)))
  out <- list()
  state$step <- state$step + 1L
  for (sid in ids) {
    out[[length(out) + 1L]] <- list(
      uuid = uuid8(),
      seurat_id = sid,
      n_cells = as.integer(sizes[[sid]]),
      alive = TRUE,
      origin = "init",
      parent_uuids = list(),
      born_at_step = state$step,
      died_at_step = NA_integer_,
      died_reason = NA_character_
    )
  }
  state$clusters <- out
  invisible(out)
}

resolve_idents <- function(spec) {
  # spec: a single uuid or seurat_id (string). Returns the seurat_id usable by
  # FindMarkers, or NULL if "rest"/NA. Errors if uuid unknown / cluster dead.
  if (is.null(spec) || is.na(spec) || identical(spec, "rest") ||
      identical(spec, "")) return(NULL)
  if (nchar(spec) == 8L && grepl("^[0-9a-f]{8}$", spec)) {
    c <- cluster_by_uuid(spec)
    if (is.null(c)) stop(sprintf("unknown cluster uuid: %s", spec))
    if (!isTRUE(c$alive)) stop(sprintf("cluster %s is dead", spec))
    return(as.character(c$seurat_id))
  }
  # Otherwise treat as raw seurat ident; verify it's alive.
  if (is.na(cluster_uuid_for_seurat_id(spec))) {
    stop(sprintf("seurat id '%s' is not an alive cluster", spec))
  }
  spec
}

resolve_ident2 <- function(spec) {
  # Pairwise vs: a single uuid, comma-separated uuids, or "rest"/NULL.
  if (is.null(spec) || is.na(spec) || identical(spec, "rest") ||
      identical(spec, "")) return(NULL)
  parts <- strsplit(spec, "[,\\s]+", perl = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  out <- vapply(parts, resolve_idents, character(1))
  if (length(out) == 1L) return(out)
  out
}

# ── ops ──────────────────────────────────────────────────────────────────

op_ping <- function(args) {
  list(pong = TRUE,
       phase = state$phase,
       step = state$step,
       libs_loaded = state$libs_loaded,
       autozyme_active = state$autozyme_active,
       seurat_loaded = !is.null(state$seurat))
}

op_init <- function(args) {
  path <- args$rds_path
  if (is.null(path) || !file.exists(path)) {
    stop(sprintf("rds path not found: %s", path))
  }
  use_autozyme <- isTRUE(args$autozyme)
  ensure_libs(use_autozyme = use_autozyme)
  obj <- readRDS(path)
  if (!inherits(obj, "Seurat")) {
    stop(sprintf("loaded object is not a Seurat object: class=%s",
                 paste(class(obj), collapse = ",")))
  }
  # ``init`` may be called again on a long-lived daemon.  Treat it as a fresh
  # scientific session so cluster UUIDs, caches, labels, and automation policy
  # can never leak from the previously loaded object.
  state$phase <- "INIT"
  state$step <- 0L
  state$clusters <- list()
  state$marker_cache <- new.env(parent = emptyenv())
  state$uuid_counter <- 0L
  state$session_prefix <- sprintf("%04x",
    as.integer((as.numeric(Sys.time()) * 1e3) %% 65536))
  state$policy <- cassia_default_policy()
  state$seurat <- obj
  state$rds_path <- path
  has_idents <- length(levels(Idents(obj))) > 1L ||
                !identical(as.character(levels(Idents(obj))), basename(path))
  if (length(levels(Idents(obj))) >= 1L && has_idents) {
    register_initial_clusters()
    state$phase <- "ORIENT"
  } else {
    # Raw object — preprocess deferred to a separate op (not iteration 1).
    state$phase <- "INIT"
  }
  list(
    rds_path = path,
    seurat_class = class(obj),
    n_cells = ncol(obj),
    n_features = nrow(obj),
    n_clusters_alive = length(alive_clusters()),
    phase = state$phase,
    autozyme_active = state$autozyme_active
  )
}

op_status <- function(args) {
  out <- list(
    phase = state$phase,
    step = state$step,
    rds_path = state$rds_path %||% NA_character_,
    n_cells = if (is.null(state$seurat)) 0L else ncol(state$seurat),
    n_clusters_alive = length(alive_clusters()),
    n_clusters_total = length(state$clusters),
    autozyme_active = state$autozyme_active
  )
  if (!is.null(state$seurat)) {
    out$seurat_mb <- round(as.numeric(object.size(state$seurat)) / 1024 / 1024, 2)
  }
  out
}

op_checkpoint <- function(args) {
  require_seurat()
  requested <- args$path
  if (is.null(requested) || !nzchar(requested)) {
    requested <- file.path(
      workdir, "checkpoints", sprintf("step_%06d.rds", state$step)
    )
  }
  path <- atomic_save_rds(checkpoint_payload(), requested)
  list(
    schema_version = state$checkpoint_schema_version,
    path = path,
    md5 = unname(tools::md5sum(path)),
    step = state$step,
    phase = state$phase,
    n_cells = ncol(state$seurat),
    n_clusters_alive = length(alive_clusters())
  )
}

op_restore <- function(args) {
  requested <- args$path
  if (is.null(requested) || !nzchar(requested)) stop("checkpoint path required")
  path <- normalizePath(requested, mustWork = TRUE)
  payload <- readRDS(path)
  if (!is.list(payload) ||
      !identical(payload$schema_version, state$checkpoint_schema_version)) {
    stop(sprintf("unsupported checkpoint schema in %s", path))
  }
  if (is.null(payload$seurat) || !inherits(payload$seurat, "Seurat")) {
    stop(sprintf("checkpoint does not contain a Seurat object: %s", path))
  }

  ensure_libs(use_autozyme = isTRUE(payload$autozyme_active))
  restored_clusters <- payload$clusters %||% list()
  alive_sids <- sort(vapply(
    Filter(function(c) isTRUE(c$alive), restored_clusters),
    function(c) as.character(c$seurat_id), character(1)
  ))
  object_sids <- sort(as.character(levels(Idents(payload$seurat))))
  if (!identical(alive_sids, object_sids)) {
    stop("checkpoint cluster registry does not match Seurat identities")
  }

  state$seurat <- payload$seurat
  state$rds_path <- payload$rds_path
  state$phase <- payload$phase %||% "ORIENT"
  state$step <- as.integer(payload$step %||% 0L)
  state$clusters <- restored_clusters
  state$marker_cache <- restore_marker_cache(payload$marker_cache)
  state$marker_cache_top_n <- as.integer(payload$marker_cache_top_n %||% 50L)
  state$uuid_counter <- as.integer(payload$uuid_counter %||% length(restored_clusters))
  state$session_prefix <- as.character(payload$session_prefix %||% state$session_prefix)
  state$policy <- payload$policy %||% cassia_default_policy()

  list(
    schema_version = state$checkpoint_schema_version,
    path = path,
    md5 = unname(tools::md5sum(path)),
    step = state$step,
    phase = state$phase,
    n_cells = ncol(state$seurat),
    n_clusters_alive = length(alive_clusters())
  )
}

op_policy <- function(args) {
  state$policy <- cassia_resolve_policy(args, state$policy)
  state$policy
}

op_export <- function(args) {
  require_seurat()
  requested <- args$out_dir
  if (is.null(requested) || !nzchar(requested)) {
    requested <- file.path(workdir, "artifacts")
  }
  out_dir <- normalizePath(requested, mustWork = FALSE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  partition_id <- args$partition_id
  if (is.null(partition_id) || !nzchar(partition_id)) {
    partition_id <- sprintf("p%06d", state$step)
  }

  object_sids <- as.character(Idents(state$seurat))
  object_cells <- colnames(state$seurat)
  uuid_by_sid <- setNames(
    vapply(alive_clusters(), function(c) c$uuid, character(1)),
    vapply(alive_clusters(), function(c) as.character(c$seurat_id), character(1))
  )
  cluster_ids <- unname(uuid_by_sid[object_sids])
  if (any(is.na(cluster_ids))) {
    stop("cannot export: at least one cell identity has no alive cluster UUID")
  }
  memberships <- data.frame(
    cell_id = object_cells,
    partition_id = rep(partition_id, length(object_cells)),
    cluster_id = cluster_ids,
    engine_cluster_id = object_sids,
    stringsAsFactors = FALSE
  )
  memberships <- memberships[order(memberships$cell_id), , drop = FALSE]
  membership_path <- atomic_write_csv(
    memberships, file.path(out_dir, paste0(partition_id, "_memberships.csv"))
  )
  membership_md5 <- unname(tools::md5sum(membership_path))

  cluster_rows <- lapply(alive_clusters(), function(c) list(
    partition_id = partition_id,
    cluster_id = c$uuid,
    engine_cluster_id = as.character(c$seurat_id),
    parent_cluster_ids = c$parent_uuids %||% list(),
    n_cells = c$n_cells,
    origin = c$origin,
    born_at_step = c$born_at_step,
    membership_file_md5 = membership_md5
  ))
  clusters_path <- atomic_write_json(
    list(
      schema_version = "cassia.partition.v1",
      partition_id = partition_id,
      step = state$step,
      memberships_path = membership_path,
      memberships_md5 = membership_md5,
      clusters = cluster_rows
    ),
    file.path(out_dir, paste0(partition_id, "_clusters.json"))
  )
  list(
    schema_version = "cassia.partition.v1",
    partition_id = partition_id,
    step = state$step,
    n_cells = nrow(memberships),
    n_clusters = length(cluster_rows),
    memberships_path = membership_path,
    memberships_md5 = membership_md5,
    clusters_path = clusters_path
  )
}

op_clusters <- function(args) {
  if (isTRUE(args$alive_only) || is.null(args$alive_only)) {
    cs <- alive_clusters()
  } else {
    cs <- state$clusters
  }
  # Strip empty parent_uuids list to avoid JSON ambiguity.
  lapply(cs, function(c) {
    c$parent_uuids <- if (length(c$parent_uuids) == 0) list() else c$parent_uuids
    c
  })
}

cache_markers_for_cluster <- function(uuid, marker_df) {
  # Store a slim per-cluster top-N data.frame keyed by gene so gene reverse
  # lookup can return logFC/pct stats without re-running FindMarkers.
  if (is.na(uuid) || is.null(uuid)) return(invisible(NULL))
  keep <- c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
  keep <- keep[keep %in% colnames(marker_df)]
  trimmed <- marker_df[, keep, drop = FALSE]
  trimmed <- head(trimmed, state$marker_cache_top_n)
  rownames(trimmed) <- NULL
  state$marker_cache[[uuid]] <- trimmed
  invisible(NULL)
}

cached_gene_names <- function(uuid) {
  df <- state$marker_cache[[uuid]]
  if (is.null(df)) return(character(0))
  as.character(df$gene)
}

op_markers <- function(args) {
  require_seurat()
  top_n <- if (is.null(args$top_n)) 20L else as.integer(args$top_n)
  cluster <- args$cluster
  vs <- args$vs

  if (identical(cluster, "all") || identical(cluster, "ALL")) {
    # For caching we always compute the full top-50 per cluster regardless of
    # what the agent asked for; the displayed result honors top_n.
    m <- FindAllMarkers(state$seurat, only.pos = TRUE, verbose = FALSE)
    if (nrow(m) == 0) return(list(scope = "all", n = 0L, markers = list()))
    m <- m[order(m$cluster, -m$avg_log2FC), , drop = FALSE]
    split_m <- split(m, m$cluster)
    for (sid in names(split_m)) {
      uuid <- cluster_uuid_for_seurat_id(sid)
      if (!is.na(uuid)) {
        cache_markers_for_cluster(uuid, split_m[[sid]])
      }
    }
    trimmed <- do.call(rbind, lapply(split_m, function(d) head(d, top_n)))
    rownames(trimmed) <- NULL
    trimmed$cluster_uuid <- vapply(as.character(trimmed$cluster),
                                   cluster_uuid_for_seurat_id, character(1))
    return(list(scope = "all", n = nrow(trimmed), markers = trimmed))
  }

  ident1 <- resolve_idents(cluster)
  ident2 <- resolve_ident2(vs)
  is_vs_rest <- is.null(ident2)
  m <- FindMarkers(
    state$seurat,
    ident.1 = ident1,
    ident.2 = ident2,
    only.pos = isTRUE(args$only_pos),
    verbose = FALSE
  )
  m$gene <- rownames(m)
  m <- m[order(-m$avg_log2FC), , drop = FALSE]
  ident1_uuid <- cluster_uuid_for_seurat_id(ident1)
  # Cache vs-rest results only — pairwise markers reflect a different question
  # and shouldn't be used to verify generic label claims.
  if (is_vs_rest && !is.na(ident1_uuid) && nrow(m) > 0) {
    cache_markers_for_cluster(ident1_uuid, m)
  }
  m <- head(m, top_n)
  rownames(m) <- NULL
  list(
    scope = "single",
    ident_1 = ident1,
    ident_2 = if (is.null(ident2)) "rest" else paste(ident2, collapse = ","),
    ident_1_uuid = ident1_uuid,
    n = nrow(m),
    markers = m
  )
}

.gene_lookup_one <- function(gene, alive) {
  hits <- list()
  for (c in alive) {
    df <- state$marker_cache[[c$uuid]]
    if (is.null(df)) next
    row_idx <- which(df$gene == gene)
    if (length(row_idx) == 0) next
    r <- df[row_idx[1], , drop = FALSE]
    hits[[length(hits) + 1L]] <- list(
      cluster_uuid = c$uuid,
      seurat_id = c$seurat_id,
      n_cells = c$n_cells,
      cluster_label = if (is.null(c$label)) NA_character_
                      else if (isTRUE(c$label$skip)) "[skipped]"
                      else c$label$name,
      rank = as.integer(row_idx[1]),
      cache_size = nrow(df),
      avg_log2FC = r$avg_log2FC,
      pct.1 = r$`pct.1`,
      pct.2 = r$`pct.2`,
      p_val_adj = r$p_val_adj
    )
  }
  if (length(hits) > 1) {
    ord <- order(-vapply(hits, function(h) as.numeric(h$avg_log2FC), numeric(1)))
    hits <- hits[ord]
  }
  hits
}

op_preprocess <- function(args) {
  require_seurat()
  if (isTRUE(state$policy$enabled) && isTRUE(state$policy$lock_preprocess)) {
    stop("automation policy locked global preprocessing for this partition")
  }
  if (length(state$clusters) > 0 &&
      any(vapply(state$clusters, function(c) !is.null(c$label), logical(1)))) {
    stop("preprocess would invalidate existing labels; unlabel everything or start a fresh workdir")
  }
  obj <- state$seurat
  n_var <- if (is.null(args$n_var)) 2000L else as.integer(args$n_var)
  n_pcs <- if (is.null(args$n_pcs)) 30L else as.integer(args$n_pcs)
  res <- if (is.null(args$resolution)) 0.5 else as.numeric(args$resolution)
  seed <- if (is.null(args$seed)) 17L else as.integer(args$seed)
  algorithm <- if (is.null(args$algorithm)) 1L else as.integer(args$algorithm)
  if (!(algorithm %in% 1:4)) stop("algorithm must be one of 1, 2, 3, or 4")
  do_scale <- !isTRUE(args$skip_scale)

  set.seed(seed)
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = n_var, verbose = FALSE)
  if (do_scale) {
    obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  }
  obj <- RunPCA(obj, npcs = n_pcs, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:n_pcs, verbose = FALSE)
  obj <- FindClusters(
    obj,
    resolution = res,
    algorithm = algorithm,
    random.seed = seed,
    verbose = FALSE
  )
  state$seurat <- obj
  # Rebuild cluster registry from scratch — old UUIDs do not survive.
  state$clusters <- list()
  state$marker_cache <- new.env(parent = emptyenv())
  register_initial_clusters()
  state$phase <- "ORIENT"
  list(
    n_var = n_var,
    n_pcs = n_pcs,
    resolution = res,
    seed = seed,
    algorithm = algorithm,
    n_clusters = length(alive_clusters()),
    cluster_sizes = setNames(
      lapply(alive_clusters(), function(c) c$n_cells),
      vapply(alive_clusters(), function(c) c$uuid, character(1))
    ),
    phase = state$phase
  )
}

op_gene <- function(args) {
  require_seurat()
  # Accept either a single `gene` or a `genes` list/vector. Strings can be
  # comma- or whitespace-separated for ergonomics.
  raw <- args$genes
  if (is.null(raw)) raw <- args$gene
  if (is.null(raw)) stop("gene(s) required (use `gene` or `genes`)")
  if (is.list(raw)) raw <- unlist(raw)
  raw <- as.character(raw)
  pieces <- unlist(strsplit(raw, "[,\\s]+", perl = TRUE))
  pieces <- trimws(pieces)
  genes <- unique(pieces[nzchar(pieces)])
  if (length(genes) == 0L) stop("no gene names parsed from input")

  alive <- alive_clusters()
  uncached <- sum(vapply(alive, function(c) {
    length(cached_gene_names(c$uuid)) == 0L
  }, integer(1)))

  per_gene <- list()
  for (g in genes) {
    hits <- .gene_lookup_one(g, alive)
    per_gene[[length(per_gene) + 1L]] <- list(
      gene = g,
      n_hits = length(hits),
      hits = hits
    )
  }

  warning_str <- NULL
  if (uncached > 0L) {
    warning_str <- sprintf(
      "%d of %d alive clusters have no cached markers; run `markers all` first",
      uncached, length(alive))
  }
  list(
    n_genes = length(genes),
    n_clusters_alive = length(alive),
    n_clusters_with_cache = length(alive) - uncached,
    results = per_gene,
    cache_warning = warning_str
  )
}

decisions_dir <- function() {
  d <- file.path(workdir, "decisions")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}

outputs_dir <- function() {
  d <- file.path(workdir, "outputs")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}

cluster_index_by_uuid <- function(uuid) {
  for (i in seq_along(state$clusters)) {
    if (identical(state$clusters[[i]]$uuid, uuid)) return(i)
  }
  NA_integer_
}

write_decision_file <- function(uuid) {
  idx <- cluster_index_by_uuid(uuid)
  if (is.na(idx)) return(invisible(NULL))
  c <- state$clusters[[idx]]
  cache_genes <- cached_gene_names(uuid)
  payload <- list(
    uuid = uuid,
    seurat_id = c$seurat_id,
    n_cells = c$n_cells,
    alive = isTRUE(c$alive),
    markers_top_cached = as.list(cache_genes),
    label = c$label,
    updated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  path <- file.path(decisions_dir(), paste0(uuid, ".json"))
  writeLines(toJSON(payload, auto_unbox = TRUE, null = "null", na = "null",
                    pretty = TRUE, force = TRUE), path)
  invisible(path)
}

op_merge <- function(args) {
  require_seurat()
  cassia_enforce_topology_policy(state$policy, "merge")
  raw <- args$uuids
  if (is.null(raw)) stop("uuids required")
  if (is.list(raw)) raw <- unlist(raw)
  uuids <- as.character(unlist(strsplit(as.character(raw), "[,\\s]+", perl = TRUE)))
  uuids <- unique(trimws(uuids))
  uuids <- uuids[nzchar(uuids)]
  if (length(uuids) < 2L) stop("merge requires at least 2 UUIDs")
  reason <- args$reason %||% ""
  if (!nzchar(reason)) stop("--reason required for merge (will be recorded in audit)")

  # Validate all alive
  indices <- integer(0)
  sids <- character(0)
  for (u in uuids) {
    idx <- cluster_index_by_uuid(u)
    if (is.na(idx)) stop(sprintf("unknown uuid: %s", u))
    c <- state$clusters[[idx]]
    if (!isTRUE(c$alive)) stop(sprintf("cluster %s is already dead", u))
    if (!is.null(c$label)) {
      stop(sprintf("cluster %s already has a label; unlabel before merging", u))
    }
    indices <- c(indices, idx)
    sids <- c(sids, as.character(c$seurat_id))
  }

  state$step <- state$step + 1L
  new_uuid <- uuid8()
  new_sid <- paste0("m_", new_uuid)
  obj <- state$seurat
  cells_to_move <- WhichCells(obj, idents = sids)
  Idents(obj, cells = cells_to_move) <- new_sid
  Idents(obj) <- droplevels(Idents(obj))
  state$seurat <- obj
  total_cells <- length(cells_to_move)

  # Kill old clusters, drop their cached markers (now stale).
  for (i in indices) {
    state$clusters[[i]]$alive <- FALSE
    state$clusters[[i]]$died_at_step <- state$step
    state$clusters[[i]]$died_reason <- sprintf("merged into %s", new_uuid)
    rm(list = state$clusters[[i]]$uuid, envir = state$marker_cache,
       inherits = FALSE)
  }

  # Register the new merged cluster.
  state$clusters[[length(state$clusters) + 1L]] <- list(
    uuid = new_uuid,
    seurat_id = new_sid,
    n_cells = total_cells,
    alive = TRUE,
    origin = "merge",
    parent_uuids = as.list(uuids),
    born_at_step = state$step,
    died_at_step = NA_integer_,
    died_reason = NA_character_
  )
  state$policy <- cassia_record_topology_edit(state$policy)

  audit_append(list(event = "merge", step = state$step,
                    uuids = as.list(uuids), new_uuid = new_uuid,
                    n_cells = total_cells, reason = reason,
                    ts = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3")))
  list(
    new_uuid = new_uuid,
    new_seurat_id = new_sid,
    n_cells = total_cells,
    merged_uuids = as.list(uuids),
    note = "run `markers <new_uuid>` to populate marker cache for the merged cluster"
  )
}

op_subcluster <- function(args) {
  require_seurat()
  cassia_enforce_topology_policy(state$policy, "subcluster")
  uuid <- args$uuid
  if (is.null(uuid)) stop("uuid required")
  idx <- cluster_index_by_uuid(uuid)
  if (is.na(idx)) stop(sprintf("unknown uuid: %s", uuid))
  c <- state$clusters[[idx]]
  if (!isTRUE(c$alive)) stop(sprintf("cluster %s is dead", uuid))
  if (!is.null(c$label)) {
    stop(sprintf("cluster %s already labeled; unlabel before subclustering", uuid))
  }
  res <- if (is.null(args$resolution)) 0.3 else as.numeric(args$resolution)
  algorithm <- if (is.null(args$algorithm)) 1L else as.integer(args$algorithm)
  if (!(algorithm %in% 1:4)) stop("algorithm must be one of 1, 2, 3, or 4")
  reason <- args$reason %||% ""
  if (!nzchar(reason)) stop("--reason required for subcluster")

  obj <- state$seurat
  graph_names <- Graphs(obj)
  snn <- graph_names[grepl("_snn$", graph_names)]
  if (length(snn) == 0L) {
    stop("no SNN graph found; run `preprocess` first (or load a preprocessed object)")
  }
  requested_graph <- args$graph
  snn_name <- if (is.null(requested_graph) || !nzchar(requested_graph)) snn[1] else requested_graph
  if (!(snn_name %in% graph_names)) {
    stop(sprintf("graph '%s' was not found; available graphs: %s",
                 snn_name, paste(graph_names, collapse = ", ")))
  }
  parent_sid <- as.character(c$seurat_id)
  tmp_col <- paste0("cassia_sub_", uuid8())
  obj <- FindSubCluster(obj, cluster = parent_sid, graph.name = snn_name,
                        subcluster.name = tmp_col, resolution = res,
                        algorithm = algorithm)
  new_idents <- obj[[tmp_col]][, 1]
  Idents(obj) <- new_idents
  Idents(obj) <- droplevels(Idents(obj))

  # FindSubCluster names children "<parent>_<n>". Collect all idents that
  # start with parent_sid + "_" and register each as a new cluster.
  child_ids <- grep(paste0("^", parent_sid, "_"), as.character(levels(Idents(obj))),
                    value = TRUE)
  if (length(child_ids) < 2L) {
    stop(sprintf(
      paste0("subcluster produced %d child at resolution=%s; transaction rolled ",
             "back and topology budget was not consumed. Try a higher resolution ",
             "only if marker evidence still supports a split"),
      length(child_ids), res
    ))
  }
  child_sizes <- vapply(
    child_ids, function(sid) length(WhichCells(obj, idents = sid)), integer(1)
  )
  cassia_validate_split_policy(child_sizes, state$policy)

  state$seurat <- obj
  state$step <- state$step + 1L
  state$clusters[[idx]]$alive <- FALSE
  state$clusters[[idx]]$died_at_step <- state$step
  state$clusters[[idx]]$died_reason <- sprintf("subclustered into children at res=%s", res)
  if (exists(uuid, envir = state$marker_cache, inherits = FALSE)) {
    rm(list = uuid, envir = state$marker_cache, inherits = FALSE)
  }

  children <- list()
  for (sid in child_ids) {
    n <- length(WhichCells(obj, idents = sid))
    new_uuid <- uuid8()
    rec <- list(
      uuid = new_uuid,
      seurat_id = sid,
      n_cells = n,
      alive = TRUE,
      origin = "subcluster",
      parent_uuids = list(uuid),
      born_at_step = state$step,
      died_at_step = NA_integer_,
      died_reason = NA_character_
    )
    state$clusters[[length(state$clusters) + 1L]] <- rec
    children[[length(children) + 1L]] <- list(uuid = new_uuid,
                                              seurat_id = sid,
                                              n_cells = n)
  }
  state$policy <- cassia_record_topology_edit(state$policy)

  audit_append(list(event = "subcluster", step = state$step,
                    parent_uuid = uuid, resolution = res,
                    algorithm = algorithm, graph = snn_name, reason = reason,
                    n_children = length(children),
                    children = children,
                    ts = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3")))
  list(
    parent_uuid = uuid,
    resolution = res,
    algorithm = algorithm,
    graph = snn_name,
    n_children = length(children),
    children = children,
    note = "run `markers <child_uuid>` (or `markers all`) to populate marker cache"
  )
}

op_label <- function(args) {
  require_seurat()
  uuid <- args$uuid
  if (is.null(uuid)) stop("uuid required")
  idx <- cluster_index_by_uuid(uuid)
  if (is.na(idx)) stop(sprintf("unknown cluster uuid: %s", uuid))
  if (!isTRUE(state$clusters[[idx]]$alive)) {
    stop(sprintf("cluster %s is dead, cannot label", uuid))
  }
  skip <- isTRUE(args$skip)
  reason <- args$reason %||% ""

  if (skip) {
    if (!nzchar(reason)) stop("--skip requires --reason")
    state$step <- state$step + 1L
    state$clusters[[idx]]$label <- list(
      name = NA_character_,
      confidence = NA_character_,
      markers = list(),
      reason = reason,
      skip = TRUE,
      set_at_step = state$step
    )
    write_decision_file(uuid)
    return(list(uuid = uuid, label = "(skipped)", reason = reason))
  }

  name <- args$name
  if (is.null(name) || !nzchar(name)) stop("name required (or use --skip)")
  conf <- toupper(args$confidence %||% "M")
  if (!(conf %in% c("H", "M", "L"))) stop("confidence must be H, M, or L")
  markers <- args$markers
  if (is.null(markers)) markers <- list()
  if (!is.list(markers)) markers <- as.list(markers)
  markers <- vapply(markers, as.character, character(1))
  markers <- trimws(markers)
  markers <- markers[nzchar(markers)]

  # Per-confidence minimum-marker rule (enforced at label time, fast feedback).
  min_markers <- c(H = 3L, M = 2L, L = 1L)[[conf]]
  if (length(markers) < min_markers) {
    stop(sprintf("confidence=%s requires at least %d evidence_markers; got %d",
                 conf, min_markers, length(markers)))
  }
  if (conf == "L" && !nzchar(reason)) {
    stop("confidence=L requires a non-empty --reason")
  }

  state$step <- state$step + 1L
  state$clusters[[idx]]$label <- list(
    name = name,
    confidence = conf,
    markers = as.list(markers),
    reason = reason,
    skip = FALSE,
    set_at_step = state$step
  )
  write_decision_file(uuid)
  list(uuid = uuid, name = name, confidence = conf,
       n_markers = length(markers), step = state$step)
}

op_unlabel <- function(args) {
  uuid <- args$uuid
  idx <- cluster_index_by_uuid(uuid)
  if (is.na(idx)) stop(sprintf("unknown cluster uuid: %s", uuid))
  if (is.null(state$clusters[[idx]]$label)) {
    return(list(uuid = uuid, removed = FALSE,
                msg = "cluster had no label"))
  }
  state$step <- state$step + 1L
  prev <- state$clusters[[idx]]$label
  state$clusters[[idx]]$label <- NULL
  write_decision_file(uuid)
  list(uuid = uuid, removed = TRUE, previous = prev, step = state$step)
}

run_qa_checks <- function() {
  findings <- list()  # each: list(severity = "error"|"warn", code, msg, uuid)

  alive <- alive_clusters()

  # Check 1: every alive cluster has a label or is skipped.
  for (c in alive) {
    if (is.null(c$label)) {
      findings[[length(findings) + 1L]] <- list(
        severity = "error", code = "missing_label",
        msg = sprintf("cluster %s has no label and is not skipped", c$uuid),
        uuid = c$uuid)
    }
  }

  # Check 2: marker citation. Every evidence marker must appear in the cluster's
  # cached top-50 (which came from FindMarkers vs rest or FindAllMarkers).
  for (c in alive) {
    if (is.null(c$label) || isTRUE(c$label$skip)) next
    cache_genes <- cached_gene_names(c$uuid)
    if (length(cache_genes) == 0L) {
      findings[[length(findings) + 1L]] <- list(
        severity = "error", code = "no_marker_cache",
        msg = sprintf("cluster %s labeled but has no cached markers; run `markers %s` first",
                      c$uuid, c$uuid),
        uuid = c$uuid)
      next
    }
    cited <- unlist(c$label$markers)
    missing <- setdiff(cited, cache_genes)
    if (length(missing) > 0) {
      findings[[length(findings) + 1L]] <- list(
        severity = "error", code = "hallucinated_markers",
        msg = sprintf("cluster %s cites markers not in its top-50: %s",
                      c$uuid, paste(missing, collapse = ",")),
        uuid = c$uuid,
        missing = as.list(missing))
    }
  }

  # Check 3: duplicate label names (warning only — sometimes legitimate, e.g.
  # CD8 memory vs CD8 effector both partially "CD8 T cell").
  by_name <- list()
  for (c in alive) {
    if (is.null(c$label) || isTRUE(c$label$skip)) next
    nm <- c$label$name
    by_name[[nm]] <- c(by_name[[nm]] %||% character(0), c$uuid)
  }
  for (nm in names(by_name)) {
    if (length(by_name[[nm]]) > 1) {
      findings[[length(findings) + 1L]] <- list(
        severity = "warn", code = "duplicate_label",
        msg = sprintf("label '%s' is used for %d clusters: %s",
                      nm, length(by_name[[nm]]),
                      paste(by_name[[nm]], collapse = ",")),
        uuids = as.list(by_name[[nm]]))
    }
  }

  n_err <- sum(vapply(findings, function(f) identical(f$severity, "error"),
                      logical(1)))
  n_warn <- sum(vapply(findings, function(f) identical(f$severity, "warn"),
                       logical(1)))

  list(
    pass = n_err == 0L,
    n_errors = n_err,
    n_warnings = n_warn,
    n_alive = length(alive),
    n_labeled = sum(vapply(alive, function(c) {
      !is.null(c$label) && !isTRUE(c$label$skip)
    }, logical(1))),
    n_skipped = sum(vapply(alive, function(c) {
      !is.null(c$label) && isTRUE(c$label$skip)
    }, logical(1))),
    findings = findings
  )
}

op_qa <- function(args) {
  require_seurat()
  run_qa_checks()
}

op_finalize <- function(args) {
  require_seurat()
  qa <- run_qa_checks()
  force <- isTRUE(args$force)
  if (!qa$pass && !force) {
    stop(sprintf("qa failed with %d error(s); use --force to override. First: %s",
                 qa$n_errors,
                 if (length(qa$findings) > 0) qa$findings[[1]]$msg else "?"))
  }
  out_dir <- outputs_dir()
  out_rds <- args$out_rds %||% file.path(out_dir, "annotated.rds")
  out_tsv <- args$out_tsv %||% file.path(out_dir, "annotation.tsv")
  out_md  <- args$out_md  %||% file.path(out_dir, "report.md")

  obj <- state$seurat
  # Build per-cell label vector.
  cell_idents <- as.character(Idents(obj))
  sid_to_label <- character(0)
  rows <- list()
  for (c in alive_clusters()) {
    sid <- as.character(c$seurat_id)
    if (is.null(c$label)) {
      lbl <- "[unresolved]"
    } else if (isTRUE(c$label$skip)) {
      lbl <- "[skipped]"
    } else {
      lbl <- c$label$name
    }
    sid_to_label[[sid]] <- lbl
    rows[[length(rows) + 1L]] <- data.frame(
      cluster_uuid = c$uuid,
      seurat_id = sid,
      n_cells = c$n_cells,
      label = lbl,
      confidence = if (is.null(c$label) || isTRUE(c$label$skip)) NA_character_ else c$label$confidence,
      markers = if (is.null(c$label) || isTRUE(c$label$skip)) NA_character_ else paste(unlist(c$label$markers), collapse = ","),
      reason = if (is.null(c$label)) NA_character_ else c$label$reason,
      skip = if (is.null(c$label)) FALSE else isTRUE(c$label$skip),
      stringsAsFactors = FALSE
    )
  }
  ann_table <- do.call(rbind, rows)
  cell_labels <- vapply(cell_idents,
                        function(s) sid_to_label[[s]] %||% "[unknown]",
                        character(1), USE.NAMES = FALSE)
  names(cell_labels) <- colnames(obj)
  obj <- AddMetaData(obj, metadata = cell_labels, col.name = "cassia_label")

  uuid_per_cell <- vapply(cell_idents, function(s) {
    u <- cluster_uuid_for_seurat_id(s)
    if (is.na(u)) s else u
  }, character(1), USE.NAMES = FALSE)
  names(uuid_per_cell) <- colnames(obj)
  obj <- AddMetaData(obj, metadata = uuid_per_cell, col.name = "cassia_cluster_uuid")

  saveRDS(obj, out_rds)
  write.table(ann_table, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

  # Markdown report
  md <- c(
    sprintf("# CASSIA annotation report"),
    sprintf(""),
    sprintf("- RDS: `%s`", state$rds_path),
    sprintf("- Cells: %d", ncol(obj)),
    sprintf("- Clusters (alive): %d", length(alive_clusters())),
    sprintf("- QA: %s (%d error(s), %d warning(s))",
            if (qa$pass) "PASS" else "FORCED",
            qa$n_errors, qa$n_warnings),
    sprintf(""),
    sprintf("## Labels"),
    sprintf("")
  )
  for (i in seq_len(nrow(ann_table))) {
    r <- ann_table[i, ]
    md <- c(md, sprintf("- **%s** (uuid `%s`, id %s, n=%d, conf=%s)  ",
                        r$label, r$cluster_uuid, r$seurat_id, r$n_cells,
                        r$confidence %||% "-"),
            sprintf("  markers: %s  ", r$markers %||% "-"),
            sprintf("  reason: %s", r$reason %||% "-"))
  }
  if (qa$n_warnings + qa$n_errors > 0) {
    md <- c(md, "", "## QA findings", "")
    for (f in qa$findings) {
      md <- c(md, sprintf("- [%s] %s", toupper(f$severity), f$msg))
    }
  }
  writeLines(md, out_md)

  state$phase <- "FINALIZED"
  list(
    qa = qa,
    out_rds = out_rds,
    out_tsv = out_tsv,
    out_md = out_md,
    n_labeled = qa$n_labeled,
    n_skipped = qa$n_skipped
  )
}

op_shutdown <- function(args) {
  log_msg("shutdown requested")
  state$shutdown_requested <- TRUE
  list(bye = TRUE)
}

dispatch_table <- list(
  ping = op_ping,
  init = op_init,
  status = op_status,
  checkpoint = op_checkpoint,
  restore = op_restore,
  policy = op_policy,
  export = op_export,
  clusters = op_clusters,
  markers = op_markers,
  gene = op_gene,
  preprocess = op_preprocess,
  merge = op_merge,
  subcluster = op_subcluster,
  label = op_label,
  unlabel = op_unlabel,
  qa = op_qa,
  finalize = op_finalize,
  shutdown = op_shutdown
)

`%||%` <- function(a, b) if (is.null(a)) b else a

.summarize_result <- function(op, result) {
  if (identical(op, "init")) {
    return(list(n_cells = result$n_cells,
                n_features = result$n_features,
                n_clusters_alive = result$n_clusters_alive,
                phase = result$phase))
  }
  if (identical(op, "status")) {
    return(list(phase = result$phase, step = result$step,
                n_clusters_alive = result$n_clusters_alive,
                seurat_mb = result$seurat_mb))
  }
  if (op %in% c("checkpoint", "restore")) {
    return(list(path = result$path, step = result$step,
                phase = result$phase,
                n_clusters_alive = result$n_clusters_alive))
  }
  if (identical(op, "export")) {
    return(list(partition_id = result$partition_id,
                n_cells = result$n_cells,
                n_clusters = result$n_clusters,
                memberships_md5 = result$memberships_md5))
  }
  if (identical(op, "policy")) {
    return(result)
  }
  if (identical(op, "clusters")) {
    return(list(n_returned = length(result)))
  }
  if (identical(op, "markers")) {
    rows <- result$markers
    if (is.data.frame(rows) && nrow(rows) > 0) {
      if (identical(result$scope, "all")) {
        # Per-cluster top gene list keyed by cluster id.
        by_cluster <- split(as.character(rows$gene), rows$cluster)
        top_genes_per_cluster <- lapply(by_cluster, function(g) head(g, 5))
        return(list(scope = "all", n_rows = nrow(rows),
                    top_genes = top_genes_per_cluster))
      }
      return(list(scope = "single",
                  ident_1 = result$ident_1,
                  ident_1_uuid = result$ident_1_uuid,
                  ident_2 = result$ident_2,
                  n_rows = nrow(rows),
                  top_genes = head(as.character(rows$gene), 5)))
    }
    return(list(scope = result$scope, n_rows = 0))
  }
  NULL
}

# ── one-shot transaction mode ───────────────────────────────────────────
#
# This is the default execution path used by the public CLI.  Every command
# starts a short-lived R process, restores an immutable input state, executes
# exactly one operation, and atomically writes a new state version when the
# operation changes scientific state.  The socket daemon below remains an
# optional acceleration for unusually large objects.

run_one_shot <- function(path) {
  request <- fromJSON(path, simplifyVector = FALSE)
  id <- request$id %||% "once"
  op <- request$op
  op_args <- request$args %||% list()
  input_state <- request$input_state
  output_state <- request$output_state
  fn <- dispatch_table[[op]]

  if (is.null(fn)) {
    response <- list(id = id, ok = FALSE, elapsed_s = 0,
                     error = paste("unknown op:", op))
    writeLines(toJSON(response, auto_unbox = TRUE, null = "null", na = "null",
                      force = TRUE, digits = 6), stdout())
    return(invisible(FALSE))
  }

  t0 <- Sys.time()
  execution <- tryCatch({
    if (!identical(op, "init") && !identical(op, "restore")) {
      if (is.null(input_state) || !nzchar(input_state)) {
        stop(sprintf("operation '%s' requires an input state", op))
      }
      op_restore(list(path = input_state))
    }
    value <- fn(op_args)
    if (!is.null(output_state) && nzchar(output_state)) {
      atomic_save_rds(checkpoint_payload(), output_state)
    }
    value
  }, error = function(e) e)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (inherits(execution, "error")) {
    response <- list(id = id, ok = FALSE, elapsed_s = elapsed,
                     error = conditionMessage(execution))
    summary_payload <- NULL
  } else {
    response <- list(id = id, ok = TRUE, elapsed_s = elapsed,
                     result = execution,
                     output_state = output_state %||% NA_character_)
    summary_payload <- tryCatch(.summarize_result(op, execution),
                                error = function(e) NULL)
  }

  audit_append(list(event = "op", id = id, op = op, args = op_args,
                    execution_mode = "one-shot", ok = response$ok,
                    elapsed_s = round(elapsed, 4),
                    ts = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3"),
                    summary = summary_payload,
                    error = if (!isTRUE(response$ok)) response$error else NULL))
  writeLines(toJSON(response, auto_unbox = TRUE, null = "null", na = "null",
                    force = TRUE, digits = 6), stdout())
  invisible(isTRUE(response$ok))
}

if (identical(run_mode, "once")) {
  run_one_shot(request_path)
  quit(save = "no", status = 0L, runLast = FALSE)
}

# ── main loop ────────────────────────────────────────────────────────────

srv <- serverSocket(port = port)
pid <- Sys.getpid()
log_msg(sprintf("daemon up: pid=%d port=%d workdir=%s", pid, port, workdir))

# Touch a ready file so the Python launcher knows the bind succeeded.
writeLines(as.character(pid), ready_path)

audit_append(list(event = "daemon_start", pid = pid, port = port,
                  ts = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))

handle_one <- function(con) {
  on.exit(try(close(con), silent = TRUE), add = TRUE)
  req_line <- readLines(con, n = 1, warn = FALSE)
  if (length(req_line) == 0 || !nzchar(req_line)) return(invisible(NULL))
  parsed <- tryCatch(fromJSON(req_line, simplifyVector = FALSE),
                     error = function(e) e)
  if (inherits(parsed, "error")) {
    resp <- list(id = "?", ok = FALSE, elapsed_s = 0,
                 error = paste("bad json:", conditionMessage(parsed)))
    writeLines(toJSON(resp, auto_unbox = TRUE, null = "null", na = "null",
                      force = TRUE, digits = 6), con)
    return(invisible(NULL))
  }
  id <- parsed$id %||% "?"
  op <- parsed$op
  args <- parsed$args
  fn <- dispatch_table[[op]]
  if (is.null(fn)) {
    resp <- list(id = id, ok = FALSE, elapsed_s = 0,
                 error = paste("unknown op:", op))
    writeLines(toJSON(resp, auto_unbox = TRUE, null = "null", na = "null",
                      force = TRUE, digits = 6), con)
    return(invisible(NULL))
  }
  t0 <- Sys.time()
  result <- tryCatch(fn(args), error = function(e) e)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(result, "error")) {
    resp <- list(id = id, ok = FALSE, elapsed_s = elapsed,
                 error = conditionMessage(result))
  } else {
    resp <- list(id = id, ok = TRUE, elapsed_s = elapsed, result = result)
  }
  # Build an op-specific result summary that captures the agent-visible answer
  # without dumping the full payload (which can be 80+ rows for FindAllMarkers).
  summary_payload <- NULL
  if (isTRUE(resp$ok)) {
    summary_payload <- tryCatch(.summarize_result(op, result),
                                error = function(e) NULL)
  }
  audit_append(list(event = "op", id = id, op = op,
                    args = args,
                    ok = resp$ok,
                    elapsed_s = round(elapsed, 4),
                    ts = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3"),
                    summary = summary_payload,
                    error = if (!isTRUE(resp$ok)) resp$error else NULL))
  writeLines(toJSON(resp, auto_unbox = TRUE, null = "null", na = "null",
                    force = TRUE, digits = 6), con)
}

state$shutdown_requested <- FALSE
repeat {
  con <- tryCatch(socketAccept(srv, blocking = TRUE, open = "r+",
                               timeout = 60),
                  error = function(e) e)
  if (inherits(con, "error")) {
    log_msg("accept error:", conditionMessage(con))
    next
  }
  if (is.null(con)) next  # accept timeout, loop and check shutdown flag
  handle_one(con)
  if (isTRUE(state$shutdown_requested)) {
    log_msg("exiting main loop")
    break
  }
}

try(close(srv), silent = TRUE)
try(file.remove(ready_path), silent = TRUE)
audit_append(list(event = "daemon_stop", pid = pid,
                  ts = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
log_msg("daemon exited cleanly")
