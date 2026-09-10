# Topology-policy primitives shared by one-shot and optional daemon runtimes.

cassia_value_or <- function(value, fallback) {
  if (is.null(value)) fallback else value
}

cassia_default_policy <- function() {
  list(
    enabled = FALSE,
    allowed_topology_actions = c("merge", "subcluster"),
    max_topology_edits = 0L,
    topology_edits = 0L,
    min_child_cells = 1L,
    max_children_per_split = .Machine$integer.max,
    lock_preprocess = FALSE
  )
}

cassia_resolve_policy <- function(args, current_policy) {
  enabled <- isTRUE(args$enabled)
  allowed <- cassia_value_or(args$allowed_topology_actions, list())
  allowed <- unique(as.character(unlist(allowed)))
  invalid <- setdiff(allowed, c("merge", "subcluster"))
  if (length(invalid) > 0L) {
    stop(sprintf("unknown topology actions: %s", paste(invalid, collapse = ", ")))
  }
  maximum <- as.integer(cassia_value_or(args$max_topology_edits, 0L))
  if (is.na(maximum) || maximum < 0L) {
    stop("max_topology_edits must be non-negative")
  }
  min_child_cells <- as.integer(cassia_value_or(args$min_child_cells, 1L))
  max_children <- as.integer(cassia_value_or(
    args$max_children_per_split, .Machine$integer.max
  ))
  if (is.na(min_child_cells) || min_child_cells < 1L) {
    stop("min_child_cells must be at least 1")
  }
  if (is.na(max_children) || max_children < 2L) {
    stop("max_children_per_split must be at least 2")
  }
  previous_used <- as.integer(cassia_value_or(current_policy$topology_edits, 0L))
  list(
    enabled = enabled,
    allowed_topology_actions = allowed,
    max_topology_edits = maximum,
    topology_edits = if (isTRUE(args$reset_counter)) 0L else previous_used,
    min_child_cells = min_child_cells,
    max_children_per_split = max_children,
    lock_preprocess = isTRUE(args$lock_preprocess)
  )
}

cassia_enforce_topology_policy <- function(policy, action) {
  if (is.null(policy) || !isTRUE(policy$enabled)) return(invisible(NULL))
  allowed <- as.character(unlist(cassia_value_or(
    policy$allowed_topology_actions, list()
  )))
  if (!(action %in% allowed)) {
    stop(sprintf("automation policy forbids topology action '%s'", action))
  }
  used <- as.integer(cassia_value_or(policy$topology_edits, 0L))
  maximum <- as.integer(cassia_value_or(policy$max_topology_edits, 0L))
  if (used >= maximum) {
    stop(sprintf("automation topology-edit budget exhausted (%d/%d)", used, maximum))
  }
  invisible(NULL)
}

cassia_record_topology_edit <- function(policy) {
  if (isTRUE(policy$enabled)) {
    policy$topology_edits <- as.integer(cassia_value_or(
      policy$topology_edits, 0L
    )) + 1L
  }
  policy
}

cassia_validate_split_policy <- function(child_sizes, policy) {
  min_child_cells <- as.integer(cassia_value_or(policy$min_child_cells, 1L))
  max_children <- as.integer(cassia_value_or(
    policy$max_children_per_split, .Machine$integer.max
  ))
  if (length(child_sizes) > max_children) {
    stop(sprintf(
      paste0("subcluster produced %d children (sizes: %s), above the policy ",
             "maximum of %d; transaction rolled back and topology budget was ",
             "not consumed. Try a lower resolution"),
      length(child_sizes), paste(child_sizes, collapse = ","), max_children
    ))
  }
  if (any(child_sizes < min_child_cells)) {
    stop(sprintf(
      paste0("subcluster child sizes %s violate the policy minimum of %d cells; ",
             "transaction rolled back and topology budget was not consumed. ",
             "Try a lower resolution"),
      paste(child_sizes, collapse = ","), min_child_cells
    ))
  }
  invisible(NULL)
}
