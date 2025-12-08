pkgname <- "CASSIA"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "CASSIA-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('CASSIA')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("calculate_evaluation_metrics")
### * calculate_evaluation_metrics

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calculate_evaluation_metrics
### Title: Calculate Evaluation Metrics
### Aliases: calculate_evaluation_metrics

### ** Examples

## Not run: 
##D # Calculate metrics from evaluation results
##D metrics <- calculate_evaluation_metrics(
##D   eval_df = evaluation_results,
##D   score_col = "evaluation_score"
##D )
##D 
##D # Print key metrics
##D cat("Mean Score:", metrics$mean_score, "\n")
##D cat("Perfect Predictions:", sprintf("%.1f%%", metrics$perfect_ratio * 100), "\n")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calculate_evaluation_metrics", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("generate_html_report")
### * generate_html_report

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: generate_html_report
### Title: Generate HTML Evaluation Report
### Aliases: generate_html_report

### ** Examples

## Not run: 
##D # Generate evaluation report
##D generate_html_report(
##D   result_df = evaluation_results,
##D   gold_col = "Cluster ID",
##D   pred_col = "Predicted General Cell Type",
##D   score_col = "evaluation_score",
##D   reasoning_col = "explanation",
##D   html_report_path = "evaluation_report.html",
##D   model_name = "Claude-3.5-Sonnet"
##D )
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("generate_html_report", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("generate_subclustering_report")
### * generate_subclustering_report

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: generate_subclustering_report
### Title: Generate Subclustering HTML Report
### Aliases: generate_subclustering_report

### ** Examples

## Not run: 
##D # Generate a subclustering report from CSV results
##D generate_subclustering_report(
##D   csv_path = "subclustering_results.csv",
##D   html_report_path = "subclustering_report.html",
##D   model_name = "Gemini-2.5-Flash"
##D )
##D 
##D # Auto-generate HTML filename and model name
##D generate_subclustering_report("subclustering_results.csv")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("generate_subclustering_report", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_cost_tier_models")
### * get_cost_tier_models

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_cost_tier_models
### Title: Get Cost Tier Models
### Aliases: get_cost_tier_models

### ** Examples

## Not run: 
##D # Get models by cost tier
##D get_cost_tier_models("very_low")
##D get_cost_tier_models("low")
##D get_cost_tier_models("medium")
##D get_cost_tier_models("high")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_cost_tier_models", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_logger")
### * get_logger

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_logger
### Title: Get CASSIA Logger
### Aliases: get_logger

### ** Examples

## Not run: 
##D # Get the default CASSIA logger
##D logger <- get_logger()
##D 
##D # Get a named logger for custom module
##D my_logger <- get_logger("my_analysis")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_logger", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_model_aliases")
### * get_model_aliases

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_model_aliases
### Title: Get Model Aliases
### Aliases: get_model_aliases

### ** Examples

## Not run: 
##D # Get aliases for models
##D get_model_aliases("gpt-4o")
##D get_model_aliases("claude-3-5-sonnet-latest")
##D get_model_aliases("google/gemini-2.5-flash")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_model_aliases", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_model_info")
### * get_model_info

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_model_info
### Title: Get Model Information
### Aliases: get_model_info

### ** Examples

## Not run: 
##D # Get info for specific models
##D get_model_info("gpt-4o")
##D get_model_info("claude-3-5-sonnet-latest")
##D get_model_info("google/gemini-2.5-flash")
##D 
##D # Using aliases
##D get_model_info("gemini")
##D get_model_info("deepseek")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_model_info", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_recommended_model")
### * get_recommended_model

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_recommended_model
### Title: Get Recommended Model
### Aliases: get_recommended_model

### ** Examples

## Not run: 
##D # Get overall best model
##D get_recommended_model()
##D 
##D # Get best model for specific use case
##D get_recommended_model(use_case = "annotation")
##D get_recommended_model(use_case = "scoring")
##D get_recommended_model(use_case = "annotation_boost")
##D 
##D # Get recommended model for specific provider
##D get_recommended_model(provider = "openai")
##D get_recommended_model(provider = "anthropic")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_recommended_model", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_use_case_recommendations")
### * get_use_case_recommendations

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_use_case_recommendations
### Title: Get Use Case Recommendations
### Aliases: get_use_case_recommendations

### ** Examples

## Not run: 
##D # Get recommendations for different use cases
##D get_use_case_recommendations("annotation")
##D get_use_case_recommendations("scoring")
##D get_use_case_recommendations("annotation_boost")
##D get_use_case_recommendations("merging")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_use_case_recommendations", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("listAvailableMarkers")
### * listAvailableMarkers

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: listAvailableMarkers
### Title: List Available Built-in Marker Sets
### Aliases: listAvailableMarkers

### ** Examples

## Not run: 
##D available_markers <- listAvailableMarkers()
##D print(available_markers)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("listAvailableMarkers", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("list_models")
### * list_models

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: list_models
### Title: List Available Models
### Aliases: list_models

### ** Examples

## Not run: 
##D # List all models
##D list_models()
##D 
##D # Filter by provider
##D list_models(provider = "openai")
##D list_models(provider = "anthropic")
##D list_models(provider = "openrouter")
##D 
##D # Filter by cost tier
##D list_models(cost_tier = "low")
##D list_models(cost_tier = "high")
##D 
##D # Filter by use case
##D list_models(use_case = "annotation")
##D list_models(use_case = "scoring")
##D 
##D # Combine filters
##D list_models(provider = "openrouter", cost_tier = "low")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("list_models", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("loadBuiltinMarkers")
### * loadBuiltinMarkers

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: loadBuiltinMarkers
### Title: Load Built-in Marker Data
### Aliases: loadBuiltinMarkers

### ** Examples

## Not run: 
##D markers <- loadBuiltinMarkers()
##D head(markers)
##D 
##D subcluster_results <- loadBuiltinMarkers("subcluster_results")
##D head(subcluster_results)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("loadBuiltinMarkers", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("loadExampleMarkers")
### * loadExampleMarkers

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: loadExampleMarkers
### Title: Load Example Marker Data
### Aliases: loadExampleMarkers

### ** Examples

markers <- loadExampleMarkers()
head(markers)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("loadExampleMarkers", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("loadExampleMarkers_subcluster")
### * loadExampleMarkers_subcluster

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: loadExampleMarkers_subcluster
### Title: Load Example Subcluster Results
### Aliases: loadExampleMarkers_subcluster

### ** Examples

subcluster_data <- loadExampleMarkers_subcluster()
head(subcluster_data)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("loadExampleMarkers_subcluster", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("merge_annotations")
### * merge_annotations

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: merge_annotations
### Title: Merge Cell Cluster Annotations
### Aliases: merge_annotations

### ** Examples

## Not run: 
##D # Basic usage - merge annotations with broad groupings
##D result <- merge_annotations("annotation_results.csv", "merged_results.csv")
##D 
##D # Use detailed groupings
##D result <- merge_annotations(
##D   csv_path = "annotation_results.csv",
##D   output_path = "merged_detailed.csv",
##D   detail_level = "detailed",
##D   provider = "anthropic"
##D )
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("merge_annotations", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("merge_annotations_all")
### * merge_annotations_all

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: merge_annotations_all
### Title: Merge All Annotation Levels
### Aliases: merge_annotations_all

### ** Examples

## Not run: 
##D # Process all grouping levels at once
##D result <- merge_annotations_all(
##D   csv_path = "annotation_results.csv",
##D   output_path = "merged_all_levels.csv",
##D   provider = "openrouter"
##D )
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("merge_annotations_all", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("print_model_recommendations")
### * print_model_recommendations

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: print_model_recommendations
### Title: Print Model Recommendations
### Aliases: print_model_recommendations

### ** Examples

## Not run: 
##D # Print all recommendations
##D print_model_recommendations()
##D 
##D # Print recommendations for specific use case
##D print_model_recommendations("annotation")
##D print_model_recommendations("scoring")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("print_model_recommendations", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("resolve_model_name")
### * resolve_model_name

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: resolve_model_name
### Title: Resolve Model Name
### Aliases: resolve_model_name

### ** Examples

## Not run: 
##D # Using simple names
##D resolve_model_name("gpt4")  # Returns gpt-4o with openai provider
##D resolve_model_name("claude")  # Returns claude-3-5-sonnet-latest with anthropic provider
##D resolve_model_name("gemini")  # Returns google/gemini-2.5-flash with openrouter provider
##D 
##D # Using aliases
##D resolve_model_name("sonnet")  # Returns claude-3-5-sonnet-latest
##D resolve_model_name("deepseek")  # Returns deepseek/deepseek-chat-v3-0324
##D 
##D # With specific provider
##D resolve_model_name("gpt-4o", "openai")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("resolve_model_name", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("runCASSIA_generate_score_report")
### * runCASSIA_generate_score_report

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: runCASSIA_generate_score_report
### Title: Generate HTML Reports from Scored Results
### Aliases: runCASSIA_generate_score_report

### ** Examples

## Not run: 
##D runCASSIA_generate_score_report("path/to/scored_results.csv")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("runCASSIA_generate_score_report", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("set_log_level")
### * set_log_level

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: set_log_level
### Title: Set CASSIA Log Level
### Aliases: set_log_level

### ** Examples

## Not run: 
##D # Enable verbose debugging
##D set_log_level("DEBUG")
##D 
##D # Show only errors (quiet mode)
##D set_log_level("ERROR")
##D 
##D # Suppress all messages
##D set_log_level("QUIET")
##D 
##D # Reset to default
##D set_log_level("INFO")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("set_log_level", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("symphonyCompare")
### * symphonyCompare

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: symphonyCompare
### Title: Symphony Compare - Advanced Multi-Model Cell Type Comparison
###   with Consensus Building
### Aliases: symphonyCompare

### ** Examples

## Not run: 
##D # Basic usage - let Symphony Compare handle everything
##D results <- symphonyCompare(
##D   tissue = "peripheral blood",
##D   celltypes = c("T cell", "B cell", "NK cell", "Monocyte"),
##D   marker_set = c("CD3", "CD4", "CD8", "CD19", "CD20", "CD16", "CD56", "CD14"),
##D   species = "human"
##D )
##D 
##D # Access the results
##D cat("Consensus:", results$consensus, "\n")
##D cat("Confidence:", sprintf("%.1f%%", results$confidence * 100), "\n")
##D 
##D # Advanced usage with custom settings
##D results <- symphonyCompare(
##D   tissue = "brain",
##D   celltypes = c("Neuron", "Astrocyte", "Microglia", "Oligodendrocyte"),
##D   marker_set = c("RBFOX3", "GFAP", "IBA1", "OLIG2", "MAP2", "S100B", "CD11B", "MBP"),
##D   species = "mouse",
##D   model_preset = "quartet",  # Use 4 models instead of 3
##D   enable_discussion = TRUE,  # Enable automatic discussion rounds
##D   max_discussion_rounds = 3,  # Allow up to 3 discussion rounds
##D   consensus_threshold = 0.75,  # 75% of models must agree
##D   output_dir = "./symphony_results",
##D   verbose = TRUE
##D )
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("symphonyCompare", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("warn_user")
### * warn_user

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: warn_user
### Title: Show a Warning to the User
### Aliases: warn_user

### ** Examples

## Not run: 
##D warn_user("This operation may take a long time with large datasets.")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("warn_user", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
