#!/usr/bin/env Rscript
# Test script for the new symphonyCompare function in R

cat("🎼 CASSIA Symphony Compare - R Package Test\n")
cat("===========================================\n\n")

# Try to load the CASSIA package
tryCatch({
  library(CASSIA)
  cat("✅ Successfully loaded CASSIA package\n")
}, error = function(e) {
  cat("❌ Failed to load CASSIA package:", e$message, "\n")
  cat("💡 Make sure the package is installed and Python environment is set up\n")
  quit(status = 1)
})

# Check if symphonyCompare function is available
if (exists("symphonyCompare")) {
  cat("✅ symphonyCompare function is available\n")
} else {
  cat("❌ symphonyCompare function not found\n")
  quit(status = 1)
}

# Display function information
cat("\n📖 Function Documentation:\n")
cat("==========================\n")

# Show function signature
cat("Function signature:\n")
print(args(symphonyCompare))

cat("\n📝 Example Usage (without API key):\n")
cat("===================================\n")

example_code <- '
# Basic usage - let Symphony Compare handle everything
results <- symphonyCompare(
  tissue = "peripheral blood",
  celltypes = c("T cell", "B cell", "NK cell", "Monocyte"),
  marker_set = c("CD3", "CD4", "CD8", "CD19", "CD20", "CD16", "CD56", "CD14"),
  species = "human"
)

# Access the results
cat("Consensus:", results$consensus, "\\n")
cat("Confidence:", sprintf("%.1f%%", results$confidence * 100), "\\n")

# Advanced usage with custom settings
results <- symphonyCompare(
  tissue = "brain cortex",
  celltypes = c("Neuron", "Astrocyte", "Microglia"),
  marker_set = c("RBFOX3", "GFAP", "IBA1", "MAP2", "S100B", "P2RY12"),
  species = "human",
  model_preset = "quartet",  # Use 4 models instead of 3
  enable_discussion = TRUE,  # Enable automatic discussion rounds
  max_discussion_rounds = 3,  # Allow up to 3 discussion rounds
  consensus_threshold = 0.75, # 75% of models must agree
  output_dir = "./symphony_results",
  verbose = TRUE
)

# Work with the results dataframe
df <- results$dataframe
head(df[, c("model", "researcher", "round", "Neuron_score", "Astrocyte_score")])

# Access summary statistics
summary_stats <- results$summary
cat("Total rounds:", summary_stats$total_rounds, "\\n")
cat("Models used:", summary_stats$models_used, "\\n")
'

cat(example_code)

# Check for API key
api_key <- Sys.getenv("OPENROUTER_API_KEY")
if (nchar(api_key) == 0) {
  cat("\n⚠️  OPENROUTER_API_KEY not set in environment\n")
  cat("To test with real API calls, set your API key:\n")
  cat('Sys.setenv(OPENROUTER_API_KEY = "your-key-here")\n')
  cat("or use: setOpenRouterApiKey(\"your-key-here\")\n\n")
} else {
  cat("\n🔑 API key detected - ready for live testing!\n")
  cat("You can now run symphonyCompare with real API calls.\n\n")
  
  # Optionally run a quick test (commented out to avoid API usage)
  cat("# Example test (uncomment to run):\n")
  cat("# results <- symphonyCompare(\n")
  cat("#   tissue = \"blood\",\n")
  cat("#   celltypes = c(\"T cell\", \"B cell\"),\n")
  cat("#   marker_set = c(\"CD3\", \"CD19\"),\n")
  cat("#   model_preset = \"budget\",  # Use budget models for testing\n")
  cat("#   verbose = TRUE\n")
  cat("# )\n")
}

cat("\n📊 Feature Comparison:\n")
cat("======================\n")
comparison_table <- '
┌─────────────────────────┬──────────────────────┬──────────────────────┐
│ Feature                 │ compareCelltypes     │ symphonyCompare      │
├─────────────────────────┼──────────────────────┼──────────────────────┤
│ Function Name           │ compareCelltypes     │ symphonyCompare ✨   │
│ Multi-model Analysis    │ ✅ Yes               │ ✅ Yes               │
│ Parallel Processing     │ ✅ Yes               │ ✅ Yes               │
│ Auto Discussion Rounds  │ ✅ Yes (optional)    │ ✅ Yes (default on)  │
│ Consensus Detection     │ ✅ Basic             │ ✅ Advanced          │
│ Confidence Scores       │ ❌ No                │ ✅ Yes               │
│ Model Presets           │ ✅ 2 presets         │ ✅ 3+ presets        │
│ Custom Models           │ ✅ Yes               │ ✅ Yes               │
│ HTML Reports            │ ✅ Yes               │ ✅ Enhanced          │
│ CSV Output              │ ✅ Yes               │ ✅ Yes               │
│ Progress Tracking       │ ❌ Basic             │ ✅ Detailed          │
│ Summary Statistics      │ ❌ No                │ ✅ Yes               │
│ Output Organization     │ ❌ Basic             │ ✅ Advanced          │
│ R Integration           │ ✅ Good              │ ✅ Enhanced          │
└─────────────────────────┴──────────────────────┴──────────────────────┘
'
cat(comparison_table)

cat("\n✨ Key Benefits of symphonyCompare:\n")
cat("===================================\n")
cat("• 🎼 Elegant musical naming theme (Symphony, Movements, etc.)\n")
cat("• 📊 Returns confidence scores and detailed statistics\n")
cat("• 🗂️ Better output file organization with timestamps\n")
cat("• 🎯 Configurable consensus threshold (default 80%)\n")
cat("• 📈 Detailed cell type score statistics (mean, std, range)\n")
cat("• 🎨 Enhanced progress messages and formatting\n")
cat("• 🔄 Automatic consensus building with discussion rounds\n")
cat("• 📋 Rich R data.frame output for further analysis\n")

cat("\n🎵 Symphony Compare is ready to orchestrate your cell type analysis!\n")
cat("Use ?symphonyCompare in R for detailed documentation.\n")