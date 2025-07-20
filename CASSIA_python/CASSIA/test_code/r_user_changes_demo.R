# CASSIA R Package: Model Settings Changes for Users
# ===================================================

# This script demonstrates what changes for R users with the new model settings system

library(CASSIA)

# =============================================================================
# BEFORE (Old Way) - Complex Model Names
# =============================================================================

# Old way: Users had to remember complex, provider-specific model names
old_way_examples <- list(
  # OpenAI models
  'runCASSIA_batch(marker = my_markers, model = "gpt-4o", provider = "openai")',
  'runCASSIA_batch(marker = my_markers, model = "gpt-4o-mini", provider = "openai")',
  
  # Anthropic models  
  'runCASSIA_batch(marker = my_markers, model = "claude-3-5-sonnet-latest", provider = "anthropic")',
  'runCASSIA_batch(marker = my_markers, model = "claude-3-5-haiku-latest", provider = "anthropic")',
  
  # OpenRouter models
  'runCASSIA_batch(marker = my_markers, model = "google/gemini-2.5-flash", provider = "openrouter")',
  'runCASSIA_batch(marker = my_markers, model = "deepseek/deepseek-chat-v3-0324", provider = "openrouter")',
  'runCASSIA_batch(marker = my_markers, model = "anthropic/claude-3-5-sonnet", provider = "openrouter")'
)

cat("=== OLD WAY (Complex Names) ===\n")
for (example in old_way_examples) {
  cat(example, "\n")
}

# =============================================================================
# AFTER (New Way) - Simple Model Names + Shortcuts
# =============================================================================

cat("\n=== NEW WAY (Simple Names) ===\n")

# New way: Simple, memorable names that work across providers
new_way_examples <- list(
  # Simple names
  'runCASSIA_batch(marker = my_markers, model = "gpt4", provider = "openai")',
  'runCASSIA_batch(marker = my_markers, model = "claude", provider = "anthropic")', 
  'runCASSIA_batch(marker = my_markers, model = "gemini", provider = "openrouter")',
  'runCASSIA_batch(marker = my_markers, model = "deepseek", provider = "openrouter")',
  
  # Quality shortcuts - different models per provider!
  'runCASSIA_batch(marker = my_markers, model = "best", provider = "openai")',     # → gpt-4o
  'runCASSIA_batch(marker = my_markers, model = "best", provider = "anthropic")',  # → claude-3-5-sonnet-latest
  'runCASSIA_batch(marker = my_markers, model = "best", provider = "openrouter")', # → google/gemini-2.5-flash
  
  'runCASSIA_batch(marker = my_markers, model = "cheap", provider = "openai")',     # → gpt-3.5-turbo
  'runCASSIA_batch(marker = my_markers, model = "cheap", provider = "anthropic")',  # → claude-3-5-haiku-latest  
  'runCASSIA_batch(marker = my_markers, model = "cheap", provider = "openrouter")', # → deepseek/deepseek-chat-v3-0324
  
  'runCASSIA_batch(marker = my_markers, model = "fast", provider = "openrouter")'   # → google/gemini-2.5-flash
)

for (example in new_way_examples) {
  cat(example, "\n")
}

# =============================================================================
# NEW R HELPER FUNCTIONS
# =============================================================================

cat("\n=== NEW R HELPER FUNCTIONS ===\n")

# 1. Check what a simple name resolves to
cat("# Check model resolution:\n")
cat('resolve_model_name("gemini", "openrouter")\n')
# Result: list(model = "google/gemini-2.5-flash", provider = "openrouter")

cat('resolve_model_name("claude")\n') 
# Will auto-detect provider from model name patterns

# 2. Get recommended models
cat("\n# Get recommendations:\n")
cat('get_recommended_model(provider = "openai")\n')
cat('get_recommended_model(provider = "anthropic")\n')
cat('get_recommended_model(use_case = "annotation")\n')

# 3. List available models
cat("\n# List models:\n")
cat('list_models()\n')                                    # All models
cat('list_models(provider = "openrouter")\n')             # By provider
cat('list_models(cost_tier = "low")\n')                   # By cost
cat('list_models(provider = "openrouter", cost_tier = "low")\n')  # Combined filters

# 4. Get model information
cat("\n# Get model info:\n")
cat('get_model_info("gemini", "openrouter")\n')
cat('get_model_info("deepseek")\n')

# 5. Print recommendations
cat("\n# Print recommendations:\n")
cat('print_model_recommendations()\n')                    # All recommendations
cat('print_model_recommendations("annotation")\n')       # For specific use case

# =============================================================================
# PRACTICAL EXAMPLES FOR R USERS
# =============================================================================

cat("\n=== PRACTICAL R EXAMPLES ===\n")

practical_examples <- c(
  "# Load your marker data",
  "markers <- read.csv('your_markers.csv')",
  "",
  "# 1. Quick analysis with simple names:",
  "runCASSIA_batch(marker = markers, model = 'gemini', provider = 'openrouter')",
  "",
  "# 2. Cost-conscious analysis:",
  "runCASSIA_batch(marker = markers, model = 'cheap', provider = 'openrouter')",
  "",
  "# 3. High-quality analysis:",
  "runCASSIA_batch(marker = markers, model = 'best', provider = 'anthropic')",
  "",
  "# 4. Multi-step pipeline with different models:",
  "runCASSIA_pipeline(",
  "  marker = markers,",
  "  annotation_model = 'gemini',",
  "  annotation_provider = 'openrouter',",
  "  score_model = 'claude', ", 
  "  score_provider = 'anthropic',",
  "  annotationboost_model = 'best',",
  "  annotationboost_provider = 'openai'",
  ")",
  "",
  "# 5. Check what models you're actually using:",
  "resolve_model_name('gemini', 'openrouter')  # Shows: google/gemini-2.5-flash",
  "resolve_model_name('cheap', 'openrouter')   # Shows: deepseek/deepseek-chat-v3-0324"
)

for (example in practical_examples) {
  cat(example, "\n")
}

# =============================================================================
# SUMMARY OF CHANGES FOR R USERS
# =============================================================================

cat("\n=== SUMMARY: WHAT CHANGES FOR R USERS ===\n")

changes <- c(
  "✅ EASIER MODEL SELECTION:",
  "   - Use 'gpt4', 'claude', 'gemini' instead of complex names",
  "   - Quality shortcuts: 'best', 'cheap', 'premium', 'fast'",
  "",
  "✅ PROVIDER CONTROL:",
  "   - Same shortcut gives different models per provider",
  "   - 'cheap' + 'openai' → gpt-3.5-turbo", 
  "   - 'cheap' + 'openrouter' → deepseek (much cheaper!)",
  "",
  "✅ NEW R HELPER FUNCTIONS:",
  "   - resolve_model_name() - see what you're actually using",
  "   - get_recommended_model() - get best model for provider/use case",
  "   - list_models() - browse available options",
  "   - get_model_info() - detailed model information",
  "",
  "✅ BACKWARD COMPATIBILITY:",
  "   - All existing code continues to work",
  "   - Can mix old and new approaches",
  "   - Gradual migration possible",
  "",
  "✅ COST CONTROL:",
  "   - Easy to switch to cheaper alternatives",
  "   - Clear cost tier information",
  "   - Provider-specific cost optimization",
  "",
  "⚠️  PROVIDER STILL REQUIRED:",
  "   - Must specify provider for API key control",
  "   - Prevents accidental API switching",
  "   - Ensures you use the right API key"
)

for (change in changes) {
  cat(change, "\n")
}

cat("\n🎉 The model settings system makes CASSIA easier to use while giving you more control!\n")