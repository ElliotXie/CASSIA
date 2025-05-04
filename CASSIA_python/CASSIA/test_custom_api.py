from CASSIA.tools_function import set_custom_provider, runCASSIA, list_custom_providers

# Test data
deepseek_api_key = "sk-e967ebf66bc24069a5dfa642792bc491"
deepseek_model = "deepseek-chat"
deepseek_url = "https://api.deepseek.com"

# Sample marker list for testing
test_markers = ["CD3D", "CD4", "CD8A", "FOXP3", "IL2RA"]

# 1. Set up the DeepSeek provider
print("Setting up DeepSeek provider...")
set_custom_provider(deepseek_api_key, "deepseek", deepseek_url)

# 2. List configured providers to verify setup
providers = list_custom_providers()
print(f"Configured providers: {providers}")

# 3. Run a simple analysis using DeepSeek
print("\nRunning analysis with DeepSeek...")
try:
    result, conversation = runCASSIA(
        model=deepseek_model,
        provider="deepseek",
        marker_list=test_markers,
        tissue="blood",
        species="human",
        temperature=0.2  # Adding some temperature for variation
    )
    
    # 4. Print results
    print("\nAnalysis complete!")
    print(f"Main cell type: {result.get('main_cell_type', 'Not found')}")
    print(f"Sub cell types: {result.get('sub_cell_types', [])}")
    print(f"Number of conversation turns: {len(conversation) if conversation else 0}")
    
except Exception as e:
    print(f"Error running analysis: {str(e)}")

# 5. Test another provider without resetting keys
print("\nSetting up a second provider for testing...")
set_custom_provider("dummy-key", "test_provider", "https://example.com/api")

# 6. List providers again to verify both are stored
providers = list_custom_providers()
print(f"Updated providers: {providers}")

# 7. Verify we can still use DeepSeek without resetting the key
print("\nRunning second analysis with DeepSeek to verify key persistence...")
try:
    result2, _ = runCASSIA(
        model=deepseek_model,
        provider="deepseek",
        marker_list=test_markers[:2],  # Using fewer markers for a quicker test
        tissue="lung",
        species="mouse",
        temperature=0
    )
    
    print(f"Second analysis main cell type: {result2.get('main_cell_type', 'Not found')}")
    
except Exception as e:
    print(f"Error running second analysis: {str(e)}")