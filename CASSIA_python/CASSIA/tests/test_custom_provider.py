"""
This test script verifies that the circular import issue has been fixed
and that the custom provider functionality works correctly.
"""

from CASSIA import (
    set_custom_provider, 
    list_custom_providers,
    get_custom_base_url,
    get_custom_api_key,
    set_custom_base_url,
    set_custom_api_key
)

print("1. Testing custom provider setup...")

# Test setting up multiple custom providers
providers = {
    "deepseek": {
        "api_key": "deepseek-api-key-1234",
        "base_url": "https://api.deepseek.com"
    },
    "anthropic_proxy": {
        "api_key": "proxy-key-5678",
        "base_url": "https://my-anthropic-proxy.example.com"
    }
}

for provider, config in providers.items():
    set_custom_provider(
        api_key=config["api_key"],
        provider_name=provider,
        base_url=config["base_url"]
    )
    print(f"  - Set up provider '{provider}'")

print("\n2. Testing list_custom_providers()...")
configured_providers = list_custom_providers()
print(f"  - Configured providers: {configured_providers}")

print("\n3. Testing individual accessor functions...")
for provider in providers:
    base_url = get_custom_base_url(provider)
    api_key = get_custom_api_key(provider)
    print(f"  - Provider '{provider}':")
    print(f"    - Base URL: {base_url}")
    print(f"    - API Key: {api_key[:4]}... (masked)")

print("\n4. Testing direct set/get functions...")
# Test setting another provider using individual functions
test_provider = "test_provider"
test_api_key = "test-api-key-9012"
test_base_url = "https://test-provider.example.com"

set_custom_api_key(test_provider, test_api_key)
set_custom_base_url(test_provider, test_base_url)

print(f"  - Set provider '{test_provider}' with direct functions")
print(f"  - Base URL: {get_custom_base_url(test_provider)}")
print(f"  - API Key: {get_custom_api_key(test_provider)}")

print("\n5. Verifying final provider list...")
final_providers = list_custom_providers()
print(f"  - Final configured providers: {final_providers}")

if len(final_providers) == 3 and all(p in final_providers for p in [*providers.keys(), test_provider]):
    print("\n✅ All tests passed! The circular import issue is fixed and custom provider functionality works correctly.")
else:
    print("\n❌ Test failed! Some providers weren't properly configured.") 