import os
import pandas as pd
import sys
import argparse

# Create test data
marker_data = pd.DataFrame({
    'gene': ['CD3D', 'CD3E', 'CD4', 'IL7R', 'CD8A', 'GZMB', 'CD19', 'MS4A1'],
    'p_val': [0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001],
    'avg_log2FC': [3.5, 3.2, 2.5, 2.0, -5.0, -3.0, -4.0, -4.5],
    'pct.1': [0.92, 0.90, 0.80, 0.75, 0.05, 0.03, 0.02, 0.01],
    'pct.2': [0.05, 0.06, 0.08, 0.10, 0.78, 0.65, 0.88, 0.82]
}).set_index('gene')

major_cluster_info = "Human PBMC"
comma_separated_genes = "CD3D, CD3E, CD4, IL7R"
annotation_history = "Previous analysis indicates T cell markers."

def check_env_var(var_name):
    """Check if environment variable exists and is not empty."""
    return var_name in os.environ and os.environ[var_name]

def test_openai():
    """Test OpenAI API call."""
    print("\n--- Testing OpenAI API ---")
    if not check_env_var("OPENAI_API_KEY"):
        print("OPENAI_API_KEY not set, skipping test")
        return False
    
    try:
        from openai import OpenAI
        
        # Manual API call to better control parameters
        client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))
        
        prompt = f"""
        You are an expert in single-cell annotation analysis. Your task is to evaluate and try to help finalize the single-cell annotation results, and generate next step for the excecuter to check. You can ask the excecuter to check certain group of genes expression, you can check for positive marker or negative marker. Provide your detailed reasoning. Note that you can also mention other possible cell types that are missed by the annotation. Note that mixed celltype is possible.

        context: the analylized cluster is from {major_cluster_info}, and has the following highly expressed markers:
        {comma_separated_genes}

        Below is the annotation analysis history:
        {annotation_history}
        """
        
        print("Sending request to OpenAI API...")
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[{"role": "user", "content": prompt}],
            max_tokens=3000
        )
        
        result = response.choices[0].message.content
        
        print("\nOpenAI API Result:")
        print(result)
        
        print("\nOpenAI test completed successfully!")
        return True
        
    except Exception as e:
        print(f"Error during OpenAI test: {str(e)}")
        return False

def test_anthropic():
    """Test Anthropic API call."""
    print("\n--- Testing Anthropic API ---")
    if not check_env_var("ANTHROPIC_API_KEY"):
        print("ANTHROPIC_API_KEY not set, skipping test")
        return False
    
    try:
        from anthropic import Anthropic
        
        # Manual API call to better control parameters
        client = Anthropic(api_key=os.environ.get("ANTHROPIC_API_KEY"))
        
        prompt = f"""
        You are an expert in single-cell annotation analysis. Your task is to evaluate and try to help finalize the single-cell annotation results, and generate next step for the excecuter to check. You can ask the excecuter to check certain group of genes expression, you can check for positive marker or negative marker. Provide your detailed reasoning. Note that you can also mention other possible cell types that are missed by the annotation. Note that mixed celltype is possible.

        context: the analylized cluster is from {major_cluster_info}, and has the following highly expressed markers:
        {comma_separated_genes}

        Below is the annotation analysis history:
        {annotation_history}
        """
        
        print("Sending request to Anthropic API...")
        response = client.messages.create(
            model="claude-3-haiku-20240307",
            max_tokens=1000,
            messages=[{"role": "user", "content": prompt}]
        )
        
        result = response.content[0].text
        
        print("\nAnthropic API Result:")
        print(result)
        
        print("\nAnthropic test completed successfully!")
        return True
        
    except Exception as e:
        print(f"Error during Anthropic test: {str(e)}")
        return False

def test_openrouter():
    """Test OpenRouter API call."""
    print("\n--- Testing OpenRouter API ---")
    if not check_env_var("OPENROUTER_API_KEY"):
        print("OPENROUTER_API_KEY not set, skipping test")
        return False
    
    try:
        import requests
        
        prompt = f"""
        You are an expert in single-cell annotation analysis. Your task is to evaluate and try to help finalize the single-cell annotation results, and generate next step for the excecuter to check. You can ask the excecuter to check certain group of genes expression, you can check for positive marker or negative marker. Provide your detailed reasoning. Note that you can also mention other possible cell types that are missed by the annotation. Note that mixed celltype is possible.

        context: the analylized cluster is from {major_cluster_info}, and has the following highly expressed markers:
        {comma_separated_genes}

        Below is the annotation analysis history:
        {annotation_history}
        """
        
        print("Sending request to OpenRouter API...")
        response = requests.post(
            url="https://openrouter.ai/api/v1/chat/completions",
            headers={
                "Authorization": f"Bearer {os.environ.get('OPENROUTER_API_KEY')}",
                "Content-Type": "application/json"
            },
            json={
                "model": "openai/gpt-3.5-turbo",
                "messages": [{"role": "user", "content": prompt}],
                "max_tokens": 1000
            }
        )
        
        if response.status_code == 200:
            result = response.json()["choices"][0]["message"]["content"]
            
            print("\nOpenRouter API Result:")
            print(result)
            
            print("\nOpenRouter test completed successfully!")
            return True
        else:
            print(f"Error: {response.status_code} - {response.text}")
            return False
        
    except Exception as e:
        print(f"Error during OpenRouter test: {str(e)}")
        return False

def test_all_providers():
    """Test all available providers."""
    results = {}
    
    if check_env_var("OPENAI_API_KEY"):
        results["openai"] = test_openai()
    else:
        print("Skipping OpenAI test - API key not set")
    
    if check_env_var("ANTHROPIC_API_KEY"):
        results["anthropic"] = test_anthropic()
    else:
        print("Skipping Anthropic test - API key not set")
    
    if check_env_var("OPENROUTER_API_KEY"):
        results["openrouter"] = test_openrouter()
    else:
        print("Skipping OpenRouter test - API key not set")
    
    # Print summary
    print("\n--- Test Results Summary ---")
    for provider, success in results.items():
        print(f"{provider}: {'PASS' if success else 'FAIL'}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Test LLM providers for iterative marker analysis")
    parser.add_argument("provider", nargs="?", choices=["openai", "anthropic", "openrouter", "all"], default="all", 
                        help="Which provider to test (default: all)")
    
    args = parser.parse_args()
    
    if args.provider == "all":
        test_all_providers()
    elif args.provider == "openai":
        test_openai()
    elif args.provider == "anthropic":
        test_anthropic()
    elif args.provider == "openrouter":
        test_openrouter() 