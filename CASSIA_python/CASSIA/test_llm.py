#!/usr/bin/env python
"""
Test script for call_llm function from merging_annotation.py
"""

import os
import sys
from dotenv import load_dotenv
from merging_annotation import call_llm

# Load environment variables from .env file if it exists
load_dotenv()


def test_openai():
    """Test OpenAI API call."""
    try:
        print("\n=== Testing OpenAI ===")
        response = call_llm(
            prompt="say hi",
            provider="openai",
            model="gpt-3.5-turbo",
            temperature=0.7
        )
        print(f"Response from OpenAI:\n{response}")
        return True
    except Exception as e:
        print(f"Error testing OpenAI: {str(e)}")
        return False

def test_anthropic():
    """Test Anthropic API call."""
    try:
        print("\n=== Testing Anthropic ===")
        response = call_llm(
            prompt="say hi",
            provider="anthropic",
            model="claude-3-sonnet-20240229",  # Using an available Claude model
            temperature=0.7
        )
        print(f"Response from Anthropic:\n{response}")
        return True
    except Exception as e:
        print(f"Error testing Anthropic: {str(e)}")
        return False

def test_openrouter():
    """Test OpenRouter API call."""
    try:
        print("\n=== Testing OpenRouter ===")
        response = call_llm(
            prompt="say hi",
            provider="openrouter",
            model="anthropic/claude-3-haiku",  # Using Anthropic through OpenRouter
            temperature=0.7
        )
        print(f"Response from OpenRouter:\n{response}")
        return True
    except Exception as e:
        print(f"Error testing OpenRouter: {str(e)}")
        return False

def test_with_system_prompt():
    """Test using a system prompt."""
    try:
        print("\n=== Testing with System Prompt ===")
        response = call_llm(
            prompt="say hi",
            provider="openai",
            model="gpt-3.5-turbo",
            system_prompt="You are a poet that specializes in haiku about technology.",
            temperature=0.7
        )
        print(f"Response with System Prompt:\n{response}")
        return True
    except Exception as e:
        print(f"Error testing with system prompt: {str(e)}")
        return False

def main():
    """Run tests based on available API keys."""
    print("Running LLM API tests...")
    
    # Track success/failure
    results = []
    
    # Only run tests for providers with API keys set
    if os.environ.get("OPENAI_API_KEY"):
        results.append(("OpenAI", test_openai()))
        results.append(("System Prompt", test_with_system_prompt()))
    else:
        print("Skipping OpenAI test: No API key found in environment")
    
    if os.environ.get("ANTHROPIC_API_KEY"):
        results.append(("Anthropic", test_anthropic()))
    else:
        print("Skipping Anthropic test: No API key found in environment")
    
    if os.environ.get("OPENROUTER_API_KEY"):
        results.append(("OpenRouter", test_openrouter()))
    else:
        print("Skipping OpenRouter test: No API key found in environment")
    
    # Print summary
    print("\n=== Test Results ===")
    for test_name, success in results:
        print(f"{test_name}: {'✅ Success' if success else '❌ Failed'}")
    
    # Return exit code based on success
    successful_tests = sum(1 for _, success in results if success)
    total_tests = len(results)
    
    if total_tests == 0:
        print("\nNo tests were run. Make sure your API keys are set in environment variables.")
        return 1
    else:
        print(f"\n{successful_tests}/{total_tests} tests passed.")
        return 0 if successful_tests == total_tests else 1

if __name__ == "__main__":
    sys.exit(main()) 









