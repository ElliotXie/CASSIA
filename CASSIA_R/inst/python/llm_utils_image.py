#!/usr/bin/env python3
"""
LLM Image Analysis Utilities for CASSIA Clustering Agent

This module provides functions to analyze clustering images using LLMs,
specifically designed for the clustering agent functionality.

Author: CASSIA Team
"""

import os
import json
import base64
import requests
from typing import Dict, List, Optional, Union
import warnings
warnings.filterwarnings('ignore')


def encode_image_to_base64(image_path: str) -> str:
    """
    Encode an image file to base64 string.
    
    Parameters:
    -----------
    image_path : str
        Path to the image file
        
    Returns:
    --------
    str : Base64 encoded image string
    """
    with open(image_path, "rb") as image_file:
        return base64.b64encode(image_file.read()).decode('utf-8')


def create_image_url(image_path: str) -> str:
    """
    Create data URL for image.
    
    Parameters:
    -----------
    image_path : str
        Path to the image file
        
    Returns:
    --------
    str : Data URL for the image
    """
    base64_image = encode_image_to_base64(image_path)
    return f"data:image/png;base64,{base64_image}"


def call_llm_with_images(
    text_prompt: str,
    image_paths: List[str],
    model: str = "google/gemini-2.5-flash",
    provider: str = "openrouter",
    api_key: Optional[str] = None,
    temperature: float = 0.1
) -> Dict:
    """
    Call LLM with images for analysis.
    
    Parameters:
    -----------
    text_prompt : str
        Text prompt for the LLM
    image_paths : List[str]
        List of image paths to analyze
    model : str
        Model name
    provider : str
        API provider
    api_key : Optional[str]
        API key
    temperature : float
        Sampling temperature
        
    Returns:
    --------
    Dict containing LLM response
    """
    
    # Get API key
    if api_key is None:
        if provider == "openrouter":
            api_key = os.environ.get('OPENROUTER_API_KEY')
        elif provider == "openai":
            api_key = os.environ.get('OPENAI_API_KEY')
        elif provider == "anthropic":
            api_key = os.environ.get('ANTHROPIC_API_KEY')
    
    if not api_key:
        raise ValueError(f"API key not found for provider {provider}")
    
    # Prepare message content
    message_content = [{"type": "text", "text": text_prompt}]
    
    # Add images
    for image_path in image_paths:
        if os.path.exists(image_path):
            image_url = create_image_url(image_path)
            message_content.append({
                "type": "image_url",
                "image_url": {"url": image_url}
            })
    
    # API call based on provider
    if provider == "openrouter":
        return _call_openrouter(message_content, model, api_key, temperature)
    elif provider == "openai":
        return _call_openai(message_content, model, api_key, temperature)
    elif provider == "anthropic":
        return _call_anthropic(message_content, model, api_key, temperature)
    else:
        raise ValueError(f"Unsupported provider: {provider}")


def _call_openrouter(message_content: List[Dict], model: str, api_key: str, temperature: float) -> Dict:
    """Call OpenRouter API with images."""
    
    response = requests.post(
        url="https://openrouter.ai/api/v1/chat/completions",
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
            "HTTP-Referer": "https://github.com/ElliotXie/CASSIA",
            "X-Title": "CASSIA Clustering Agent",
        },
        json={
            "model": model,
            "messages": [{"role": "user", "content": message_content}],
            "temperature": temperature,
            "max_tokens": 2000
        }
    )
    
    if response.status_code != 200:
        raise Exception(f"OpenRouter API error: {response.status_code} - {response.text}")
    
    return response.json()


def _call_openai(message_content: List[Dict], model: str, api_key: str, temperature: float) -> Dict:
    """Call OpenAI API with images."""
    
    response = requests.post(
        url="https://api.openai.com/v1/chat/completions",
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json"
        },
        json={
            "model": model,
            "messages": [{"role": "user", "content": message_content}],
            "temperature": temperature,
            "max_tokens": 2000
        }
    )
    
    if response.status_code != 200:
        raise Exception(f"OpenAI API error: {response.status_code} - {response.text}")
    
    return response.json()


def _call_anthropic(message_content: List[Dict], model: str, api_key: str, temperature: float) -> Dict:
    """Call Anthropic API with images."""
    
    response = requests.post(
        url="https://api.anthropic.com/v1/messages",
        headers={
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json",
            "anthropic-version": "2023-06-01"
        },
        json={
            "model": model,
            "messages": [{"role": "user", "content": message_content}],
            "temperature": temperature,
            "max_tokens": 2000
        }
    )
    
    if response.status_code != 200:
        raise Exception(f"Anthropic API error: {response.status_code} - {response.text}")
    
    return response.json()


def analyze_clustering_images(
    image_paths: List[str],
    combined_image_path: Optional[str],
    results_data: Union[List[Dict], Dict],
    analysis_type: str,
    model: str = "google/gemini-2.5-flash",
    provider: str = "openrouter",
    api_key: Optional[str] = None,
    max_images_per_call: int = 6
) -> Dict:
    """
    Analyze clustering images using LLM.
    
    Parameters:
    -----------
    image_paths : List[str]
        List of individual image paths
    combined_image_path : Optional[str]
        Path to combined image
    results_data : Union[List[Dict], Dict]
        Results data associated with images
    analysis_type : str
        Type of analysis ("parameter_grid", "resolution_stability", "final_recommendation")
    model : str
        LLM model
    provider : str
        API provider
    api_key : Optional[str]
        API key
    max_images_per_call : int
        Maximum images per call
        
    Returns:
    --------
    Dict containing analysis results
    """
    
    print(f"Analyzing {len(image_paths)} images for {analysis_type}...")
    
    if analysis_type == "parameter_grid":
        return _analyze_parameter_grid(image_paths, combined_image_path, results_data, 
                                     model, provider, api_key, max_images_per_call)
    elif analysis_type == "resolution_stability":
        return _analyze_resolution_stability(image_paths, combined_image_path, results_data,
                                           model, provider, api_key, max_images_per_call)
    elif analysis_type == "final_recommendation":
        return _analyze_final_recommendation(image_paths, combined_image_path, results_data,
                                          model, provider, api_key, max_images_per_call)
    else:
        raise ValueError(f"Unknown analysis type: {analysis_type}")


def _analyze_parameter_grid(
    image_paths: List[str],
    combined_image_path: Optional[str],
    results_data: List[Dict],
    model: str,
    provider: str,
    api_key: Optional[str],
    max_images_per_call: int
) -> Dict:
    """Analyze parameter grid images."""
    
    # Create detailed prompt for parameter analysis
    prompt = """
You are an expert in single-cell RNA sequencing data analysis and clustering. 
I will show you UMAP plots generated with different clustering parameters.

Please analyze these clustering visualizations and provide:

1. **Visual Quality Assessment**: Rate each plot on cluster separation, compactness, and overall structure (1-10 scale)

2. **Parameter Recommendations**: Based on the plots, which parameter combination produces:
   - Well-separated clusters
   - Appropriate cluster sizes (not too many tiny clusters, not too few large ones)
   - Biologically meaningful structure
   - Good overall balance

3. **Top 3 Choices**: Rank your top 3 parameter combinations with reasoning

4. **Specific Observations**: Note any artifacts, over-clustering, or under-clustering

For each plot, the parameters are shown in the title:
- n_pcs: Number of principal components used
- min_dist: UMAP min_dist parameter (lower = tighter clusters)
- n_neighbors: Number of neighbors for graph construction

Please provide your analysis in JSON format:
{
    "analysis": "Your detailed analysis here",
    "top_3_recommendations": [
        {
            "rank": 1,
            "parameters": {"n_pcs": X, "min_dist": Y, "n_neighbors": Z},
            "score": X.X,
            "reasoning": "Why this is the best choice"
        },
        ... (for ranks 2 and 3)
    ],
    "best_parameters": {"n_pcs": X, "min_dist": Y, "n_neighbors": Z},
    "confidence": X.X
}
"""

    # Use combined image if available, otherwise analyze individual images
    images_to_analyze = [combined_image_path] if combined_image_path else image_paths[:max_images_per_call]
    
    try:
        response = call_llm_with_images(prompt, images_to_analyze, model, provider, api_key)
        
        # Extract content based on provider
        if provider == "openrouter" or provider == "openai":
            content = response['choices'][0]['message']['content']
        elif provider == "anthropic":
            content = response['content'][0]['text']
        else:
            content = str(response)
        
        # Try to parse JSON response
        try:
            import re
            json_match = re.search(r'\{.*\}', content, re.DOTALL)
            if json_match:
                analysis_result = json.loads(json_match.group())
            else:
                analysis_result = {"analysis": content, "best_parameters": {}}
        except (json.JSONDecodeError, AttributeError):
            analysis_result = {"analysis": content, "best_parameters": {}}
        
        return analysis_result
        
    except Exception as e:
        print(f"Error in LLM analysis: {e}")
        return _fallback_parameter_analysis(results_data)


def _analyze_resolution_stability(
    image_paths: List[str],
    combined_image_path: Optional[str],
    results_data: List[Dict],
    model: str,
    provider: str,
    api_key: Optional[str],
    max_images_per_call: int
) -> Dict:
    """Analyze resolution stability images."""
    
    # Create stability analysis prompt
    prompt = """
You are analyzing clustering stability across different resolutions in single-cell data.

I will show you UMAP plots at different clustering resolutions, along with stability metrics.

Please analyze:

1. **Stability vs Resolution Trade-off**: How does clustering stability change with resolution?

2. **Optimal Resolution**: Which resolution provides the best balance of:
   - High stability (consistent clustering across iterations)
   - Meaningful cluster granularity
   - Biological interpretability

3. **Visual Assessment**: Look for:
   - Over-clustering (too many tiny clusters)
   - Under-clustering (large heterogeneous clusters)
   - Stability artifacts

Each plot shows:
- Resolution value
- Number of clusters (±standard deviation)
- Stability score (ARI) (±standard deviation)

Higher stability scores (closer to 1.0) indicate more consistent clustering.

Please provide your analysis in JSON format:
{
    "stability_analysis": "Your detailed analysis",
    "recommended_resolution": X.X,
    "reasoning": "Why this resolution is optimal",
    "confidence": X.X,
    "trade_off_assessment": "Analysis of stability vs granularity trade-off"
}
"""

    images_to_analyze = [combined_image_path] if combined_image_path else image_paths[:max_images_per_call]
    
    try:
        response = call_llm_with_images(prompt, images_to_analyze, model, provider, api_key)
        
        # Extract content
        if provider == "openrouter" or provider == "openai":
            content = response['choices'][0]['message']['content']
        elif provider == "anthropic":
            content = response['content'][0]['text']
        else:
            content = str(response)
        
        # Parse JSON
        try:
            import re
            json_match = re.search(r'\{.*\}', content, re.DOTALL)
            if json_match:
                analysis_result = json.loads(json_match.group())
            else:
                analysis_result = {"stability_analysis": content, "recommended_resolution": 0.7}
        except (json.JSONDecodeError, AttributeError):
            analysis_result = {"stability_analysis": content, "recommended_resolution": 0.7}
        
        return analysis_result
        
    except Exception as e:
        print(f"Error in stability analysis: {e}")
        return _fallback_stability_analysis(results_data)


def _analyze_final_recommendation(
    image_paths: List[str],
    combined_image_path: Optional[str],
    results_data: Dict,
    model: str,
    provider: str,
    api_key: Optional[str],
    max_images_per_call: int
) -> Dict:
    """Analyze final recommendation combining parameter grid and stability results."""
    
    prompt = """
You are making final recommendations for single-cell clustering parameters.

I will show you:
1. Parameter grid results (different UMAP parameter combinations)
2. Resolution stability analysis

Please provide final recommendations that combine:
- Best visual clustering from parameter grid
- Optimal resolution from stability analysis

Your final recommendation should balance:
- Cluster separation and compactness
- Appropriate granularity
- Stability across iterations
- Biological interpretability

Please provide your analysis in JSON format:
{
    "final_analysis": "Your comprehensive analysis",
    "recommended_parameters": {
        "n_pcs": X,
        "min_dist": X.X,
        "n_neighbors": X,
        "resolution": X.X
    },
    "reasoning": "Why these parameters are optimal",
    "confidence": X.X,
    "expected_outcomes": "What to expect from these parameters"
}
"""

    try:
        response = call_llm_with_images(prompt, image_paths, model, provider, api_key)
        
        # Extract content
        if provider == "openrouter" or provider == "openai":
            content = response['choices'][0]['message']['content']
        elif provider == "anthropic":
            content = response['content'][0]['text']
        else:
            content = str(response)
        
        # Parse JSON
        try:
            import re
            json_match = re.search(r'\{.*\}', content, re.DOTALL)
            if json_match:
                analysis_result = json.loads(json_match.group())
            else:
                analysis_result = {
                    "final_analysis": content,
                    "recommended_parameters": {"n_pcs": 30, "min_dist": 0.2, "n_neighbors": 30, "resolution": 0.7}
                }
        except (json.JSONDecodeError, AttributeError):
            analysis_result = {
                "final_analysis": content,
                "recommended_parameters": {"n_pcs": 30, "min_dist": 0.2, "n_neighbors": 30, "resolution": 0.7}
            }
        
        return analysis_result
        
    except Exception as e:
        print(f"Error in final analysis: {e}")
        return {
            "final_analysis": f"Error occurred: {e}",
            "recommended_parameters": {"n_pcs": 30, "min_dist": 0.2, "n_neighbors": 30, "resolution": 0.7}
        }


def _fallback_parameter_analysis(results_data: List[Dict]) -> Dict:
    """Fallback analysis for parameter grid."""
    
    # Find moderate cluster count
    cluster_counts = [r['n_clusters'] for r in results_data]
    median_clusters = np.median(cluster_counts)
    best_idx = np.argmin([abs(count - median_clusters) for count in cluster_counts])
    best_result = results_data[best_idx]
    
    return {
        "analysis": "Fallback analysis: selected parameters producing moderate cluster count",
        "best_parameters": {
            "n_pcs": best_result['n_pcs'],
            "min_dist": best_result['min_dist'],
            "n_neighbors": best_result['n_neighbors']
        },
        "confidence": 0.5
    }


def _fallback_stability_analysis(results_data: List[Dict]) -> Dict:
    """Fallback analysis for resolution stability."""
    
    # Find most stable resolution
    import numpy as np
    stabilities = [r['mean_stability'] for r in results_data]
    best_idx = np.argmax(stabilities)
    best_resolution = results_data[best_idx]['resolution']
    
    return {
        "stability_analysis": "Fallback analysis: selected most stable resolution",
        "recommended_resolution": best_resolution,
        "confidence": 0.5
    }


if __name__ == "__main__":
    print("LLM Image Analysis Utilities for CASSIA Clustering Agent")
    print("This module provides image analysis functions for clustering optimization.")