#!/usr/bin/env python3
"""
CASSIA Symphony Compare - Multi-Model Cell Type Comparison with AI Consensus Building

This module provides an advanced cell type comparison function that orchestrates multiple AI models
to analyze and compare cell types based on marker expression patterns. It features:

- Parallel multi-model analysis for diverse perspectives
- Automatic consensus detection and discussion rounds when models disagree
- Beautiful interactive HTML reports with score progression tracking
- Structured CSV output for downstream analysis
- Support for custom model configurations
"""

import pandas as pd
import os
import requests
import re
import json
from typing import Dict, List, Tuple, Optional, Union
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from collections import Counter


def symphonyCompare(
    tissue: str,
    celltypes: List[str],
    marker_set: str,
    species: str = "human",
    model_preset: str = "budget",
    custom_models: Optional[List[str]] = None,
    output_dir: Optional[str] = None,
    output_basename: Optional[str] = None,
    enable_discussion: bool = True,
    max_discussion_rounds: int = 3,
    consensus_threshold: float = 2/3,
    generate_report: bool = True,
    api_key: Optional[str] = None,
    verbose: bool = True
) -> Dict:
    """
    Symphony Compare - Orchestrate multiple AI models to compare cell types with consensus building.
    
    This function conducts a comprehensive cell type comparison using multiple AI models in parallel,
    automatically triggering discussion rounds when models disagree on the best matching cell type.
    Think of it as a virtual panel of expert biologists debating and reaching consensus.
    
    Args:
        tissue (str): The tissue type being analyzed (e.g., "blood", "brain", "liver")
        celltypes (List[str]): List of 2-4 cell types to compare
        marker_set (str): Comma-separated string of gene markers to analyze
        species (str): Species being analyzed (default: "human")
        model_preset (str): Preset model configuration (default: "budget"). Options:
            - "budget": Cost-effective models (DeepSeek, Grok 4 Fast, Kimi K2, Gemini Flash)
            - "premium": High-performance ensemble (Gemini 3 Pro, Claude Sonnet 4.5, GPT-5.1, Grok 4)
            - "custom": Use custom_models list
        custom_models (List[str]): Custom list of models to use (when model_preset="custom")
        output_dir (str): Directory to save results (default: current directory)
        output_basename (str): Base name for output files (auto-generated if None)
        enable_discussion (bool): Enable automatic discussion rounds when no consensus (default: True)
        max_discussion_rounds (int): Maximum discussion rounds to perform (default: 3)
        consensus_threshold (float): Fraction of models that must agree for consensus (default: 2/3)
        generate_report (bool): Generate interactive HTML report (default: True)
        api_key (str): OpenRouter API key (uses environment variable if None)
        verbose (bool): Print progress messages (default: True)
    
    Returns:
        Dict containing:
            - 'results': List of all model responses and scores
            - 'consensus': The consensus cell type (if reached)
            - 'confidence': Confidence level of the consensus
            - 'csv_file': Path to the generated CSV file
            - 'html_file': Path to the generated HTML report (if enabled)
            - 'summary': Summary statistics of the comparison
    
    Raises:
        ValueError: If API key not set or invalid parameters provided
        
    Example:
        >>> results = symphonyCompare(
        ...     tissue="peripheral blood",
        ...     celltypes=["T cell", "B cell", "NK cell", "Monocyte"],
        ...     marker_set="CD3, CD4, CD8, CD19, CD20, CD16, CD56, CD14",
        ...     model_preset="symphony",
        ...     enable_discussion=True
        ... )
        >>> print(f"Consensus: {results['consensus']} (confidence: {results['confidence']:.1%})")
    """
    
    # Import the core comparison functionality
    try:
        # Try CASSIA package import first
        from CASSIA.evaluation.cell_type_comparison import (
            extract_celltype_scores,
            extract_discussion,
            generate_comparison_html_report,
            _call_model
        )
    except ImportError:
        try:
            # Try relative import (for Python package usage)
            from .cell_type_comparison import (
                extract_celltype_scores,
                extract_discussion,
                generate_comparison_html_report,
                _call_model
            )
        except ImportError:
            # Fall back to absolute import (for R/reticulate usage)
            import sys
            current_dir = os.path.dirname(os.path.abspath(__file__))
            if current_dir not in sys.path:
                sys.path.insert(0, current_dir)

            from cell_type_comparison import (
                extract_celltype_scores,
                extract_discussion,
                generate_comparison_html_report,
                _call_model
            )
    
    # Get API key
    if api_key is None:
        api_key = os.environ.get('OPENROUTER_API_KEY')
    if not api_key:
        raise ValueError(
            "OPENROUTER_API_KEY not found.\n"
            "Symphony Compare uses models from multiple providers (OpenAI, Anthropic, Google, etc.),\n"
            "so an OpenRouter API key is required to access all models through a unified endpoint.\n"
            "Get your API key at: https://openrouter.ai/keys\n"
            "Then set it as an environment variable: export OPENROUTER_API_KEY='your-key'\n"
            "Or pass it directly: symphonyCompare(..., api_key='your-key')"
        )
    
    # Input validation
    if not celltypes or len(celltypes) < 2 or len(celltypes) > 4:
        raise ValueError("Please provide 2-4 cell types to compare")
    
    # Set up output directory and filenames
    if output_dir is None:
        output_dir = os.getcwd()
    os.makedirs(output_dir, exist_ok=True)
    
    if output_basename is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        celltype_str = '_vs_'.join(ct.replace(' ', '_').replace('/', '_') for ct in celltypes)
        output_basename = f"symphony_compare_{species}_{tissue.replace(' ', '_')}_{celltype_str}_{timestamp}"
    
    csv_file = os.path.join(output_dir, f"{output_basename}.csv")
    html_file = os.path.join(output_dir, f"{output_basename}_report.html") if generate_report else None
    
    # Load model presets and personas from JSON config file
    config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'model_config.json')

    # Default fallback configuration
    default_presets = {
        "premium": [
            "google/gemini-3-pro-preview",
            "anthropic/claude-sonnet-4.5",
            "openai/gpt-5.1",
            "x-ai/grok-4"
        ],
        "budget": [
            "deepseek/deepseek-v3.2",
            "x-ai/grok-4-fast",
            "moonshotai/kimi-k2-thinking",
            "google/gemini-2.5-flash"
        ]
    }

    default_personas = {
        "google/gemini-3-pro-preview": "Dr. Emmy Noether",
        "anthropic/claude-sonnet-4.5": "Dr. Claude Shannon",
        "openai/gpt-5.1": "Dr. Albert Einstein",
        "x-ai/grok-4": "Dr. Marie Curie",
        "deepseek/deepseek-v3.2": "Dr. Alan Turing",
        "x-ai/grok-4-fast": "Dr. Nikola Tesla",
        "moonshotai/kimi-k2-thinking": "Dr. Ada Lovelace",
        "google/gemini-2.5-flash": "Dr. Rosalind Franklin"
    }

    # Try to load from JSON config file
    try:
        with open(config_file, 'r') as f:
            config = json.load(f)
        model_presets = config.get('presets', default_presets)
        model_personas = config.get('personas', default_personas)
        if verbose:
            print(f"  Loaded model configuration from: {config_file}")
    except (FileNotFoundError, json.JSONDecodeError) as e:
        model_presets = default_presets
        model_personas = default_personas
        if verbose:
            print(f"  Using default model configuration (config file not found or invalid)")
    
    # Select models based on preset or custom list
    if model_preset == "custom" and custom_models:
        model_list = custom_models
    elif model_preset in model_presets:
        model_list = model_presets[model_preset]
    else:
        if verbose:
            print(f"Warning: Unknown preset '{model_preset}'. Using 'budget' preset.")
        model_list = model_presets["budget"]
    
    # Get persona names
    model_to_persona = {m: model_personas.get(m, f"Researcher_{m.split('/')[-1]}") for m in model_list}
    
    if verbose:
        print(f"\n🎼 CASSIA Symphony Compare - Orchestrating {len(model_list)} AI Models")
        print(f"{'='*60}")
        print(f"📍 Tissue: {species} {tissue}")
        print(f"🔬 Comparing: {', '.join(celltypes)}")
        print(f"🧬 Markers: {marker_set}")
        print(f"🤖 Models: {', '.join([model_to_persona[m].split()[-1] for m in model_list])}")
        if enable_discussion:
            print(f"💬 Discussion: Enabled (max {max_discussion_rounds} rounds)")
        if model_preset == "budget":
            print(f"💡 Tip: For better performance, use model_preset='premium'")
        print(f"{'='*60}\n")
    
    # Construct initial prompt
    celltypes_list_str = "\n".join([f"- {ct}" for ct in celltypes])
    initial_prompt = f"""You are a professional biologist. Your task is to analyze how well a given marker set matches a list of cell types from {species} {tissue}.

For EACH of the following cell types, you must provide your analysis in a specific structured format.
The cell types to analyze are:
{celltypes_list_str}

The required output format for EACH cell type is:
<celltype>cell type name</celltype>
<reasoning>
Your detailed reasoning for the match, considering each marker's relevance.
</reasoning>
<score>A score from 0-100 indicating the match quality.</score>

Please provide a complete block of <celltype>, <reasoning>, and <score> for every cell type listed above.

Ranked marker set: {marker_set}"""
    
    # Initialize results storage
    all_results = []
    current_results = []
    rounds_performed = 0
    consensus_reached = False
    final_consensus = None
    
    # --- Initial Analysis Round ---
    if verbose:
        print("🎵 Movement I: Initial Analysis (Parallel Processing)")
    
    with ThreadPoolExecutor(max_workers=len(model_list)) as executor:
        future_to_model = {
            executor.submit(_call_model, model, initial_prompt, tissue, species, celltypes, 'initial', api_key, is_discussion_round=False): model 
            for model in model_list
        }
        for future in as_completed(future_to_model):
            result = future.result()
            result['researcher'] = model_to_persona.get(result['model'], result['model'])
            current_results.append(result)
            all_results.append(result)
    
    # Check for initial consensus
    winners = []
    valid_results = [r for r in current_results if r['status'] == 'success' and r['extracted_scores']]
    
    for result in valid_results:
        scores = {}
        for celltype, data in result['extracted_scores'].items():
            try:
                scores[celltype] = float(data['score'])
            except (ValueError, TypeError):
                scores[celltype] = -1
        if scores:
            winner = max(scores, key=scores.get)
            winners.append(winner)
    
    # Calculate consensus
    if winners:
        from collections import Counter
        winner_counts = Counter(winners)
        most_common = winner_counts.most_common(1)[0]
        consensus_ratio = most_common[1] / len(valid_results)
        
        if consensus_ratio >= consensus_threshold:
            consensus_reached = True
            final_consensus = most_common[0]
            if verbose:
                print(f"\n✅ Consensus reached! {consensus_ratio:.0%} agree on: {final_consensus}")
    
    # --- Discussion Rounds ---
    if enable_discussion and not consensus_reached and len(model_list) > 1 and max_discussion_rounds > 0:
        if verbose:
            print(f"\n🎵 Movement II: Discussion & Debate")
        
        for round_num in range(max_discussion_rounds):
            if consensus_reached:
                break
                
            rounds_performed = round_num + 1
            if verbose:
                print(f"\n  📢 Discussion Round {rounds_performed}/{max_discussion_rounds}")
            
            # Prepare discussion prompt
            all_responses = ""
            for res in current_results:
                researcher = res.get('researcher', res['model'])
                all_responses += f"\n--- Analysis from {researcher} ---\n"
                all_responses += f"{res['response']}\n"
            
            discussion_prompt_template = """You are a professional biologist participating in a panel discussion.
You are {persona_name}. Your colleagues' analyses are provided below. Review their arguments critically.

First, provide a brief critique of each colleague's analysis in a <discussion> block.
Then, provide your refined analysis for each cell type.

Original request:
{original_prompt}

Colleague analyses:
{all_responses}

Start with <discussion>, then provide your analysis for each cell type using <celltype>, <reasoning>, <score> format."""
            
            # Run discussion round
            discussion_results = []
            with ThreadPoolExecutor(max_workers=len(model_list)) as executor:
                round_name = f'discussion_{rounds_performed}'
                future_to_model = {}
                
                for model in model_list:
                    persona_name = model_to_persona[model]
                    this_prompt = discussion_prompt_template.format(
                        persona_name=persona_name,
                        original_prompt=initial_prompt,
                        all_responses=all_responses
                    )
                    future = executor.submit(_call_model, model, this_prompt, tissue, species, celltypes, round_name, api_key, is_discussion_round=True)
                    future_to_model[future] = model
                
                for future in as_completed(future_to_model):
                    result = future.result()
                    result['researcher'] = model_to_persona.get(result['model'], result['model'])
                    discussion_results.append(result)
                    all_results.append(result)
            
            current_results = discussion_results
            
            # Check consensus again
            winners = []
            valid_results = [r for r in current_results if r['status'] == 'success' and r['extracted_scores']]
            
            for result in valid_results:
                scores = {}
                for celltype, data in result['extracted_scores'].items():
                    try:
                        scores[celltype] = float(data['score'])
                    except (ValueError, TypeError):
                        scores[celltype] = -1
                if scores:
                    winner = max(scores, key=scores.get)
                    winners.append(winner)
            
            if winners:
                winner_counts = Counter(winners)
                most_common = winner_counts.most_common(1)[0]
                consensus_ratio = most_common[1] / len(valid_results)
                
                if consensus_ratio >= consensus_threshold:
                    consensus_reached = True
                    final_consensus = most_common[0]
                    if verbose:
                        print(f"\n  ✅ Consensus reached! {consensus_ratio:.0%} agree on: {final_consensus}")
                elif verbose:
                    print(f"  ⚡ No consensus yet ({consensus_ratio:.0%} for {most_common[0]})")
    
    # --- Generate Summary Statistics ---
    summary = {
        'total_rounds': 1 + rounds_performed,
        'models_used': len(model_list),
        'consensus_reached': consensus_reached,
        'consensus_celltype': final_consensus,
        'consensus_confidence': 0.0
    }
    
    # Calculate final scores for each cell type
    celltype_final_scores = {}
    for celltype in celltypes:
        scores = []
        final_round = f'discussion_{rounds_performed}' if rounds_performed > 0 else 'initial'
        for result in [r for r in all_results if r.get('round') == final_round]:
            if result['status'] == 'success' and celltype in result.get('extracted_scores', {}):
                try:
                    score = float(result['extracted_scores'][celltype]['score'])
                    scores.append(score)
                except (ValueError, TypeError):
                    pass
        if scores:
            celltype_final_scores[celltype] = {
                'mean': sum(scores) / len(scores),
                'min': min(scores),
                'max': max(scores),
                'std': pd.Series(scores).std() if len(scores) > 1 else 0
            }
    
    summary['celltype_scores'] = celltype_final_scores
    
    if final_consensus and final_consensus in celltype_final_scores:
        summary['consensus_confidence'] = celltype_final_scores[final_consensus]['mean'] / 100.0
    
    # --- Save Results ---
    if verbose:
        print(f"\n🎵 Movement III: Synthesis & Documentation")
    
    # Create CSV
    csv_data = []
    for result in all_results:
        base_row = {
            'model': result['model'],
            'researcher': result.get('researcher', result['model']),
            'tissue': result['tissue'],
            'species': result['species'], 
            'round': result.get('round', 'initial'),
            'status': result['status']
        }
        
        for celltype in celltypes:
            if celltype in result.get('extracted_scores', {}):
                base_row[f'{celltype}_score'] = result['extracted_scores'][celltype]['score']
                base_row[f'{celltype}_reasoning'] = result['extracted_scores'][celltype]['reasoning']
            else:
                base_row[f'{celltype}_score'] = 'N/A'
                base_row[f'{celltype}_reasoning'] = 'N/A'
        
        base_row['discussion'] = result.get('discussion', 'N/A')
        csv_data.append(base_row)
    
    df = pd.DataFrame(csv_data)
    df.to_csv(csv_file, index=False)
    
    if verbose:
        print(f"  📊 Results saved to: {csv_file}")
    
    # Generate HTML report
    if generate_report:
        html_content = generate_comparison_html_report(all_results, html_file)
        if verbose:
            print(f"  🎨 Interactive report: {html_file}")
    
    # --- Final Summary ---
    if verbose:
        print(f"\n🎼 Symphony Complete!")
        print(f"{'='*60}")
        print(f"📈 Performance Summary:")
        print(f"  • Models: {summary['models_used']} experts consulted")
        print(f"  • Rounds: {summary['total_rounds']} total (1 initial + {rounds_performed} discussion)")
        print(f"  • Consensus: {'✅ Yes' if consensus_reached else '❌ No'}")
        if final_consensus:
            print(f"  • Winner: {final_consensus} (confidence: {summary['consensus_confidence']:.1%})")
        print(f"\n📊 Detailed scores are available in the generated reports.")
    
    return {
        'results': all_results,
        'consensus': final_consensus,
        'confidence': summary['consensus_confidence'],
        'csv_file': csv_file,
        'html_file': html_file,
        'summary': summary,
        'dataframe': df
    }