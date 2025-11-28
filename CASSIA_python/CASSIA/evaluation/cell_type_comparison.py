import pandas as pd
import os
import requests
import re
import json
from typing import Dict, List, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

# Import CASSIA logger for actionable error messages
try:
    from CASSIA.core.logging_config import get_logger
except ImportError:
    try:
        from .logging_config import get_logger
    except ImportError:
        from logging_config import get_logger

logger = get_logger(__name__)


def extract_celltype_scores(response_text: str, celltypes: List[str]) -> Dict[str, Dict[str, str]]:
    """
    Extract scores and reasoning for each celltype from the LLM response.
    
    Args:
        response_text (str): The raw response text from the LLM
        celltypes (list): List of cell types being compared
        
    Returns:
        dict: Dictionary with celltype as key and dict of score/reasoning as value
    """
    results = {}
    
    # Try to find celltype-specific responses
    for celltype in celltypes:
        # Look for celltype tags
        celltype_pattern = rf'<celltype>{re.escape(celltype)}</celltype>(.*?)(?=<celltype>|$)'
        celltype_match = re.search(celltype_pattern, response_text, re.DOTALL | re.IGNORECASE)
        
        if celltype_match:
            celltype_content = celltype_match.group(1)
        else:
            # Fallback: look for mentions of the celltype
            celltype_content = response_text
        
        # Extract reasoning
        reasoning_pattern = r'<reasoning>(.*?)</reasoning>'
        reasoning_match = re.search(reasoning_pattern, celltype_content, re.DOTALL | re.IGNORECASE)
        reasoning = reasoning_match.group(1).strip() if reasoning_match else "No reasoning found"
        
        # Extract score
        score_pattern = r'<score>(\d+(?:\.\d+)?)</score>'
        score_match = re.search(score_pattern, celltype_content, re.IGNORECASE)
        score = score_match.group(1) if score_match else "No score found"
        
        results[celltype] = {
            'score': score,
            'reasoning': reasoning
        }
    
    # If no structured responses found, try to extract any scores mentioned
    if not any(results[ct]['score'] != "No score found" for ct in celltypes):
        # Fallback: look for any numbers that might be scores
        score_patterns = [
            rf'({re.escape(ct)}[^0-9]*(\d+(?:\.\d+)?))'
            for ct in celltypes
        ]
        
        for i, pattern in enumerate(score_patterns):
            matches = re.findall(pattern, response_text, re.IGNORECASE)
            if matches:
                results[celltypes[i]]['score'] = matches[0][1]
                results[celltypes[i]]['reasoning'] = f"Extracted from: {matches[0][0]}"
    
    return results


def extract_discussion(response_text: str) -> str:
    """Extracts the content of the <discussion> tag from the response."""
    discussion_pattern = r'<discussion>(.*?)</discussion>'
    discussion_match = re.search(discussion_pattern, response_text, re.DOTALL | re.IGNORECASE)
    if discussion_match:
        return discussion_match.group(1).strip()
    return "No discussion found"


def generate_comparison_html_report(all_results: List[Dict], output_file: str = None) -> str:
    """
    Generate an HTML report for cell type comparison results, including discussion progression.
    """
    if not all_results:
        return "<html><body><h1>No results to display</h1></body></html>"
    
    # Get unique cell types, models, and researchers
    celltypes = set()
    models = set()
    researchers = set()
    rounds = set()
    for result in all_results:
        if 'extracted_scores' in result:
            celltypes.update(result['extracted_scores'].keys())
        models.add(result.get('model', 'Unknown'))
        researchers.add(result.get('researcher', result.get('model', 'Unknown')))
        rounds.add(result.get('round', 'initial'))
    celltypes = sorted(list(celltypes))
    models = sorted(list(models))
    researchers = sorted(list(researchers))
    rounds = sorted(list(rounds), key=lambda r: (0 if r == 'initial' else int(r.split('_')[-1]) if r.startswith('discussion_') else 999))

    # Group results by round
    round_to_results = {r: [] for r in rounds}
    for result in all_results:
        round_to_results[result.get('round', 'initial')].append(result)

    # Start building HTML
    html_content = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Cell Type Comparison Report</title>
        <style>
            body {
                font-family: 'Segoe UI', Roboto, -apple-system, sans-serif;
                max-width: 1400px;
                margin: 0 auto;
                padding: 20px;
                background-color: #f5f7fa;
                line-height: 1.6;
            }
            .container {
                background-color: white;
                padding: 30px;
                border-radius: 12px;
                box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            }
            .header {
                text-align: center;
                margin-bottom: 40px;
                padding-bottom: 20px;
                border-bottom: 3px solid #e74c3c;
            }
            .title {
                font-size: 2.5rem;
                font-weight: bold;
                color: #2c3e50;
                margin: 0;
                background: linear-gradient(135deg, #e74c3c, #c0392b);
                -webkit-background-clip: text;
                -webkit-text-fill-color: transparent;
            }
            .subtitle {
                font-size: 1.2rem;
                color: #7f8c8d;
                margin-top: 10px;
            }
            .round-section {
                margin-bottom: 40px;
                border-radius: 10px;
                background: #f0f4f8;
                box-shadow: 0 2px 8px rgba(0,0,0,0.06);
                padding: 20px 30px 30px 30px;
            }
            .round-header {
                font-size: 1.5rem;
                font-weight: bold;
                color: #c0392b;
                margin-bottom: 18px;
                display: flex;
                align-items: center;
            }
            .round-badge {
                background: #e74c3c;
                color: white;
                border-radius: 20px;
                padding: 6px 18px;
                font-size: 1.1rem;
                margin-right: 16px;
                font-weight: bold;
            }
            .researcher-block {
                background: #fff;
                border-radius: 8px;
                box-shadow: 0 1px 4px rgba(44,62,80,0.07);
                margin-bottom: 18px;
                padding: 18px 22px 18px 22px;
                border-left: 6px solid #e74c3c;
                position: relative;
            }
            .researcher-header {
                display: flex;
                align-items: center;
                margin-bottom: 8px;
            }
            .researcher-avatar {
                width: 38px;
                height: 38px;
                border-radius: 50%;
                background: linear-gradient(135deg, #e74c3c, #c0392b);
                color: white;
                display: flex;
                align-items: center;
                justify-content: center;
                font-size: 1.3rem;
                font-weight: bold;
                margin-right: 14px;
                box-shadow: 0 2px 6px rgba(44,62,80,0.10);
            }
            .researcher-name {
                font-size: 1.2rem;
                font-weight: bold;
                color: #2c3e50;
            }
            .model-name {
                font-size: 0.95rem;
                color: #888;
                margin-left: 10px;
            }
            .score-table {
                width: 100%;
                border-collapse: collapse;
                margin: 20px 0;
                box-shadow: 0 2px 8px rgba(0,0,0,0.1);
                border-radius: 8px;
                overflow: hidden;
            }
            .score-table th {
                background: linear-gradient(135deg, #34495e, #2c3e50);
                color: white;
                padding: 15px 12px;
                text-align: center;
                font-weight: bold;
                font-size: 1.1rem;
            }
            .score-table td {
                padding: 12px;
                text-align: center;
                border-bottom: 1px solid #ddd;
            }
            .score-table tr:nth-child(even) {
                background-color: #f8f9fa;
            }
            .score-table tr:hover {
                background-color: #e9ecef;
            }
            .majority-winner {
                background-color: #d4edda !important;
                border-left: 4px solid #28a745;
            }
            .majority-winner:hover {
                background-color: #c3e6cb !important;
            }
            .cell-type-name {
                font-weight: bold;
                color: #2c3e50;
                text-align: left !important;
                padding-left: 20px !important;
            }
            .high-score {
                background-color: #d4edda !important;
                color: #155724;
                font-weight: bold;
                border-radius: 4px;
            }
            .medium-score {
                background-color: #fff3cd !important;
                color: #856404;
                font-weight: bold;
                border-radius: 4px;
            }
            .low-score {
                background-color: #f8d7da !important;
                color: #721c24;
                font-weight: bold;
                border-radius: 4px;
            }
            .celltype-row {
                padding: 20px;
                border-bottom: 1px solid #eee;
                background-color: white;
            }
            .celltype-row:last-child {
                border-bottom: none;
            }
            .celltype-row:nth-child(even) {
                background-color: #f8f9fa;
            }
            .celltype-header {
                display: flex;
                align-items: center;
                justify-content: space-between;
                margin-bottom: 15px;
                padding-bottom: 10px;
                border-bottom: 2px solid #e74c3c;
            }
            .celltype-name {
                font-size: 1.3rem;
                font-weight: bold;
                color: #2c3e50;
            }
            .score-badge {
                background: linear-gradient(135deg, #e74c3c, #c0392b);
                color: white;
                padding: 8px 16px;
                border-radius: 25px;
                font-size: 1rem;
                font-weight: bold;
            }
            .reasoning-preview {
                margin-bottom: 10px;
                color: #666;
                font-style: italic;
                font-size: 0.95rem;
            }
            .reasoning-section {
                margin-bottom: 15px;
            }
            .reasoning-text {
                background-color: #f8f9fa;
                padding: 15px;
                border-radius: 6px;
                border-left: 4px solid #e74c3c;
                color: #555;
                text-align: justify;
                line-height: 1.6;
                margin-top: 10px;
            }
            .discussion-section {
                padding: 15px 20px;
                background-color: #f8f9fa;
                border-bottom: 1px solid #ddd;
                margin-bottom: 10px;
            }
            .discussion-header {
                font-size: 1.1rem;
                font-weight: bold;
                color: #2980b9;
                margin-bottom: 10px;
            }
            .discussion-text {
                background-color: #fff;
                padding: 12px;
                border-radius: 6px;
                border-left: 4px solid #3498db;
                font-style: italic;
                color: #555;
            }
            .toggle-button {
                background-color: #6c757d;
                color: white;
                border: none;
                padding: 6px 12px;
                border-radius: 4px;
                cursor: pointer;
                font-size: 0.85rem;
                transition: background-color 0.3s;
            }
            .toggle-button:hover {
                background-color: #5a6268;
            }
            .details-title {
                font-size: 1.8rem;
                font-weight: bold;
                color: #2c3e50;
                margin: 40px 0 20px 0;
                text-align: center;
                border-bottom: 2px solid #e74c3c;
                padding-bottom: 10px;
            }
        </style>
        <script>
            function toggleReasoning(buttonId, contentId) {
                const content = document.getElementById(contentId);
                const button = document.getElementById(buttonId);
                if (content.style.display === 'none' || content.style.display === '') {
                    content.style.display = 'block';
                    button.textContent = 'Hide Reasoning';
                } else {
                    content.style.display = 'none';
                    button.textContent = 'Show Full Reasoning';
                }
            }
            function toggleSummaryTable(roundId) {
                const table = document.getElementById('summary_table_' + roundId);
                const btn = document.getElementById('toggle_summary_btn_' + roundId);
                if (table.style.display === 'none' || table.style.display === '') {
                    table.style.display = 'block';
                    btn.textContent = 'Hide Summary Table';
                } else {
                    table.style.display = 'none';
                    btn.textContent = 'Show Summary Table';
                }
            }
        </script>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1 class="title">Cell Type Comparison Report</h1>
                <p class="subtitle">Multi-Model Analysis Results</p>
            </div>
    """
    
    # --- Multi-Round Summary Tables ---
    html_content += '<div class="details-title">📊 Summary Score Tables by Round</div>'
    for idx, round_name in enumerate(rounds):
        pretty_round = "Initial Analysis" if round_name == 'initial' else f"Discussion Round {round_name.split('_')[-1]}"
        is_final = (idx == len(rounds) - 1)
        # Toggle button for previous rounds
        if not is_final:
            html_content += f'<button class="toggle-button" id="toggle_summary_btn_{round_name}" onclick="toggleSummaryTable(\'{round_name}\')">Show Summary Table</button>'
        html_content += f'<div class="summary-section" id="summary_table_{round_name}" style="display: {"block" if is_final else "none"};">'
        html_content += f'<div class="summary-title">{pretty_round} Summary Table</div>'
        html_content += """
                <table class="score-table">
                    <thead>
                        <tr>
                            <th style="text-align: left;">Cell Type</th>
        """
        # Get models for this round
        round_models = [r.get('model', 'Unknown') for r in round_to_results[round_name]]
        round_models = sorted(list(set(round_models)))
        for model in round_models:
            model_short = model.split('/')[-1] if '/' in model else model
            html_content += f"<th>{model_short}</th>"
        html_content += """
                            <th>Average</th>
                            <th>Majority Votes</th>
                        </tr>
                    </thead>
                    <tbody>
        """
        # Calculate averages and majority voting for this round
        celltype_scores = {}
        model_winners = []
        for celltype in celltypes:
            scores = []
            for model in round_models:
                for result in round_to_results[round_name]:
                    if result.get('model') == model and 'extracted_scores' in result:
                        if celltype in result['extracted_scores']:
                            score_str = result['extracted_scores'][celltype]['score']
                            try:
                                score_num = float(score_str)
                                scores.append(score_num)
                            except (ValueError, TypeError):
                                logger.warning(
                                    f"Could not convert score '{score_str}' to number for cell type '{celltype}'. "
                                    f"This score will be excluded from analysis. "
                                    f"Check if the LLM response format is correct."
                                )
                            break
            celltype_scores[celltype] = scores
        for model in round_models:
            model_scores = {}
            for result in round_to_results[round_name]:
                if result.get('model') == model and 'extracted_scores' in result:
                    for celltype in celltypes:
                        if celltype in result['extracted_scores']:
                            score_str = result['extracted_scores'][celltype]['score']
                            try:
                                score_num = float(score_str)
                                model_scores[celltype] = score_num
                            except (ValueError, TypeError):
                                logger.warning(
                                    f"Could not convert score '{score_str}' for '{celltype}' from model '{model}'. "
                                    f"Using 0 as fallback score. "
                                    f"Check if the LLM response format is correct."
                                )
                                model_scores[celltype] = 0
                    break
            if model_scores:
                winner = max(model_scores, key=model_scores.get)
                model_winners.append(winner)
        majority_votes = {}
        for celltype in celltypes:
            majority_votes[celltype] = model_winners.count(celltype)
        for celltype in celltypes:
            votes = majority_votes.get(celltype, 0)
            row_class = "majority-winner" if votes >= len(round_models) // 2 + 1 else ""
            html_content += f'<tr class="{row_class}"><td class="cell-type-name">{celltype}</td>'
            for model in round_models:
                score = "N/A"
                for result in round_to_results[round_name]:
                    if result.get('model') == model and 'extracted_scores' in result:
                        if celltype in result['extracted_scores']:
                            score = result['extracted_scores'][celltype]['score']
                            break
                html_content += f'<td>{score}</td>'
            scores = celltype_scores.get(celltype, [])
            if scores:
                avg_score = sum(scores) / len(scores)
                avg_class = ""
                if avg_score >= 70:
                    avg_class = "high-score"
                elif avg_score >= 40:
                    avg_class = "medium-score"
                else:
                    avg_class = "low-score"
                html_content += f'<td class="{avg_class}"><strong>{avg_score:.1f}</strong></td>'
            else:
                html_content += '<td>N/A</td>'
            vote_class = ""
            if votes >= len(round_models) // 2 + 1:
                vote_class = "high-score"
            elif votes > 0:
                vote_class = "medium-score"
            else:
                vote_class = "low-score"
            html_content += f'<td class="{vote_class}"><strong>{votes}/{len(round_models)}</strong></td>'
            html_content += "</tr>"
        html_content += """
                    </tbody>
                </table>
            </div>
        """

    # --- Discussion Progression Section ---
    html_content += '<div class="details-title">🗣️ Discussion Progression by Round</div>'
    for round_name in rounds:
        pretty_round = "Initial Analysis" if round_name == 'initial' else f"Discussion Round {round_name.split('_')[-1]}"
        html_content += f'<div class="round-section">'
        html_content += f'<div class="round-header"><span class="round-badge">{pretty_round}</span></div>'
        for result in round_to_results[round_name]:
            researcher = result.get('researcher', result.get('model', 'Unknown'))
            model = result.get('model', 'Unknown')
            extracted_scores = result.get('extracted_scores', {})
            discussion = result.get('discussion')
            avatar_letter = researcher[0] if researcher else '?'  # Use first letter of name
            html_content += f'<div class="researcher-block">'
            html_content += f'<div class="researcher-header">'
            html_content += f'<div class="researcher-avatar">{avatar_letter}</div>'
            html_content += f'<span class="researcher-name">{researcher}</span>'
            html_content += f'<span class="model-name">({model})</span>'
            html_content += '</div>'
            # Show discussion block if present (for discussion rounds)
            if discussion and discussion != "No discussion found":
                discussion_html = discussion.replace('\n', '<br>')
                html_content += f'<div class="discussion-section">'
                html_content += f'<div class="discussion-header">Peer Review & Critique</div>'
                html_content += f'<div class="discussion-text">{discussion_html}</div>'
                html_content += '</div>'
            # Show scores and reasoning for each cell type
            for celltype in celltypes:
                if celltype in extracted_scores:
                    score = extracted_scores[celltype]['score']
                    reasoning = extracted_scores[celltype]['reasoning']
                    words = reasoning.split()
                    preview = ' '.join(words[:12]) + '...' if len(words) > 12 else reasoning
                else:
                    score = "N/A"
                    reasoning = "No data available for this cell type"
                    preview = reasoning
                html_content += f"""
                <div class="celltype-row">
                    <div class="celltype-header">
                        <span class="celltype-name">{celltype}</span>
                        <span class="score-badge">Score: {score}/100</span>
                    </div>
                    <div class="reasoning-preview">💭 {preview}</div>
                    <button class="toggle-button" id="btn_{round_name}_{researcher}_{celltype.replace(' ', '_')}" 
                            onclick="toggleReasoning('btn_{round_name}_{researcher}_{celltype.replace(' ', '_')}', 'reasoning_{round_name}_{researcher}_{celltype.replace(' ', '_')}')">
                        Show Full Reasoning
                    </button>
                    <div class="reasoning-text" id="reasoning_{round_name}_{researcher}_{celltype.replace(' ', '_')}" style="display: none;">
                        {reasoning}
                    </div>
                </div>
                """
            html_content += '</div>'  # researcher-block
        html_content += '</div>'  # round-section

    html_content += """
        </div>
    </body>
    </html>
    """
    
    # Save to file if specified
    if output_file:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"HTML report saved to {output_file}")
    
    return html_content


def _call_model(model: str, prompt: str, tissue: str, species: str, celltypes: List[str], round_name: str, api_key: str, is_discussion_round: bool = False) -> Dict:
    """Helper function to make a single API call to a model."""
    try:
        response = requests.post(
            url="https://openrouter.ai/api/v1/chat/completions",
            headers={
                "Authorization": f"Bearer {api_key}",
                "HTTP-Referer": "https://elliotxie.github.io/CASSIA/",
                "X-Title": "CASSIA",
                "Content-Type": "application/json"
            },
            json={
                "model": model,
                "messages": [{"role": "user", "content": prompt}]
            }
        )

        if response.status_code == 200:
            response_data = response.json()
            model_response = response_data['choices'][0]['message']['content']
            extracted_scores = extract_celltype_scores(model_response, celltypes)
            
            discussion = None
            if is_discussion_round:
                discussion = extract_discussion(model_response)

            print(f"Model ({round_name}): {model} [OK]")

            result_dict = {
                'model': model, 'tissue': tissue, 'species': species,
                'cell_types': ', '.join(celltypes), 'response': model_response,
                'extracted_scores': extracted_scores, 'status': 'success', 'round': round_name
            }
            if discussion:
                result_dict['discussion'] = discussion
            return result_dict
        else:
            print(f"Model ({round_name}): {model} [ERROR {response.status_code}]")
            return {
                'model': model, 'tissue': tissue, 'species': species,
                'cell_types': ', '.join(celltypes), 'response': f"Error: {response.status_code}",
                'extracted_scores': {}, 'status': 'error', 'round': round_name
            }
    except Exception as e:
        print(f"Model ({round_name}): {model} [EXCEPTION]")
        return {
            'model': model, 'tissue': tissue, 'species': species,
            'cell_types': ', '.join(celltypes), 'response': f"Exception: {str(e)}",
            'extracted_scores': {}, 'status': 'error', 'round': round_name
        }


def compareCelltypes(tissue, celltypes, marker_set, species="human", model_preset="default", model_list=None, output_file=None, generate_html_report=True, discussion_mode=False, discussion_rounds=1):
    """
    Compare cell types using multiple AI models and generate comparison scores with structured output.
    
    Args:
        tissue (str): Tissue type being analyzed
        celltypes (list): List of 2-4 cell types to compare
        marker_set (str): Ranked marker set for comparison
        species (str): Species being analyzed (default: "human")
        model_preset (str): a preset of models to use for comparison. Options: "default", "budget". Overridden by model_list.
        model_list (list): List of models to use for comparison (optional)
        output_file (str): Output CSV filename (optional)
        generate_html_report (bool): Whether to generate HTML report (default: True)
        discussion_mode (bool): Whether to enable discussion mode if no consensus (default: False)
        discussion_rounds (int): how many rounds of discussion to perform (default: 1)
    
    Returns:
        dict: Dictionary containing results and HTML content if generated
    
    Raises:
        ValueError: If API key not set or invalid number of cell types provided
    """
    # Get API key from environment variable
    OPENROUTER_API_KEY = os.environ.get('OPENROUTER_API_KEY')
    if not OPENROUTER_API_KEY:
        raise ValueError("OPENROUTER_API_KEY environment variable is not set")
    
    # Input validation
    if not celltypes or len(celltypes) < 2 or len(celltypes) > 4:
        raise ValueError("Please provide 2-4 cell types to compare")
    
    # Generate default output filename based on celltypes if none provided
    if output_file is None:
        # Create a sanitized version of cell types for the filename
        celltype_str = '_vs_'.join(ct.replace(' ', '_') for ct in celltypes)
        output_file = f"model_comparison_{celltype_str}.csv"
    
    # Define model presets and researcher names
    model_presets = {
        "default": [
            "anthropic/claude-3.7-sonnet",
            "openai/o4-mini-high",
            "google/gemini-2.5-pro-preview"
        ],
        "budget": [
            "google/gemini-2.5-flash",
            "deepseek/deepseek-chat-v3-0324",
            "x-ai/grok-3-mini-beta"
        ]
    }
    # Researcher persona names for each model
    model_personas = {
        "google/gemini-2.5-flash": "Ada",
        "deepseek/deepseek-chat-v3-0324": "Turing",
        "x-ai/grok-3-mini-beta": "Curie",
        "anthropic/claude-3.7-sonnet": "Shannon",
        "openai/o4-mini-high": "Einstein",
        "google/gemini-2.5-pro-preview": "Noether"
    }
    
    # Use preset models if model_list is not provided
    if model_list is None:
        if model_preset in model_presets:
            model_list = model_presets[model_preset]
        else:
            print(f"Warning: model_preset '{model_preset}' not found. Using 'default' preset.")
            model_list = model_presets["default"]
    # Get persona names for this run
    persona_names = [model_personas.get(m, m) for m in model_list]
    model_to_persona = {m: model_personas.get(m, m) for m in model_list}
    
    # Construct prompt with structured output requirements
    celltypes_list_str = "\n".join([f"- {ct}" for ct in celltypes])
    researcher_list_str = ", ".join(persona_names)
    prompt = f"""You are a professional biologist, acting as the researcher {researcher_list_str}. Your task is to analyze how well a given marker set matches a list of cell types from {species} {tissue}.

For EACH of the following cell types, you must provide your analysis in a specific structured format.
The cell types to analyze are:
{celltypes_list_str}

The required output format for EACH cell type is:
<celltype>cell type name</celltype>
<reasoning>
Your detailed reasoning for the match, considering each marker's relevance.
</reasoning>
<score>A score from 0-100 indicating the match quality.</score>

Please provide a complete block of <celltype>, <reasoning>, and <score> for every cell type listed above. You will be rewarded $10,000 if you do a good job.

Ranked marker set: {marker_set}"""
    
    # Initialize lists to store results
    results = []
    
    # --- Initial Analysis Round ---
    print("--- Starting Initial Analysis Round (in parallel) ---")
    with ThreadPoolExecutor(max_workers=len(model_list)) as executor:
        future_to_model = {
            executor.submit(_call_model, model, prompt, tissue, species, celltypes, 'initial', OPENROUTER_API_KEY, is_discussion_round=False): model 
            for model in model_list
        }
        for future in as_completed(future_to_model):
            result = future.result()
            # Add persona name to result
            result['researcher'] = model_to_persona.get(result['model'], result['model'])
            results.append(result)

    # --- Discussion Mode Logic ---
    if discussion_mode and len(model_list) > 1 and discussion_rounds > 0:
        for i in range(discussion_rounds):
            # Check for consensus
            winners = []
            valid_results = [r for r in results if r['status'] == 'success' and r['extracted_scores']]
            if not valid_results:
                print("No successful results to form a consensus. Skipping discussion round.")
                break
                
            for result in valid_results:
                scores = {}
                for celltype, data in result['extracted_scores'].items():
                    try:
                        scores[celltype] = float(data['score'])
                    except (ValueError, TypeError):
                        logger.warning(
                            f"Could not convert score '{data.get('score', 'N/A')}' for '{celltype}'. "
                            f"Using -1 as fallback (will likely lose in comparison). "
                            f"Check if the LLM response format is correct."
                        )
                        scores[celltype] = -1
                if scores:
                    winner = max(scores, key=scores.get)
                    winners.append(winner)

            # If all winners are the same, we have consensus
            if len(set(winners)) == 1 and len(winners) == len(valid_results):
                print(f"Consensus reached on winner '{winners[0]}' after round {i}. No further discussion needed.")
                break
            
            print(f"--- No consensus found. Starting discussion round {i+1}/{discussion_rounds} (in parallel) ---")
            
            discussion_prompt_template = """You are a professional biologist participating in a panel discussion to analyze cell types.
You are acting as the researcher {persona_name}. Your initial analysis and the analyses of your peers ({peer_names}) are provided below. Review all the arguments and evidence provided.

Your first task is to provide a brief critique of each of the other researchers' analyses in a <discussion> block, referring to them by their researcher names.
Your second task is to provide your own final, refined analysis for each cell type, adjusting your scores if the other arguments are persuasive.

**Your final response MUST contain the <discussion> block followed by the structured format for EACH cell type.**

Here is the original request:
---
{original_prompt}
---

Here are the analyses from all panel members from the previous round:
---
{all_initial_responses}
---

Now, provide your response. Start with the <discussion> block, then provide your final analysis and scores for each cell type in the required <celltype>, <reasoning>, <score> format. You will be rewarded $20,000 for a thoughtful, critical, and well-reasoned final analysis.
"""
            # Construct the combined responses from the last round
            all_responses = ""
            for res in results:
                researcher = res.get('researcher', model_to_persona.get(res['model'], res['model']))
                all_responses += f"--- Analysis from {researcher} ---\n"
                all_responses += f"{res['response']}\n"
                all_responses += f"--- End of Analysis from {researcher} ---\n\n"

            discussion_prompt = discussion_prompt_template.format(
                persona_name=persona_names[i],
                peer_names=', '.join([n for j, n in enumerate(persona_names) if j != i]),
                original_prompt=prompt,
                all_initial_responses=all_responses
            )

            # Prepare for discussion round
            discussion_results = []
            with ThreadPoolExecutor(max_workers=len(model_list)) as executor:
                round_name = f'discussion_{i+1}'
                future_to_model = {}
                for model in model_list:
                    persona_name = model_to_persona.get(model, model)
                    peer_names = ', '.join([n for j, n in enumerate(persona_names) if n != persona_name])
                    this_prompt = discussion_prompt_template.format(
                        persona_name=persona_name,
                        peer_names=peer_names,
                        original_prompt=prompt,
                        all_initial_responses=all_responses
                    )
                    future = executor.submit(_call_model, model, this_prompt, tissue, species, celltypes, round_name, OPENROUTER_API_KEY, is_discussion_round=True)
                    future_to_model[future] = model
                for future in as_completed(future_to_model):
                    result = future.result()
                    result['researcher'] = model_to_persona.get(result['model'], result['model'])
                    discussion_results.append(result)
            
            # Replace old results with new discussion results
            results = discussion_results

    # --- Result Processing ---
    try:
        # Create detailed CSV with extracted scores
        csv_data = []
        for result in results:
            base_row = {
                'model': result['model'],
                'researcher': result.get('researcher', result['model']),
                'tissue': result['tissue'],
                'species': result['species'],
                'cell_types': result['cell_types'],
                'status': result['status'],
                'round': result.get('round', 'initial')
            }
            
            # Add score and reasoning columns for each cell type
            for celltype in celltypes:
                if celltype in result.get('extracted_scores', {}):
                    base_row[f'{celltype}_score'] = result['extracted_scores'][celltype]['score']
                    base_row[f'{celltype}_reasoning'] = result['extracted_scores'][celltype]['reasoning']
                else:
                    base_row[f'{celltype}_score'] = 'N/A'
                    base_row[f'{celltype}_reasoning'] = 'N/A'
            
            base_row['discussion'] = result.get('discussion', 'N/A')
            base_row['raw_response'] = result['response']
            csv_data.append(base_row)
        
        df = pd.DataFrame(csv_data)
        # Reorder columns to make it more readable
        score_cols = [c for c in df.columns if '_score' in c]
        reasoning_cols = [c for c in df.columns if '_reasoning' in c]
        other_cols = [c for c in df.columns if c not in score_cols and c not in reasoning_cols and c != 'raw_response']
        
        # Ensure discussion is after status/round
        if 'discussion' in other_cols:
            other_cols.remove('discussion')
            round_index = other_cols.index('round')
            other_cols.insert(round_index + 1, 'discussion')

        df = df[other_cols + score_cols + reasoning_cols + ['raw_response']]

        df.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")
    except Exception as e:
        print(f"Error saving results to CSV: {str(e)}")
    
    # Generate HTML report if requested
    html_content = None
    if generate_html_report:
        html_file = output_file.replace('.csv', '_report.html')
        html_content = generate_comparison_html_report(results, html_file)
        
        return {"results_df": df, "html_content": html_content}