import pandas as pd
import os
import requests
import re
import json
from typing import Dict, List, Tuple


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


def generate_comparison_html_report(all_results: List[Dict], output_file: str = None) -> str:
    """
    Generate an HTML report for cell type comparison results.
    
    Args:
        all_results (list): List of result dictionaries from all models
        output_file (str): Optional output file path
        
    Returns:
        str: HTML content
    """
    if not all_results:
        return "<html><body><h1>No results to display</h1></body></html>"
    
    # Get unique cell types and models
    celltypes = set()
    models = set()
    for result in all_results:
        if 'extracted_scores' in result:
            celltypes.update(result['extracted_scores'].keys())
        models.add(result.get('model', 'Unknown'))
    
    celltypes = sorted(list(celltypes))
    models = sorted(list(models))
    
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
            .summary-section {
                margin-bottom: 40px;
                padding: 20px;
                background: linear-gradient(135deg, #f8f9fa, #e9ecef);
                border-radius: 8px;
            }
            .summary-title {
                font-size: 1.8rem;
                font-weight: bold;
                color: #2c3e50;
                margin-bottom: 20px;
                text-align: center;
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
            .model-section {
                margin-bottom: 40px;
                border: 1px solid #ddd;
                border-radius: 8px;
                overflow: hidden;
                box-shadow: 0 2px 8px rgba(0,0,0,0.1);
            }
            .model-header {
                background: linear-gradient(135deg, #3498db, #2980b9);
                color: white;
                padding: 15px 20px;
                font-size: 1.4rem;
                font-weight: bold;
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
        </script>
    </head>
    <body>
        <div class="container">
            <div class="header">
                <h1 class="title">Cell Type Comparison Report</h1>
                <p class="subtitle">Multi-Model Analysis Results</p>
            </div>
    """
    
    # Add summary table FIRST
    html_content += """
            <div class="summary-section">
                <div class="summary-title">📊 Summary Score Table</div>
                <table class="score-table">
                    <thead>
                        <tr>
                            <th style="text-align: left;">Cell Type</th>
    """
    
    for model in models:
        model_short = model.split('/')[-1] if '/' in model else model
        html_content += f"<th>{model_short}</th>"
    
    html_content += """
                            <th>Average</th>
                            <th>Majority Votes</th>
                        </tr>
                    </thead>
                    <tbody>
    """
    
    # Calculate averages and majority voting
    celltype_scores = {}
    model_winners = []  # Track which celltype won for each model
    
    # Collect all scores for each celltype
    for celltype in celltypes:
        scores = []
        for model in models:
            for result in all_results:
                if result.get('model') == model and 'extracted_scores' in result:
                    if celltype in result['extracted_scores']:
                        score_str = result['extracted_scores'][celltype]['score']
                        try:
                            score_num = float(score_str)
                            scores.append(score_num)
                        except (ValueError, TypeError):
                            pass
                        break
        celltype_scores[celltype] = scores
    
    # Determine winner for each model
    for model in models:
        model_scores = {}
        for result in all_results:
            if result.get('model') == model and 'extracted_scores' in result:
                for celltype in celltypes:
                    if celltype in result['extracted_scores']:
                        score_str = result['extracted_scores'][celltype]['score']
                        try:
                            score_num = float(score_str)
                            model_scores[celltype] = score_num
                        except (ValueError, TypeError):
                            model_scores[celltype] = 0
                break
        
        # Find winner (highest score) for this model
        if model_scores:
            winner = max(model_scores, key=model_scores.get)
            model_winners.append(winner)
    
    # Count majority votes
    majority_votes = {}
    for celltype in celltypes:
        majority_votes[celltype] = model_winners.count(celltype)
    
    for celltype in celltypes:
        # Check if this celltype has majority votes
        votes = majority_votes.get(celltype, 0)
        row_class = "majority-winner" if votes >= len(models) // 2 + 1 else ""
        
        html_content += f'<tr class="{row_class}"><td class="cell-type-name">{celltype}</td>'
        
        # Individual model scores (no highlighting)
        for model in models:
            score = "N/A"
            for result in all_results:
                if result.get('model') == model and 'extracted_scores' in result:
                    if celltype in result['extracted_scores']:
                        score = result['extracted_scores'][celltype]['score']
                        break
            
            html_content += f'<td>{score}</td>'
        
        # Average score (keep highlighting)
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
        
        # Majority votes (keep highlighting)
        vote_class = ""
        if votes >= len(models) // 2 + 1:  # Majority
            vote_class = "high-score"
        elif votes > 0:
            vote_class = "medium-score"
        else:
            vote_class = "low-score"
        
        html_content += f'<td class="{vote_class}"><strong>{votes}/{len(models)}</strong></td>'
        html_content += "</tr>"
    
    html_content += """
                    </tbody>
                </table>
            </div>
    """
    
    # Add detailed reasoning by model SECOND
    html_content += '<div class="details-title">📝 Detailed Reasoning by Model</div>'
    
    for i, result in enumerate(all_results):
        model = result.get('model', 'Unknown')
        extracted_scores = result.get('extracted_scores', {})
        
        html_content += f"""
            <div class="model-section">
                <div class="model-header">🤖 {model}</div>
        """
        
        for celltype in celltypes:
            if celltype in extracted_scores:
                score = extracted_scores[celltype]['score']
                reasoning = extracted_scores[celltype]['reasoning']
                
                # Create short preview (first 10-12 words)
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
                    <button class="toggle-button" id="btn_{i}_{celltype.replace(' ', '_')}" 
                            onclick="toggleReasoning('btn_{i}_{celltype.replace(' ', '_')}', 'reasoning_{i}_{celltype.replace(' ', '_')}')">
                        Show Full Reasoning
                    </button>
                    <div class="reasoning-text" id="reasoning_{i}_{celltype.replace(' ', '_')}" style="display: none;">
                        {reasoning}
                    </div>
                </div>
            """
        
        html_content += """
            </div>
        """
    
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


def compareCelltypes(tissue, celltypes, marker_set, species="human", model_list=None, output_file=None, generate_html_report=True):
    """
    Compare cell types using multiple AI models and generate comparison scores with structured output.
    
    Args:
        tissue (str): Tissue type being analyzed
        celltypes (list): List of 2-4 cell types to compare
        marker_set (str): Ranked marker set for comparison
        species (str): Species being analyzed (default: "human")
        model_list (list): List of models to use for comparison (optional)
        output_file (str): Output CSV filename (optional)
        generate_html_report (bool): Whether to generate HTML report (default: True)
    
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
    
    # Use default models if none provided
    if model_list is None:
        model_list = [
            "anthropic/claude-3.7-sonnet",
            "openai/o4-mini-high",
            "google/gemini-2.5-pro-preview"
        ]
    
    # Construct prompt with structured output requirements
    celltypes_formatted = ', '.join(celltypes)
    prompt = f"""You are a professional biologist. Based on the ranked marker set from {species} {tissue}, analyze and score each of the following cell types: {celltypes_formatted}.

For each cell type, provide your analysis in the following structured format:

<celltype>{celltypes[0]}</celltype>
<reasoning>
Provide detailed reasoning for why this marker set does or does not match this cell type. Consider the specific markers and their relevance to this cell type.
</reasoning>
<score>X</score>

<celltype>{celltypes[1] if len(celltypes) > 1 else celltypes[0]}</celltype>
<reasoning>
Provide detailed reasoning for why this marker set does or does not match this cell type. Consider the specific markers and their relevance to this cell type.
</reasoning>
<score>X</score>

{f'''<celltype>{celltypes[2]}</celltype>
<reasoning>
Provide detailed reasoning for why this marker set does or does not match this cell type. Consider the specific markers and their relevance to this cell type.
</reasoning>
<score>X</score>''' if len(celltypes) > 2 else ''}

{f'''<celltype>{celltypes[3]}</celltype>
<reasoning>
Provide detailed reasoning for why this marker set does or does not match this cell type. Consider the specific markers and their relevance to this cell type.
</reasoning>
<score>X</score>''' if len(celltypes) > 3 else ''}

Score each option from 0-100 based on how well the marker set matches that cell type. You will be rewarded $10,000 if you do a good job.

Ranked marker set: {marker_set}"""
    
    # Initialize lists to store results
    results = []
    processed_models = set()  # Track which models we've already processed
    
    for model in model_list:
        # Skip if we've already processed this model
        if model in processed_models:
            continue
            
        try:
            response = requests.post(
                url="https://openrouter.ai/api/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {OPENROUTER_API_KEY}",
                    "HTTP-Referer": "https://elliotxie.github.io/CASSIA/",
                    "X-Title": "CASSIA",
                    "Content-Type": "application/json"
                },
                json={
                    "model": model,
                    "messages": [
                        {
                            "role": "user",
                            "content": prompt
                        }
                    ]
                }
            )
            
            if response.status_code == 200:
                response_data = response.json()
                model_response = response_data['choices'][0]['message']['content']
                
                # Extract structured scores and reasoning
                extracted_scores = extract_celltype_scores(model_response, celltypes)
                
                # Store result with metadata
                results.append({
                    'model': model,
                    'tissue': tissue,
                    'species': species,
                    'cell_types': ', '.join(celltypes),
                    'response': model_response,
                    'extracted_scores': extracted_scores,
                    'status': 'success'
                })
                print(f"Model: {model}")
                print(f"Extracted scores: {extracted_scores}")
                print("-" * 50)
            else:
                results.append({
                    'model': model,
                    'tissue': tissue,
                    'species': species,
                    'cell_types': ', '.join(celltypes),
                    'response': f"Error: {response.status_code}",
                    'extracted_scores': {},
                    'status': 'error'
                })
                
            processed_models.add(model)  # Mark this model as processed
                
        except Exception as e:
            if model not in processed_models:  # Only add error result if we haven't processed this model
                results.append({
                    'model': model,
                    'tissue': tissue,
                    'species': species,
                    'cell_types': ', '.join(celltypes),
                    'response': f"Exception: {str(e)}",
                    'extracted_scores': {},
                    'status': 'error'
                })
                processed_models.add(model)
    
    # Convert results to DataFrame and save to CSV
    try:
        # Create detailed CSV with extracted scores
        csv_data = []
        for result in results:
            base_row = {
                'model': result['model'],
                'tissue': result['tissue'],
                'species': result['species'],
                'cell_types': result['cell_types'],
                'status': result['status']
            }
            
            # Add score and reasoning columns for each cell type
            for celltype in celltypes:
                if celltype in result.get('extracted_scores', {}):
                    base_row[f'{celltype}_score'] = result['extracted_scores'][celltype]['score']
                    base_row[f'{celltype}_reasoning'] = result['extracted_scores'][celltype]['reasoning']
                else:
                    base_row[f'{celltype}_score'] = 'N/A'
                    base_row[f'{celltype}_reasoning'] = 'N/A'
            
            base_row['raw_response'] = result['response']
            csv_data.append(base_row)
        
        df = pd.DataFrame(csv_data)
        df.to_csv(output_file, index=False)
        print(f"Results saved to {output_file}")
    except Exception as e:
        print(f"Error saving results to CSV: {str(e)}")
    
    # Generate HTML report if requested
    html_content = None
    if generate_html_report:
        html_file = output_file.replace('.csv', '_report.html')
        html_content = generate_comparison_html_report(results, html_file)