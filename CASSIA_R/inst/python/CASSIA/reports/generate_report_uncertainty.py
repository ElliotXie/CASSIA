"""
CASSIA Uncertainty Quantification Report Generator
==================================================
Generates HTML reports for uncertainty quantification results from
runCASSIA_n_times_similarity_score().
"""

import os
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
from collections import Counter

# Try to import matplotlib for pie chart generation
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    import base64
    from io import BytesIO
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

import re


def _parse_unified_results(unified_str: str) -> Counter:
    """
    Parse unified results string and return Counter of main cell types.

    Args:
        unified_str: String like "result1:(Plasma Cell, IgA-producing Plasma Cell),result2:..."

    Returns:
        Counter of main cell types
    """
    if not unified_str:
        return Counter()

    # Pattern: result1:(MainType, SubType) or result1:('MainType', 'SubType')
    pattern = r"result\d+:\(?['\"]?([^,'\")]+)"
    matches = re.findall(pattern, unified_str)
    return Counter(matches)


def _parse_unified_results_subtypes(unified_str: str) -> Counter:
    """
    Parse unified results string and return Counter of sub cell types.

    Args:
        unified_str: String like "result1:(Plasma Cell, IgA-producing Plasma Cell),result2:..."

    Returns:
        Counter of sub cell types
    """
    if not unified_str:
        return Counter()

    # Pattern: result1:(MainType, SubType) or result1:('MainType', 'SubType')
    # Capture the sub type (second element after the comma)
    pattern = r"result\d+:\(?['\"]?[^,'\")]+['\"]?,\s*['\"]?([^)'\")]+)"
    matches = re.findall(pattern, unified_str)
    return Counter([m.strip() for m in matches])


def _parse_unified_results_to_table(unified_str: str) -> List[Tuple[str, str, str]]:
    """
    Parse unified results string into table data.

    Args:
        unified_str: String like "result1:(Plasma Cell, IgA-producing),result2:..."

    Returns:
        List of tuples: (iteration_num, main_type, sub_type)
    """
    if not unified_str:
        return []

    results = []
    # Pattern: result1:(MainType, SubType) or result1:('MainType', 'SubType')
    pattern = r"result(\d+):\(?['\"]?([^,'\")]+)['\"]?,\s*['\"]?([^)'\")]+)"
    matches = re.findall(pattern, unified_str)
    for match in matches:
        iter_num, main_type, sub_type = match
        results.append((iter_num, main_type.strip(), sub_type.strip()))
    return results


def _get_score_class(score: float) -> Tuple[str, str, str]:
    """
    Get CSS class and interpretation for a similarity score.

    Args:
        score: Similarity score between 0.0 and 1.0

    Returns:
        Tuple of (css_class, color, interpretation)
    """
    if score >= 0.8:
        return 'score-high', '#22c55e', 'High Confidence'
    elif score >= 0.5:
        return 'score-medium', '#eab308', 'Moderate Confidence'
    else:
        return 'score-low', '#ef4444', 'Low Confidence'


def _check_agreement(main_type: str, sub_type: str,
                     consensus_main: str, consensus_sub: str) -> str:
    """
    Check agreement between a round's result and consensus.

    Returns:
        'full' - both main and sub match
        'partial' - only main matches
        'none' - neither matches
    """
    main_match = main_type.lower().strip() == consensus_main.lower().strip()
    sub_match = sub_type.lower().strip() == consensus_sub.lower().strip()

    if main_match and sub_match:
        return 'full'
    elif main_match:
        return 'partial'
    else:
        return 'none'


def _calculate_summary_stats(original_results: List,
                             consensus_main: str,
                             consensus_sub: str) -> Dict[str, Any]:
    """
    Calculate summary statistics from per-round results.

    Args:
        original_results: List of (main_type, sub_type) tuples
        consensus_main: Consensus main cell type
        consensus_sub: Consensus sub cell type

    Returns:
        Dictionary with summary statistics
    """
    main_types = []
    sub_types = []
    full_agreements = 0
    partial_agreements = 0

    for item in original_results:
        if isinstance(item, tuple) and len(item) >= 2:
            main_type = str(item[0]) if item[0] else 'Unknown'
            sub_type = str(item[1]) if item[1] else 'Unknown'
        else:
            main_type = 'Unknown'
            sub_type = 'Unknown'

        main_types.append(main_type)
        sub_types.append(sub_type)

        agreement = _check_agreement(main_type, sub_type, consensus_main, consensus_sub)
        if agreement == 'full':
            full_agreements += 1
        elif agreement == 'partial':
            partial_agreements += 1

    n_rounds = len(original_results)
    unique_main = len(set(main_types))
    unique_sub = len(set(sub_types))

    # Calculate agreement percentages
    main_agreement = sum(1 for m in main_types
                        if m.lower().strip() == consensus_main.lower().strip())
    sub_agreement = sum(1 for s in sub_types
                       if s.lower().strip() == consensus_sub.lower().strip())

    return {
        'n_rounds': n_rounds,
        'unique_main_types': unique_main,
        'unique_sub_types': unique_sub,
        'main_agreement_pct': (main_agreement / n_rounds * 100) if n_rounds > 0 else 0,
        'sub_agreement_pct': (sub_agreement / n_rounds * 100) if n_rounds > 0 else 0,
        'full_agreements': full_agreements,
        'partial_agreements': partial_agreements,
        'main_type_counts': Counter(main_types),
        'sub_type_counts': Counter(sub_types)
    }


def _generate_pie_chart(type_counts: Counter, title: str) -> str:
    """
    Generate a pie chart as base64 encoded image.

    Args:
        type_counts: Counter of cell type frequencies
        title: Chart title

    Returns:
        Base64 encoded PNG image string, or empty string if matplotlib unavailable
    """
    if not HAS_MATPLOTLIB or not type_counts:
        return ''

    try:
        fig, ax = plt.subplots(figsize=(6, 4))

        labels = list(type_counts.keys())
        sizes = list(type_counts.values())

        # Color palette
        colors = ['#22c55e', '#3b82f6', '#f59e0b', '#ef4444', '#8b5cf6',
                  '#ec4899', '#06b6d4', '#84cc16']

        # Create pie chart
        wedges, texts, autotexts = ax.pie(
            sizes,
            labels=labels,
            autopct='%1.0f%%',
            colors=colors[:len(sizes)],
            startangle=90,
            textprops={'fontsize': 10}
        )

        ax.set_title(title, fontsize=12, fontweight='bold')

        # Convert to base64
        buf = BytesIO()
        plt.tight_layout()
        plt.savefig(buf, format='png', dpi=100, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        plt.close(fig)
        buf.seek(0)

        return base64.b64encode(buf.read()).decode('utf-8')
    except Exception:
        return ''


def _build_rounds_table(original_results: List,
                        consensus_main: str,
                        consensus_sub: str) -> str:
    """
    Build HTML table for per-round results.

    Args:
        original_results: List of results per round
        consensus_main: Consensus main cell type
        consensus_sub: Consensus sub cell type

    Returns:
        HTML table string
    """
    rows = []

    for i, item in enumerate(original_results, 1):
        # Extract data based on result format
        if isinstance(item, tuple):
            if len(item) >= 2:
                main_type = str(item[0]) if item[0] else 'Unknown'
                sub_type = str(item[1]) if item[1] else 'Unknown'
            else:
                main_type = str(item[0]) if item else 'Unknown'
                sub_type = 'Unknown'
        else:
            main_type = str(item) if item else 'Unknown'
            sub_type = 'Unknown'

        # Check agreement
        agreement = _check_agreement(main_type, sub_type, consensus_main, consensus_sub)

        # Set row styling based on agreement
        if agreement == 'full':
            row_class = ''
            icon = '<span style="color: #22c55e; font-size: 18px;">✓</span>'
        elif agreement == 'partial':
            row_class = 'class="row-partial"'
            icon = '<span style="color: #f59e0b; font-size: 18px;">⚠</span>'
        else:
            row_class = 'class="row-disagree"'
            icon = '<span style="color: #ef4444; font-size: 18px;">✗</span>'

        rows.append(f'''
        <tr {row_class}>
            <td style="text-align: center; font-weight: bold;">{i}</td>
            <td>{main_type}</td>
            <td>{sub_type}</td>
        </tr>
        ''')

    return '\n'.join(rows)


def generate_uq_html_report(
    results: Dict[str, Any],
    output_path: str,
    tissue: str = None,
    species: str = None,
    model: str = None,
    n_iterations: int = None,
    marker_list: List[str] = None
) -> str:
    """
    Generate an HTML report for uncertainty quantification results.

    Args:
        results: Output dictionary from runCASSIA_n_times_similarity_score()
        output_path: Path to save the HTML report
        tissue: Tissue type used in analysis
        species: Species used in analysis
        model: Model name used
        n_iterations: Number of iterations run
        marker_list: List of markers analyzed

    Returns:
        Path to the generated HTML report
    """
    # Extract data from results
    # Use LLM Generated Consensus Score (0-100) instead of similarity_score (0-1)
    llm_consensus_score = results.get('llm_generated_consensus_score_llm')
    # Handle None or missing values
    if llm_consensus_score is None or llm_consensus_score == 0:
        similarity_score = results.get('consensus_score_llm') or results.get('similarity_score') or 0
        # Ensure similarity_score is a number
        if similarity_score is None:
            similarity_score = 0
        llm_consensus_score = similarity_score * 100 if similarity_score <= 1 else similarity_score
    # Ensure llm_consensus_score is a number
    if llm_consensus_score is None:
        llm_consensus_score = 0
    consensus_types = results.get('consensus_types', ('Unknown', 'Unknown'))
    general_celltype = results.get('general_celltype_llm', 'Unknown')
    sub_celltype = results.get('sub_celltype_llm', 'Unknown')
    mixed_types = results.get('Possible_mixed_celltypes_llm', [])
    original_results = results.get('original_results', [])
    llm_reasoning = results.get('llm_response', '')  # LLM's reasoning for the consensus score
    # Try unified_results_llm first (batch), then unified_results (single)
    unified_results = results.get('unified_results_llm', '') or results.get('unified_results', '')

    # Handle consensus types
    if isinstance(consensus_types, tuple) and len(consensus_types) >= 2:
        consensus_main = consensus_types[0] or general_celltype or 'Unknown'
        consensus_sub = consensus_types[1] or sub_celltype or 'Unknown'
    else:
        consensus_main = general_celltype or 'Unknown'
        consensus_sub = sub_celltype or 'Unknown'

    # Get score styling - llm_consensus_score is 0-100, divide by 100 for _get_score_class
    score_class, score_color, score_interpretation = _get_score_class(llm_consensus_score / 100)
    score_pct = int(llm_consensus_score)  # Already 0-100

    # Calculate summary stats
    stats = _calculate_summary_stats(original_results, consensus_main, consensus_sub)

    # Generate pie charts using unified results
    unified_main_counts = _parse_unified_results(unified_results)
    unified_sub_counts = _parse_unified_results_subtypes(unified_results)

    # Main cell type pie chart
    pie_chart_main_base64 = _generate_pie_chart(
        unified_main_counts if unified_main_counts else stats['main_type_counts'],
        'Main Cell Type Distribution'
    )
    # Sub cell type pie chart
    pie_chart_sub_base64 = _generate_pie_chart(
        unified_sub_counts if unified_sub_counts else stats['sub_type_counts'],
        'Sub Cell Type Distribution'
    )

    # Build pie charts HTML (side by side)
    if pie_chart_main_base64 or pie_chart_sub_base64:
        pie_chart_html = '<div style="display: flex; flex-wrap: wrap; gap: 20px; justify-content: center;">'
        if pie_chart_main_base64:
            pie_chart_html += f'''
            <div class="chart-container">
                <img src="data:image/png;base64,{pie_chart_main_base64}" alt="Main Cell Type Distribution">
            </div>
            '''
        if pie_chart_sub_base64:
            pie_chart_html += f'''
            <div class="chart-container">
                <img src="data:image/png;base64,{pie_chart_sub_base64}" alt="Sub Cell Type Distribution">
            </div>
            '''
        pie_chart_html += '</div>'
    else:
        pie_chart_html = '<p style="color: #666; font-style: italic;">Pie charts unavailable (matplotlib not installed)</p>'

    # Build rounds table (original results)
    rounds_table = _build_rounds_table(original_results, consensus_main, consensus_sub)

    # Build unified results table
    unified_results_data = _parse_unified_results_to_table(unified_results)
    unified_rows = []
    for iter_num, main_type, sub_type in unified_results_data:
        unified_rows.append(f'''
        <tr>
            <td style="text-align: center; font-weight: bold;">{iter_num}</td>
            <td>{main_type}</td>
            <td>{sub_type}</td>
        </tr>
        ''')
    unified_rounds_table = '\n'.join(unified_rows) if unified_rows else '<tr><td colspan="3">No unified results data</td></tr>'

    # Format mixed types
    mixed_types_html = ', '.join(mixed_types) if mixed_types else '<span style="color: #666;">None detected</span>'

    # Format markers preview
    if marker_list:
        markers_preview = ', '.join(marker_list[:10])
        if len(marker_list) > 10:
            markers_preview += f' ... (+{len(marker_list) - 10} more)'
    else:
        markers_preview = 'Not specified'

    # Format LLM reasoning - clean up and make HTML-safe
    if llm_reasoning:
        # Escape HTML and preserve line breaks
        import html as html_module
        llm_reasoning_html = html_module.escape(str(llm_reasoning)).replace('\n', '<br>')
    else:
        llm_reasoning_html = '<span style="color: #999;">No reasoning available</span>'

    # Build HTML
    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Uncertainty Quantification Report</title>
    <style>
        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }}
        body {{
            font-family: 'Segoe UI', Roboto, -apple-system, sans-serif;
            background: linear-gradient(135deg, #f5f7fa 0%, #e4e8ec 100%);
            min-height: 100vh;
            padding: 20px;
            line-height: 1.6;
        }}
        .container {{
            max-width: 1000px;
            margin: 0 auto;
            background: white;
            border-radius: 16px;
            box-shadow: 0 4px 24px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            text-align: center;
        }}
        .header h1 {{
            font-size: 28px;
            margin-bottom: 10px;
            font-weight: 700;
        }}
        .header-meta {{
            font-size: 14px;
            opacity: 0.9;
        }}
        .header-meta span {{
            margin: 0 10px;
        }}
        .section {{
            padding: 25px 30px;
            border-bottom: 1px solid #e5e7eb;
        }}
        .section:last-child {{
            border-bottom: none;
        }}
        .section-title {{
            font-size: 18px;
            font-weight: 600;
            color: #1f2937;
            margin-bottom: 15px;
            display: flex;
            align-items: center;
            gap: 8px;
        }}
        .top-row {{
            display: flex;
            gap: 20px;
            flex-wrap: wrap;
        }}
        .score-card {{
            flex: 1;
            min-width: 200px;
            background: linear-gradient(145deg, #f8fafc, #e2e8f0);
            border-radius: 12px;
            padding: 25px;
            text-align: center;
        }}
        .score-badge {{
            width: 120px;
            height: 120px;
            border-radius: 50%;
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            margin: 0 auto 15px;
            color: white;
            font-weight: bold;
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
        }}
        .score-high {{
            background: linear-gradient(135deg, #22c55e, #16a34a);
        }}
        .score-medium {{
            background: linear-gradient(135deg, #eab308, #ca8a04);
        }}
        .score-low {{
            background: linear-gradient(135deg, #ef4444, #dc2626);
        }}
        .score-value {{
            font-size: 36px;
            line-height: 1;
        }}
        .score-label {{
            font-size: 12px;
            opacity: 0.9;
            margin-top: 5px;
        }}
        .score-interpretation {{
            font-size: 14px;
            color: #4b5563;
            margin-top: 10px;
        }}
        .consensus-card {{
            flex: 2;
            min-width: 300px;
            background: #f0f9ff;
            border-radius: 12px;
            padding: 25px;
            border-left: 4px solid #3b82f6;
        }}
        .consensus-item {{
            margin-bottom: 12px;
        }}
        .consensus-label {{
            font-size: 12px;
            text-transform: uppercase;
            color: #6b7280;
            font-weight: 600;
            letter-spacing: 0.5px;
        }}
        .consensus-value {{
            font-size: 18px;
            color: #1f2937;
            font-weight: 500;
            margin-top: 4px;
        }}
        .chart-container {{
            text-align: center;
            padding: 20px;
        }}
        .chart-container img {{
            max-width: 100%;
            height: auto;
            border-radius: 8px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 10px;
        }}
        th {{
            background: #f1f5f9;
            padding: 12px;
            text-align: left;
            font-weight: 600;
            color: #374151;
            font-size: 13px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        td {{
            padding: 12px;
            border-bottom: 1px solid #e5e7eb;
            color: #4b5563;
        }}
        tr:hover {{
            background: #f8fafc;
        }}
        .row-partial {{
            background-color: #fffbeb !important;
            border-left: 4px solid #f59e0b;
        }}
        .row-disagree {{
            background-color: #fef2f2 !important;
            border-left: 4px solid #ef4444;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
        }}
        .summary-item {{
            background: #f8fafc;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
        }}
        .summary-value {{
            font-size: 24px;
            font-weight: 700;
            color: #1f2937;
        }}
        .summary-label {{
            font-size: 12px;
            color: #6b7280;
            margin-top: 5px;
        }}
        .progress-bar {{
            height: 8px;
            background: #e5e7eb;
            border-radius: 4px;
            margin-top: 8px;
            overflow: hidden;
        }}
        .progress-fill {{
            height: 100%;
            border-radius: 4px;
            transition: width 0.3s ease;
        }}
        .footer {{
            text-align: center;
            padding: 20px;
            color: #9ca3af;
            font-size: 12px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>Uncertainty Quantification Report</h1>
            <div class="header-meta">
                <span><strong>Tissue:</strong> {tissue or 'N/A'}</span>
                <span>|</span>
                <span><strong>Species:</strong> {species or 'N/A'}</span>
                <span>|</span>
                <span><strong>Model:</strong> {model or 'N/A'}</span>
                <span>|</span>
                <span><strong>Iterations:</strong> {n_iterations or stats['n_rounds']}</span>
            </div>
        </div>

        <div class="section">
            <div class="top-row">
                <div class="score-card">
                    <div class="score-badge {score_class}">
                        <span class="score-value">{score_pct}%</span>
                        <span class="score-label">LLM Consensus</span>
                    </div>
                    <div class="score-interpretation"><strong>LLM Consensus Score:</strong> {score_interpretation}</div>
                </div>

                <div class="consensus-card">
                    <h3 class="section-title">🎯 Consensus Result</h3>
                    <div class="consensus-item">
                        <div class="consensus-label">Main Cell Type</div>
                        <div class="consensus-value">{general_celltype}</div>
                    </div>
                    <div class="consensus-item">
                        <div class="consensus-label">Sub Cell Type</div>
                        <div class="consensus-value">{sub_celltype}</div>
                    </div>
                    <div class="consensus-item">
                        <div class="consensus-label">Possible Mixed Types</div>
                        <div class="consensus-value">{mixed_types_html}</div>
                    </div>
                </div>
            </div>
        </div>

        <div class="section">
            <h3 class="section-title">🤖 LLM Consensus Reasoning</h3>
            <div style="background: #f8fafc; border-radius: 8px; padding: 15px; border-left: 4px solid #667eea; font-size: 14px; line-height: 1.7; color: #374151; max-height: 400px; overflow-y: auto;">
                {llm_reasoning_html}
            </div>
        </div>

        <div class="section">
            <h3 class="section-title">📊 Cell Type Distribution (from Unified Results)</h3>
            {pie_chart_html}
        </div>

        <div class="section">
            <h3 class="section-title">🔄 Per-Round Results (Original)</h3>
            <table>
                <thead>
                    <tr>
                        <th style="width: 60px;">Round</th>
                        <th>Main Cell Type</th>
                        <th>Sub Cell Type</th>
                    </tr>
                </thead>
                <tbody>
                    {rounds_table}
                </tbody>
            </table>
        </div>

        <div class="section">
            <h3 class="section-title">🔄 Per-Round Results (Unified LLM)</h3>
            <table>
                <thead>
                    <tr>
                        <th style="width: 60px;">Round</th>
                        <th>Main Cell Type</th>
                        <th>Sub Cell Type</th>
                    </tr>
                </thead>
                <tbody>
                    {unified_rounds_table}
                </tbody>
            </table>
        </div>

        <div class="section">
            <h3 class="section-title">🧬 Markers Analyzed</h3>
            <p style="color: #4b5563; font-size: 14px;">{markers_preview}</p>
        </div>

        <div class="footer">
            Generated by CASSIA on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
        </div>
    </div>
</body>
</html>
'''

    # Write to file
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)

    print(f"Uncertainty quantification report saved to: {output_path}")
    return output_path


def generate_uq_batch_html_report(
    processed_results: Dict[str, Any],
    organized_results: Dict[str, List],
    output_path: str,
    model: str = None,
    provider: str = None
) -> str:
    """
    Generate an HTML report for batch uncertainty quantification results.

    Args:
        processed_results: Dictionary mapping cluster IDs to their processed results
        organized_results: Dictionary mapping cluster IDs to their original iteration results
        output_path: Path to save the HTML report
        model: Model name used
        provider: Provider used

    Returns:
        Path to the generated HTML report
    """
    # Calculate summary statistics
    total_clusters = len(processed_results)
    scores = []

    for cluster_id, result in processed_results.items():
        if isinstance(result, dict) and 'error' not in result:
            # Use llm_generated_consensus_score_llm (already 0-100 percentage)
            score = result.get('llm_generated_consensus_score_llm', 0)
            if isinstance(score, (int, float)):
                scores.append(score)

    avg_score = sum(scores) / len(scores) if scores else 0
    avg_score_class, _, avg_interpretation = _get_score_class(avg_score / 100)  # Divide by 100 for class

    # Build cluster rows
    cluster_rows = []
    for cluster_id, result in processed_results.items():
        if isinstance(result, dict) and 'error' not in result:
            main_type = result.get('general_celltype_llm', 'Unknown')
            sub_type = result.get('sub_celltype_llm', 'Unknown')
            # Use llm_generated_consensus_score_llm (already 0-100 percentage)
            score = result.get('llm_generated_consensus_score_llm', 0)
            unified_llm = result.get('unified_results_llm', 'N/A')
            unified_oncology = result.get('unified_results_oncology', 'N/A')

            # Get score styling (divide by 100 for class since it expects 0-1)
            score_class, score_color, _ = _get_score_class(score / 100 if score else 0)
            score_pct = int(score)  # Already 0-100, no multiplication needed

            # Build iteration details (original results)
            original = organized_results.get(cluster_id, [])
            iteration_rows = []
            for i, item in enumerate(original, 1):
                if isinstance(item, tuple) and len(item) >= 2:
                    iter_main = str(item[0]) if item[0] else 'N/A'
                    iter_sub = str(item[1]) if item[1] else 'N/A'
                else:
                    iter_main = 'N/A'
                    iter_sub = 'N/A'
                iteration_rows.append(f'<tr><td>{i}</td><td>{iter_main}</td><td>{iter_sub}</td></tr>')

            iterations_table = '\n'.join(iteration_rows) if iteration_rows else '<tr><td colspan="3">No iteration data</td></tr>'

            # Parse unified results into table format
            unified_llm_data = _parse_unified_results_to_table(unified_llm)
            unified_oncology_data = _parse_unified_results_to_table(unified_oncology)

            # Build LLM unified table rows
            unified_llm_rows = []
            for iter_num, main, sub in unified_llm_data:
                unified_llm_rows.append(f'<tr><td>{iter_num}</td><td>{main}</td><td>{sub}</td></tr>')
            unified_llm_table = '\n'.join(unified_llm_rows) if unified_llm_rows else '<tr><td colspan="3">No data</td></tr>'

            # Build Oncology unified table rows
            unified_oncology_rows = []
            for iter_num, main, sub in unified_oncology_data:
                unified_oncology_rows.append(f'<tr><td>{iter_num}</td><td>{main}</td><td>{sub}</td></tr>')
            unified_oncology_table = '\n'.join(unified_oncology_rows) if unified_oncology_rows else '<tr><td colspan="3">No data</td></tr>'

            detail_id = f"cluster_detail_{cluster_id.replace(' ', '_').replace('.', '_')}"

            cluster_rows.append(f'''
            <tr class="cluster-row">
                <td style="font-weight: 600;">{cluster_id}</td>
                <td>{main_type}</td>
                <td>{sub_type}</td>
                <td style="text-align: center;">
                    <span class="score-pill {score_class}" style="background: {score_color};">{score_pct}%</span>
                </td>
                <td style="text-align: center;">
                    <button class="toggle-btn" onclick="toggleClusterDetail('{detail_id}')">Show</button>
                </td>
            </tr>
            <tr id="{detail_id}" class="detail-row" style="display: none;">
                <td colspan="5">
                    <div class="iteration-box">
                        <strong>Per-Iteration Results (Original):</strong>
                        <table class="iteration-table">
                            <thead><tr><th>Iter</th><th>Main Type</th><th>Sub Type</th></tr></thead>
                            <tbody>{iterations_table}</tbody>
                        </table>
                    </div>
                    <div class="iteration-box" style="margin-top: 15px;">
                        <strong>Unified Per-Iteration Results (LLM):</strong>
                        <table class="iteration-table">
                            <thead><tr><th>Iter</th><th>Main Type</th><th>Sub Type</th></tr></thead>
                            <tbody>{unified_llm_table}</tbody>
                        </table>
                    </div>
                    <div class="iteration-box" style="margin-top: 15px;">
                        <strong>Unified Per-Iteration Results (Oncology):</strong>
                        <table class="iteration-table">
                            <thead><tr><th>Iter</th><th>Main Type</th><th>Sub Type</th></tr></thead>
                            <tbody>{unified_oncology_table}</tbody>
                        </table>
                    </div>
                </td>
            </tr>
            ''')
        else:
            # Handle error cases
            error_msg = result.get('error', 'Unknown error') if isinstance(result, dict) else 'Processing failed'
            cluster_rows.append(f'''
            <tr class="cluster-row row-error">
                <td style="font-weight: 600;">{cluster_id}</td>
                <td colspan="3" style="color: #ef4444;">{error_msg}</td>
                <td></td>
            </tr>
            ''')

    cluster_table = '\n'.join(cluster_rows)

    # Build HTML
    html = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Batch Uncertainty Quantification Report</title>
    <style>
        * {{
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }}
        body {{
            font-family: 'Segoe UI', Roboto, -apple-system, sans-serif;
            background: linear-gradient(135deg, #f5f7fa 0%, #e4e8ec 100%);
            min-height: 100vh;
            padding: 20px;
            line-height: 1.6;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 16px;
            box-shadow: 0 4px 24px rgba(0,0,0,0.1);
            overflow: hidden;
        }}
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            text-align: center;
        }}
        .header h1 {{
            font-size: 28px;
            margin-bottom: 10px;
            font-weight: 700;
        }}
        .header-meta {{
            font-size: 14px;
            opacity: 0.9;
        }}
        .header-meta span {{
            margin: 0 10px;
        }}
        .section {{
            padding: 25px 30px;
            border-bottom: 1px solid #e5e7eb;
        }}
        .section:last-child {{
            border-bottom: none;
        }}
        .section-title {{
            font-size: 18px;
            font-weight: 600;
            color: #1f2937;
            margin-bottom: 15px;
            display: flex;
            align-items: center;
            gap: 8px;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 15px;
        }}
        .summary-card {{
            background: linear-gradient(145deg, #f8fafc, #e2e8f0);
            border-radius: 12px;
            padding: 20px;
            text-align: center;
        }}
        .summary-value {{
            font-size: 32px;
            font-weight: 700;
            color: #1f2937;
        }}
        .summary-label {{
            font-size: 12px;
            color: #6b7280;
            margin-top: 5px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        .score-badge {{
            width: 80px;
            height: 80px;
            border-radius: 50%;
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            margin: 0 auto 10px;
            color: white;
            font-weight: bold;
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
        }}
        .score-high {{
            background: linear-gradient(135deg, #22c55e, #16a34a);
        }}
        .score-medium {{
            background: linear-gradient(135deg, #eab308, #ca8a04);
        }}
        .score-low {{
            background: linear-gradient(135deg, #ef4444, #dc2626);
        }}
        .score-pill {{
            display: inline-block;
            padding: 4px 12px;
            border-radius: 12px;
            color: white;
            font-size: 13px;
            font-weight: 600;
        }}
        .confidence-bar {{
            display: flex;
            height: 24px;
            border-radius: 12px;
            overflow: hidden;
            margin-top: 10px;
        }}
        .confidence-segment {{
            display: flex;
            align-items: center;
            justify-content: center;
            color: white;
            font-size: 11px;
            font-weight: 600;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 10px;
        }}
        th {{
            background: #f1f5f9;
            padding: 12px;
            text-align: left;
            font-weight: 600;
            color: #374151;
            font-size: 13px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        td {{
            padding: 12px;
            border-bottom: 1px solid #e5e7eb;
            color: #4b5563;
        }}
        .cluster-row:hover {{
            background: #f8fafc;
        }}
        .row-error {{
            background-color: #fef2f2 !important;
        }}
        .toggle-btn {{
            background: #e0e7ff;
            color: #4338ca;
            border: none;
            padding: 6px 12px;
            border-radius: 6px;
            cursor: pointer;
            font-size: 12px;
            font-weight: 500;
            transition: all 0.2s;
        }}
        .toggle-btn:hover {{
            background: #c7d2fe;
        }}
        .detail-row {{
            background: #f8fafc;
        }}
        .iteration-box {{
            background: white;
            padding: 15px;
            border-radius: 8px;
            border: 1px solid #e5e7eb;
        }}
        .iteration-table {{
            margin-top: 10px;
            font-size: 13px;
        }}
        .iteration-table th {{
            background: #f8fafc;
            padding: 8px;
        }}
        .iteration-table td {{
            padding: 8px;
        }}
        .footer {{
            text-align: center;
            padding: 20px;
            color: #9ca3af;
            font-size: 12px;
        }}
    </style>
    <script>
        function toggleClusterDetail(id) {{
            var row = document.getElementById(id);
            var btn = row.previousElementSibling.querySelector('.toggle-btn');
            if (row.style.display === 'none') {{
                row.style.display = 'table-row';
                btn.textContent = 'Hide';
            }} else {{
                row.style.display = 'none';
                btn.textContent = 'Show';
            }}
        }}
    </script>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>Batch Uncertainty Quantification Report</h1>
            <div class="header-meta">
                <span><strong>Model:</strong> {model or 'N/A'}</span>
                <span>|</span>
                <span><strong>Provider:</strong> {provider or 'N/A'}</span>
                <span>|</span>
                <span><strong>Total Clusters:</strong> {total_clusters}</span>
            </div>
        </div>

        <div class="section">
            <h3 class="section-title">Summary Statistics</h3>
            <div class="summary-grid">
                <div class="summary-card">
                    <div class="summary-value">{total_clusters}</div>
                    <div class="summary-label">Total Clusters</div>
                </div>
                <div class="summary-card">
                    <div class="score-badge {avg_score_class}">
                        <span style="font-size: 24px;">{int(avg_score)}%</span>
                    </div>
                    <div class="summary-label">Avg LLM Consensus Score</div>
                </div>
            </div>
        </div>

        <div class="section">
            <h3 class="section-title">Cluster Results</h3>
            <table>
                <thead>
                    <tr>
                        <th>Cluster ID</th>
                        <th>Main Cell Type</th>
                        <th>Sub Cell Type</th>
                        <th style="text-align: center;">LLM Consensus Score</th>
                        <th style="text-align: center;">Details</th>
                    </tr>
                </thead>
                <tbody>
                    {cluster_table}
                </tbody>
            </table>
        </div>

        <div class="footer">
            Generated by CASSIA on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
        </div>
    </div>
</body>
</html>
'''

    # Write to file
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html)

    print(f"Batch uncertainty quantification report saved to: {output_path}")
    return output_path
