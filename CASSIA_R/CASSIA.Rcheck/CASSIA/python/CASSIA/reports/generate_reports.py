import os
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
import seaborn as sns
import base64
from io import BytesIO
import re
import glob
from typing import Dict, Any, List, Optional

def calculate_evaluation_metrics(eval_df: pd.DataFrame, score_col: str = 'score') -> Dict[str, float]:
    """
    Calculate metrics from batch evaluation results.
    
    Args:
        eval_df (pd.DataFrame): DataFrame with evaluation results
        score_col (str): Column name for evaluation scores (0-5 or 0-100 scale)
        
    Returns:
        Dict[str, float]: Dictionary with evaluation metrics
    """
    max_score = eval_df[score_col].max()
    is_similarity_scale = max_score > 10
    metrics = {
        'mean_score': eval_df[score_col].mean(),
        'median_score': eval_df[score_col].median(),
        'min_score': eval_df[score_col].min(),
        'max_score': eval_df[score_col].max(),
        'std_score': eval_df[score_col].std(),
        'count': len(eval_df),
    }
    if not is_similarity_scale:
        metrics.update({
            'perfect_ratio': (eval_df[score_col] == 5).mean(),
            'very_good_ratio': (eval_df[score_col] == 4).mean(),
            'good_ratio': (eval_df[score_col] == 3).mean(),
            'partial_ratio': (eval_df[score_col] == 2).mean(),
            'poor_ratio': (eval_df[score_col] == 1).mean(),
            'nonsensical_ratio': (eval_df[score_col] == 0).mean(),
        })
    return metrics

def generate_html_report(result_df: pd.DataFrame, 
                         gold_col: str, 
                         pred_col: str, 
                         score_col: str = "score", 
                         reasoning_col: str = "reasoning",
                         metrics: Optional[Dict[str, float]] = None, 
                         html_report_path: str = "report.html",
                         model_name: str = None) -> None:
    """Generate an HTML report for the evaluation results."""
    
    # Get score range info to determine reporting approach
    min_score = result_df[score_col].min()
    max_score = result_df[score_col].max()
    
    # Determine if we're using rule-based (0-5) or similarity (0-100) scoring
    is_similarity_scale = max_score > 10  # Simple heuristic to detect similarity scale
    
    # Generate a larger, more beautiful histogram and encode as base64 (do not show on screen)
    buf1 = BytesIO()
    plt.figure(figsize=(14, 7))
    sns.set_theme(style="whitegrid")
    palette = sns.color_palette("crest", 6)
    
    if is_similarity_scale:
        # For continuous 0-100 similarity scale
        bins = [0, 20, 40, 60, 80, 100]
        bin_labels = ['0-20', '20-40', '40-60', '60-80', '80-100']
        
        # Create histogram with continuous bins
        ax = sns.histplot(result_df[score_col], bins=bins, kde=False, color=palette[3], edgecolor='black')
        plt.title(f"Distribution of Similarity Scores" + (f" - {model_name}" if model_name else ""), 
                  fontsize=22, fontweight='bold')
        plt.xlabel("Similarity Score", fontsize=18)
        plt.ylabel("Count", fontsize=18)
        plt.xticks([10, 30, 50, 70, 90], bin_labels, fontsize=16)
    else:
        # For discrete 0-5 scale
        ax = sns.histplot(result_df[score_col], bins=[-0.5,0.5,1.5,2.5,3.5,4.5,5.5], 
                          kde=False, discrete=True, color=palette[3], edgecolor='black')
        plt.title(f"Distribution of Evaluation Scores" + (f" - {model_name}" if model_name else ""), 
                  fontsize=22, fontweight='bold')
        plt.xlabel("Score", fontsize=18)
        plt.ylabel("Count", fontsize=18)
        plt.xticks([0,1,2,3,4,5], fontsize=16)
    
    plt.yticks(fontsize=16)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    for patch in ax.patches:
        ax.annotate(int(patch.get_height()), (patch.get_x() + patch.get_width() / 2, patch.get_height()),
                    ha='center', va='bottom', fontsize=14, color='black', fontweight='bold')
    plt.tight_layout()
    plt.savefig(buf1, format="png")
    plt.close()
    buf1.seek(0)
    hist_img = base64.b64encode(buf1.read()).decode("utf-8")

    # Calculate metrics if not provided
    if metrics is None:
        metrics = calculate_evaluation_metrics(result_df, score_col=score_col)
    
    # Calculate score ratio - adjust for similarity scale if needed
    score_sum = result_df[score_col].sum()
    score_count = result_df[score_col].count()
    max_possible = 100 if is_similarity_scale else 5
    score_ratio = (score_sum / (score_count * max_possible)) if score_count > 0 else 0
    
    # Format metrics HTML
    if is_similarity_scale:
        metrics_html = "<ul>" + "".join([
            f"<li><b>{k.replace('_',' ').capitalize()}:</b> {v:.3f}</li>" for k, v in metrics.items() 
            if k in ['mean_score','median_score','min_score','max_score','std_score','count']
        ]) + f"<li><b>Score Ratio:</b> {score_ratio*100:.1f}%</li></ul>"
    else:
        metrics_html = "<ul>" + "".join([
            f"<li><b>{k.replace('_',' ').capitalize()}:</b> {v:.3f}</li>" for k, v in metrics.items()
        ]) + f"<li><b>Score Ratio:</b> {score_ratio*100:.1f}%</li></ul>"

    # Add tissue and species breakdown if available
    additional_breakdowns = ""
    if 'Tissue' in result_df.columns and 'Species' in result_df.columns:
        # Create tissue+species breakdown
        grouped = result_df.groupby(['Tissue', 'Species'])[score_col].mean().reset_index()
        grouped_html = "<table border='1' cellpadding='4' cellspacing='0' style='border-collapse:collapse;'>"
        grouped_html += "<tr><th>Tissue</th><th>Species</th><th>Average Score</th></tr>"
        for _, row in grouped.iterrows():
            grouped_html += f"<tr><td>{row['Tissue']}</td><td>{row['Species']}</td><td>{row[score_col]:.2f}</td></tr>"
        grouped_html += "</table>"
        additional_breakdowns = f"""
        <div class="section">
            <h2>Breakdown by Tissue and Species</h2>
            {grouped_html}
        </div>
        """

    # Sample results: handle differently for similarity vs discrete scale
    sample_rows = []
    if is_similarity_scale:
        # For similarity, show samples from each bin
        bins = [0, 20, 40, 60, 80, 100]
        for i in range(len(bins)-1):
            bin_low, bin_high = bins[i], bins[i+1]
            bin_df = result_df[(result_df[score_col] >= bin_low) & (result_df[score_col] < bin_high)]
            n = min(5, len(bin_df))
            if n > 0:
                bin_desc = f"{bin_low}-{bin_high}"
                for _, row in bin_df.head(n).iterrows():
                    gold = row[gold_col] if gold_col in row else ''
                    pred = row[pred_col] if pred_col in row else ''
                    scr = row[score_col] if score_col in row else ''
                    expl = row[reasoning_col] if reasoning_col in row else ''
                    sample_rows.append(f'<tr><td>{gold}</td><td>{pred}</td><td>{scr:.1f} (bin: {bin_desc})</td><td>{expl[:200]}...</td></tr>')
    else:
        # For discrete, show samples for each score value
        for score in range(6):
            score_df = result_df[result_df[score_col] == score]
            n = min(5, len(score_df))
            if n > 0:
                for _, row in score_df.head(n).iterrows():
                    gold = row[gold_col] if gold_col in row else ''
                    pred = row[pred_col] if pred_col in row else ''
                    scr = row[score_col] if score_col in row else ''
                    expl = row[reasoning_col] if reasoning_col in row else ''
                    sample_rows.append(f'<tr><td>{gold}</td><td>{pred}</td><td>{scr}</td><td>{expl[:200]}...</td></tr>')

    # HTML content
    html = f"""
    <html>
    <head>
        <title>LLM Celltype Annotation Evaluation Report{" - " + model_name if model_name else ""}</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; }}
            h1 {{ color: #2c3e50; }}
            h2 {{ color: #34495e; }}
            .section {{ margin-bottom: 30px; }}
            .metrics {{ background: #f8f8f8; padding: 15px; border-radius: 8px; }}
            .img-container {{ display: flex; gap: 40px; }}
            .img-container img {{ border: 1px solid #ccc; border-radius: 8px; background: #fff; }}
        </style>
    </head>
    <body>
        <h1>LLM Celltype Annotation Evaluation Report{" - " + model_name if model_name else ""}</h1>
        <div class="section">
            <b>Generated:</b> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}<br>
            <b>Total Samples:</b> {len(result_df)}<br>
            <b>Score Type:</b> {"Similarity (0-100)" if is_similarity_scale else "Discrete (0-5)"}<br>
        </div>
        <div class="section metrics">
            <h2>Summary Metrics</h2>
            {metrics_html}
        </div>
        {additional_breakdowns}
        <div class="section">
            <h2>Score Distribution</h2>
            <div class="img-container">
                <div><img src="data:image/png;base64,{hist_img}" width="800"><br>Histogram</div>
            </div>
        </div>
        <div class="section">
            <h2>Sample Results</h2>
            <table border="1" cellpadding="4" cellspacing="0" style="border-collapse:collapse;">
                <tr>
                    <th>Gold Standard</th>
                    <th>Prediction</th>
                    <th>Score</th>
                    <th>Explanation</th>
                </tr>
                {''.join(sample_rows)}
            </table>
            <p><i>Showing up to 5 examples for each {"score bin" if is_similarity_scale else "score"} ({"0-20, 20-40, etc." if is_similarity_scale else "0-5"}).</i></p>
        </div>
    </body>
    </html>
    """
    with open(html_report_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"HTML report saved to {html_report_path}")

def find_evaluation_csvs(root_dir: str) -> List[str]:
    """Find all evaluation CSV files in the given directory and its subdirectories."""
    pattern = os.path.join(root_dir, "**", "*_evaluation.csv")
    return glob.glob(pattern, recursive=True)

def extract_model_name(file_path: str) -> str:
    """Extract model name from the file path."""
    # Extract the filename without extension
    filename = os.path.basename(file_path)
    # Remove '_evaluation.csv' suffix
    if '_evaluation.csv' in filename:
        model_name = filename.replace('_evaluation.csv', '')
        model_name = model_name.replace('combined_scores_', '')
        return model_name
    return "Unknown Model"

def generate_subclustering_report(csv_path, html_report_path=None, model_name=None):
    """
    Generate a beautiful HTML report for subclustering batch results, showing annotation, reasoning, and top marker.
    """
    df = pd.read_csv(csv_path)
    if html_report_path is None:
        html_report_path = csv_path.replace('.csv', '.html')
    if model_name is None:
        model_name = os.path.basename(csv_path).replace('.csv', '')

    # Build table rows with popup markers
    rows = []
    for idx, row in df.iterrows():
        cluster = row.get('Result ID', '')
        main_type = row.get('main_cell_type', '')
        sub_type = row.get('sub_cell_type', '')
        key_markers = row.get('key_markers', '')
        reason = row.get('reason', '')
        # Escape quotes for JavaScript
        escaped_markers = str(key_markers).replace('\\', '\\\\').replace("'", "\\'").replace('"', '&quot;').replace('\n', '<br>')
        escaped_cluster = str(cluster).replace("'", "\\'")
        rows.append(f'''
        <tr>
            <td class="cluster-col">{cluster}</td>
            <td class="annotation-col">
                <div class="main-type">{main_type}</div>
                <div class="sub-type">{sub_type}</div>
            </td>
            <td class="marker-col">
                <button class="marker-toggle" onclick="showMarkerPopup('{escaped_cluster}', '{escaped_markers}')">Show Markers</button>
            </td>
            <td class="reasoning-col"><div class="reasoning-box">{reason}</div></td>
        </tr>
        ''')

    html = f'''
    <html>
    <head>
        <title>Subclustering Annotation Report - {model_name}</title>
        <style>
            body {{ font-family: 'Segoe UI', Arial, sans-serif; background: #f7f9fa; margin: 0; }}
            .container {{ max-width: 1200px; margin: 40px auto; background: #fff; border-radius: 16px; box-shadow: 0 4px 24px rgba(0,0,0,0.08); padding: 32px 32px 24px 32px; }}
            h1 {{ color: #1976d2; margin-bottom: 8px; }}
            .meta {{ color: #555; font-size: 15px; margin-bottom: 24px; }}
            table {{ width: 100%; border-collapse: separate; border-spacing: 0; }}
            th, td {{ padding: 12px 10px; text-align: left; }}
            th {{ background: #1976d2; color: #fff; font-size: 16px; font-weight: 600; border-top-left-radius: 8px; border-top-right-radius: 8px; }}
            tr {{ background: #fff; transition: background 0.2s; }}
            tr:nth-child(even) {{ background: #f3f6fa; }}
            tr:hover {{ background: #e3f2fd; }}
            .cluster-col {{ width: 5%; font-weight: bold; color: #1976d2; font-size: 18px; }}
            .annotation-col {{ width: 20%; }}
            .main-type {{ font-weight: bold; color: #222; font-size: 16px; }}
            .sub-type {{ color: #888; font-size: 14px; margin-top: 2px; }}
            .marker-col {{ width: 20%; font-size: 13px; color: #888; }}
            .marker-toggle {{ background: #e3f2fd; color: #1976d2; border: none; border-radius: 5px; padding: 3px 10px; font-size: 13px; cursor: pointer; }}
            .marker-toggle:hover {{ background: #bbdefb; }}
            .reasoning-col {{ width: 55%; }}
            /* Modal popup styles */
            .modal-overlay {{ display: none; position: fixed; top: 0; left: 0; width: 100%; height: 100%; background: rgba(0,0,0,0.5); z-index: 1000; justify-content: center; align-items: center; }}
            .modal-overlay.active {{ display: flex; }}
            .modal-content {{ background: #fff; border-radius: 12px; padding: 24px; max-width: 600px; max-height: 80vh; overflow-y: auto; box-shadow: 0 8px 32px rgba(0,0,0,0.2); position: relative; }}
            .modal-header {{ display: flex; justify-content: space-between; align-items: center; margin-bottom: 16px; border-bottom: 1px solid #e0e0e0; padding-bottom: 12px; }}
            .modal-title {{ font-size: 18px; font-weight: 600; color: #1976d2; margin: 0; }}
            .modal-close {{ background: #f44336; color: #fff; border: none; border-radius: 50%; width: 28px; height: 28px; font-size: 18px; cursor: pointer; display: flex; align-items: center; justify-content: center; }}
            .modal-close:hover {{ background: #d32f2f; }}
            .modal-body {{ font-size: 14px; color: #444; line-height: 1.6; }}
            .reasoning-box {{ background: #fffde7; border-left: 5px solid #ffe082; border-radius: 7px; padding: 12px 16px; font-size: 15px; color: #444; box-shadow: 0 1px 4px rgba(255,193,7,0.07); }}
        </style>
        <script>
        function showMarkerPopup(cluster, markers) {{
            document.getElementById('modal-cluster').textContent = 'Cluster ' + cluster + ' - Key Markers';
            document.getElementById('modal-markers').innerHTML = markers;
            document.getElementById('marker-modal').classList.add('active');
        }}
        function closeMarkerPopup() {{
            document.getElementById('marker-modal').classList.remove('active');
        }}
        // Close modal when clicking outside
        document.addEventListener('click', function(e) {{
            var modal = document.getElementById('marker-modal');
            if (e.target === modal) {{
                closeMarkerPopup();
            }}
        }});
        // Close modal with Escape key
        document.addEventListener('keydown', function(e) {{
            if (e.key === 'Escape') {{
                closeMarkerPopup();
            }}
        }});
        </script>
    </head>
    <body>
        <div class="container">
            <h1>Subclustering Annotation Report</h1>
            <div class="meta"><b>Model:</b> {model_name} &nbsp; | &nbsp; <b>Generated:</b> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
            <table>
                <tr>
                    <th>Cluster</th>
                    <th>Annotation<br><span style="font-weight:normal;font-size:12px">(Main / Subtype)</span></th>
                    <th>Top Markers</th>
                    <th>Reasoning</th>
                </tr>
                {''.join(rows)}
            </table>
        </div>
        <!-- Modal Popup -->
        <div id="marker-modal" class="modal-overlay">
            <div class="modal-content">
                <div class="modal-header">
                    <h3 id="modal-cluster" class="modal-title">Key Markers</h3>
                    <button class="modal-close" onclick="closeMarkerPopup()">&times;</button>
                </div>
                <div id="modal-markers" class="modal-body"></div>
            </div>
        </div>
    </body>
    </html>
    '''
    with open(html_report_path, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f"Subclustering HTML report saved to {html_report_path}")

def process_evaluation_csv(csv_path: str, overwrite: bool = False, model_name: str = None) -> None:
    try:
        if not os.path.exists(csv_path):
            print(f"File not found: {csv_path}")
            return
        if model_name is None:
            model_name = extract_model_name(csv_path)
        html_path = csv_path.replace('.csv', '.html')
        if os.path.exists(html_path) and not overwrite:
            print(f"HTML report already exists for {model_name}. Skipping. Use --overwrite to regenerate.")
            return
        df = pd.read_csv(csv_path)
        # If subclustering format, use the new report
        expected_cols = {'Result ID', 'main_cell_type', 'sub_cell_type', 'key_markers', 'reason'}
        if expected_cols.issubset(set(df.columns)):
            generate_subclustering_report(csv_path, html_report_path=html_path, model_name=model_name)
            return
        # Determine column names (check new names first, then fall back to old names for backward compatibility)
        gold_col = "Cluster ID" if "Cluster ID" in df.columns else ("True Cell Type" if "True Cell Type" in df.columns else "gold_standard")
        pred_col = "Predicted Detailed Cell Type" if "Predicted Detailed Cell Type" in df.columns else ("Predicted Sub Cell Types" if "Predicted Sub Cell Types" in df.columns else "predicted_celltype")
        
        # For score and reasoning columns, check multiple possibilities
        score_col = None
        for col in ["score", "evaluation_score", "similarity_score"]:
            if col in df.columns:
                score_col = col
                break
        
        reasoning_col = None
        for col in ["reasoning", "evaluation_explanation", "similarity_reasoning", "explanation"]:
            if col in df.columns:
                reasoning_col = col
                break
        
        if not score_col:
            print(f"No score column found in {csv_path}. Skipping.")
            return
            
        if not reasoning_col:
            print(f"No reasoning column found in {csv_path}. Will generate report without reasoning.")
            reasoning_col = None
        
        # Calculate metrics
        metrics = calculate_evaluation_metrics(df, score_col=score_col)
        
        # Generate HTML report
        generate_html_report(
            result_df=df,
            gold_col=gold_col,
            pred_col=pred_col,
            score_col=score_col,
            reasoning_col=reasoning_col,
            metrics=metrics,
            html_report_path=html_path,
            model_name=model_name
        )
        
        print(f"Processed {model_name} - Mean score: {metrics['mean_score']:.2f}")
        
    except Exception as e:
        print(f"Error processing {csv_path}: {str(e)}")

def create_index_html(csv_files: List[str], output_dir: str) -> None:
    """Create an index.html file that links to all the reports."""
    reports = []
    
    for csv_file in csv_files:
        try:
            # Read the CSV to get metrics
            df = pd.read_csv(csv_file)
            
            # Determine score column
            score_col = None
            for col in ["score", "evaluation_score", "similarity_score"]:
                if col in df.columns:
                    score_col = col
                    break
            
            if not score_col:
                continue
                
            # Calculate mean score
            mean_score = df[score_col].mean()
            
            # Get model name and HTML path
            model_name = extract_model_name(csv_file)
            html_path = csv_file.replace('.csv', '.html')
            rel_path = os.path.relpath(html_path, output_dir)
            
            # Add to reports list
            reports.append({
                'model_name': model_name,
                'mean_score': mean_score,
                'html_path': rel_path,
                'count': len(df)
            })
            
        except Exception as e:
            print(f"Error processing {csv_file} for index: {str(e)}")
    
    # Sort reports by mean score (descending)
    reports.sort(key=lambda x: x['mean_score'], reverse=True)
    
    # Create index.html content
    rows = []
    for i, report in enumerate(reports):
        rows.append(f'''
        <tr>
            <td>{i+1}</td>
            <td><a href="{report['html_path']}">{report['model_name']}</a></td>
            <td>{report['mean_score']:.2f}</td>
            <td>{report['count']}</td>
        </tr>
        ''')
    
    html = f'''
    <html>
    <head>
        <title>LLM Celltype Annotation Evaluation - Summary</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; }}
            h1 {{ color: #2c3e50; }}
            table {{ border-collapse: collapse; width: 100%; }}
            th, td {{ padding: 8px; text-align: left; border-bottom: 1px solid #ddd; }}
            tr:hover {{ background-color: #f5f5f5; }}
            th {{ background-color: #4CAF50; color: white; }}
        </style>
    </head>
    <body>
        <h1>LLM Celltype Annotation Evaluation - Summary</h1>
        <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        <table>
            <tr>
                <th>Rank</th>
                <th>Model</th>
                <th>Average Score</th>
                <th>Sample Count</th>
            </tr>
            {''.join(rows)}
        </table>
    </body>
    </html>
    '''
    
    # Write index.html
    index_path = os.path.join(output_dir, 'index.html')
    with open(index_path, 'w', encoding='utf-8') as f:
        f.write(html)
    
    print(f"Index page created at {index_path}")


# =============================================================================
# Score Report HTML Generation Functions
# These functions generate HTML reports from scored annotation results
# =============================================================================

def generate_analysis_html_report(analysis_text):
    """
    Generate an HTML report from conversation analysis text.

    This parses the conversation history and creates a styled HTML report
    with sections for annotation, validation, formatting, and scoring.

    Args:
        analysis_text (str): Pipe-separated conversation text from CASSIA analysis

    Returns:
        str: Complete HTML document as a string
    """
    import json

    # Split the text into sections based on agents
    sections = analysis_text.split(" | ")

    # HTML template with CSS styling - note the double curly braces for CSS
    html_template = """
    <!DOCTYPE html>
    <html>
    <head>
        <style>
            body {{
                font-family: 'Segoe UI', Roboto, -apple-system, sans-serif;
                max-width: 1200px;
                margin: 0 auto;
                padding: 20px;
                background-color: #f0f2f5;
                line-height: 1.6;
            }}
            .container {{
                background-color: white;
                padding: 40px;
                border-radius: 16px;
                box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            }}
            .agent-section {{
                margin-bottom: 35px;
                padding: 25px;
                border-radius: 12px;
                transition: all 0.3s ease;
            }}
            .agent-section:hover {{
                transform: translateY(-2px);
                box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            }}
            .final-annotation {{
                background-color: #f0f7ff;
                border-left: 5px solid #2196f3;
            }}
            .validator {{
                background-color: #f0fdf4;
                border-left: 5px solid #22c55e;
            }}
            .formatting {{
                background: linear-gradient(145deg, #fff7ed, #ffe4c4);
                border-left: 5px solid #f97316;
                box-shadow: 0 4px 15px rgba(249, 115, 22, 0.1);
            }}
            h2 {{
                color: #1a2b3c;
                margin-top: 0;
                font-size: 1.5rem;
                font-weight: 600;
                display: flex;
                align-items: center;
                gap: 10px;
            }}
            ul {{
                margin: 15px 0;
                padding-left: 20px;
            }}
            pre {{
                background-color: #f8fafc;
                padding: 20px;
                border-radius: 8px;
                overflow-x: auto;
                font-family: 'Consolas', 'Monaco', monospace;
                font-size: 0.9rem;
                line-height: 1.5;
            }}
            .validation-result {{
                font-weight: 600;
                color: #16a34a;
                padding: 12px 20px;
                background-color: #dcfce7;
                border-radius: 8px;
                display: inline-block;
                margin: 10px 0;
            }}
            br {{
                margin-bottom: 8px;
            }}
            p {{
                margin: 12px 0;
                color: #374151;
            }}
            .summary-content {{
                display: flex;
                flex-direction: column;
                gap: 24px;
            }}
            .summary-item {{
                display: flex;
                flex-direction: column;
                gap: 8px;
                background: rgba(255, 255, 255, 0.7);
                padding: 16px;
                border-radius: 12px;
                backdrop-filter: blur(8px);
                box-shadow: 0 2px 8px rgba(0, 0, 0, 0.05);
            }}
            .summary-label {{
                font-weight: 600;
                color: #c2410c;
                font-size: 0.95rem;
                text-transform: uppercase;
                letter-spacing: 0.5px;
            }}
            .summary-value {{
                color: #1f2937;
                font-size: 1.1rem;
                padding: 8px 16px;
                background-color: rgba(255, 255, 255, 0.9);
                border-radius: 8px;
                display: inline-block;
                box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
            }}
            .summary-list {{
                margin: 0;
                padding-left: 24px;
                list-style-type: none;
            }}
            .summary-list li {{
                color: #1f2937;
                padding: 8px 0;
                position: relative;
            }}
            .summary-list li:before {{
                content: "•";
                color: #f97316;
                font-weight: bold;
                position: absolute;
                left: -20px;
            }}
            .report-header {{
                text-align: center;
                margin-bottom: 40px;
                padding-bottom: 30px;
                border-bottom: 2px solid rgba(249, 115, 22, 0.2);
            }}
            .report-title {{
                font-size: 2.5rem;
                font-weight: 800;
                color: #1a2b3c;
                margin: 0;
                padding: 0;
                background: linear-gradient(135deg, #f97316, #c2410c);
                -webkit-background-clip: text;
                -webkit-text-fill-color: transparent;
                letter-spacing: -0.5px;
            }}
            .report-subtitle {{
                font-size: 1.1rem;
                color: #64748b;
                margin-top: 8px;
                font-weight: 500;
            }}
            .scoring {{
                background: linear-gradient(145deg, #f0fdf4, #dcfce7);
                border-left: 5px solid #22c55e;
                box-shadow: 0 4px 15px rgba(34, 197, 94, 0.1);
            }}
            .scoring-content {{
                display: flex;
                flex-direction: column;
                gap: 16px;
                color: #1f2937;
                line-height: 1.8;
            }}
            .scoring-content br + br {{
                content: "";
                display: block;
                margin: 12px 0;
            }}
            .empty-list {{
                color: #6b7280;
                font-style: italic;
            }}
            .error-message {{
                color: #dc2626;
                padding: 12px;
                background-color: #fef2f2;
                border-radius: 6px;
                border-left: 4px solid #dc2626;
            }}
            .score-badge {{
                background: linear-gradient(135deg, #22c55e, #16a34a);
                color: white;
                padding: 8px 16px;
                border-radius: 12px;
                font-size: 1.5rem;
                font-weight: 700;
                display: inline-block;
                margin: 12px 0;
                box-shadow: 0 4px 12px rgba(34, 197, 94, 0.2);
                position: relative;
                top: -10px;
            }}
            .score-badge::before {{
                content: "Score:";
                font-size: 0.9rem;
                font-weight: 500;
                margin-right: 8px;
                opacity: 0.9;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="report-header">
                <h1 class="report-title">CASSIA Analysis Report</h1>
                <p class="report-subtitle">Comprehensive Cell Type Analysis and Annotation</p>
            </div>
            {0}
        </div>
    </body>
    </html>
    """

    content = []

    # Process each section
    for section in sections:
        if section.startswith("Final Annotation Agent:"):
            annotation_content = section.replace("Final Annotation Agent:", "").strip()
            content.append("""
                <div class="agent-section final-annotation">
                    <h2>Final Annotation Analysis</h2>
                    {0}
                </div>
            """.format(annotation_content.replace('\n', '<br>')))

        elif section.startswith("Coupling Validator:"):
            validator_content = section.replace("Coupling Validator:", "").strip()
            validation_result = '<div class="validation-result">VALIDATION PASSED</div>' if "VALIDATION PASSED" in validator_content else ""

            content.append("""
                <div class="agent-section validator">
                    <h2>Validation Check</h2>
                    {0}
                    {1}
                </div>
            """.format(validation_result, validator_content.replace('\n', '<br>')))

        elif section.startswith("Formatting Agent:"):
            try:
                # Get the content after "Formatting Agent:"
                json_text = section.replace("Formatting Agent:", "").strip()

                # Since the JSON is consistently formatted with newlines,
                # we can find where it ends (the last '}' followed by a newline or end of string)
                json_end = json_text.rfind('}')
                if json_end != -1:
                    json_content = json_text[:json_end + 1]
                    data = json.loads(json_content)

                    # Process the data...
                    main_cell_type = data.get('main_cell_type', 'Not specified')
                    sub_cell_types = data.get('sub_cell_types', [])
                    mixed_types = data.get('possible_mixed_cell_types', [])
                    num_markers = data.get('num_markers', 'Not specified')

                    # Format the content...
                    formatted_content = f"""
                        <div class="summary-content">
                            <div class="summary-item">
                                <span class="summary-label">Main Cell Type:</span>
                                <span class="summary-value">{main_cell_type}</span>
                            </div>

                            <div class="summary-item">
                                <span class="summary-label">Sub Cell Types:</span>
                                <ul class="summary-list">
                                    {"".join(f'<li>{item}</li>' for item in sub_cell_types) if sub_cell_types
                                     else '<li class="empty-list">No sub cell types identified</li>'}
                                </ul>
                            </div>

                            <div class="summary-item">
                                <span class="summary-label">Possible Mixed Cell Types:</span>
                                <ul class="summary-list">
                                    {"".join(f'<li>{item}</li>' for item in mixed_types) if mixed_types
                                     else '<li class="empty-list">No mixed cell types identified</li>'}
                                </ul>
                            </div>

                            <div class="summary-item">
                                <span class="summary-label">Number of Markers:</span>
                                <span class="summary-value">{num_markers}</span>
                            </div>
                        </div>
                    """

                    content.append(f"""
                        <div class="agent-section formatting">
                            <h2>Summary</h2>
                            {formatted_content}
                        </div>
                    """)
                else:
                    raise ValueError("Could not find JSON content")

            except Exception as e:
                content.append(f"""
                    <div class="agent-section formatting">
                        <h2>Summary</h2>
                        <p class="error-message">Error formatting data: {str(e)}</p>
                    </div>
                """)
        elif section.startswith("Scoring Agent:"):
            try:
                # Get the content after "Scoring Agent:"
                scoring_text = section.split("Scoring Agent:", 1)[1].strip()

                # Split the score from the main text
                main_text, score = scoring_text.rsplit("Score:", 1)
                score = score.strip()

                content.append(r"""
                    <div class="agent-section scoring">
                        <h2>Quality Assessment</h2>
                        <div class="score-badge">{0}</div>
                        <div class="scoring-content">
                            {1}
                        </div>
                    </div>
                """.format(score, main_text.replace('\n', '<br>')))
            except Exception as e:
                content.append(r"""
                    <div class="agent-section scoring">
                        <h2>Quality Assessment</h2>
                        <p class="error-message">Error formatting scoring data: {0}</p>
                    </div>
                """.format(str(e)))

    # Combine all sections
    final_html = html_template.format(''.join(content))
    return final_html


def process_single_report(text, score_reasoning, score):
    """
    Process a single report by combining text with scoring information.

    Args:
        text (str): Conversation history text
        score_reasoning (str): Reasoning from the scoring agent
        score: The numerical score

    Returns:
        str: HTML report as a string
    """
    combined = (
        f"{text}\n"
        f" | Scoring Agent: {score_reasoning}\n"
        f"Score: {score}"
    )
    return generate_analysis_html_report(combined)


def generate_score_index_page(report_files):
    """
    Generate an index HTML page linking to all score reports.

    Args:
        report_files (list): List of report filenames

    Returns:
        str: HTML index page as a string
    """
    index_template = """
    <!DOCTYPE html>
    <html>
    <head>
        <style>
            body {{
                font-family: 'Segoe UI', Roboto, -apple-system, sans-serif;
                max-width: 1200px;
                margin: 0 auto;
                padding: 20px;
                background-color: #f0f2f5;
                line-height: 1.6;
            }}
            .container {{
                background-color: white;
                padding: 40px;
                border-radius: 16px;
                box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            }}
            .report-list {{
                display: grid;
                grid-template-columns: repeat(auto-fill, minmax(250px, 1fr));
                gap: 20px;
                padding: 20px 0;
            }}
            .report-link {{
                background: white;
                padding: 20px;
                border-radius: 12px;
                text-decoration: none;
                color: #1a2b3c;
                border: 1px solid #e5e7eb;
                transition: all 0.3s ease;
                display: flex;
                align-items: center;
                gap: 10px;
            }}
            .report-link:hover {{
                transform: translateY(-2px);
                box-shadow: 0 4px 15px rgba(0,0,0,0.1);
                border-color: #f97316;
            }}
            .report-icon {{
                font-size: 24px;
            }}
            .report-header {{
                text-align: center;
                margin-bottom: 40px;
                padding-bottom: 30px;
                border-bottom: 2px solid rgba(249, 115, 22, 0.2);
            }}
            .index-title {{
                font-size: 2.5rem;
                font-weight: 800;
                color: #1a2b3c;
                margin: 0;
                padding: 0;
                background: linear-gradient(135deg, #f97316, #c2410c);
                -webkit-background-clip: text;
                -webkit-text-fill-color: transparent;
                letter-spacing: -0.5px;
            }}
            .index-subtitle {{
                font-size: 1.1rem;
                color: #64748b;
                margin-top: 8px;
                font-weight: 500;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <div class="report-header">
                <h1 class="index-title">CASSIA Reports Summary</h1>
                <p class="index-subtitle">Select a report to view detailed analysis</p>
            </div>
            <div class="report-list">
                {0}
            </div>
        </div>
    </body>
    </html>
    """

    # Generate links for each report
    links = []
    for filename in sorted(report_files):
        display_name = filename.replace('report_', '').replace('.html', '')
        links.append(f'<a href="{filename}" class="report-link"><span class="report-icon">📊</span>{display_name}</a>')

    return index_template.format('\n'.join(links))


def runCASSIA_generate_score_report(csv_path, index_name="CASSIA_reports_summary"):
    """
    Generate HTML reports from a scored CSV file and create an index page.

    Args:
        csv_path (str): Path to the CSV file containing the score results
        index_name (str): Base name for the index file (without .html extension)
    """
    # Read the CSV file
    report = pd.read_csv(csv_path)
    report_files = []

    # Determine output folder (same folder as the CSV file)
    output_folder = os.path.dirname(csv_path)
    if not output_folder:
        output_folder = "."

    # Process each row
    for index, row in report.iterrows():
        # Get the first column value for the filename
        filename = str(row.iloc[0]).strip()
        filename = "".join(c for c in filename if c.isalnum() or c in (' ', '-', '_')).strip()

        # Handle both 'Conversation History' and 'Conversation.History' column names
        history_column_options = ['Conversation History', 'Conversation.History', 'conversation_history', 'Conversation_History']
        text = None
        for col in history_column_options:
            if col in row:
                text = row[col]
                break
        if text is None:
            raise KeyError(f"Could not find conversation history column. Available columns: {list(row.index)}")

        # Handle both 'Scoring_Reasoning' and 'Scoring.Reasoning' column names
        reasoning_column_options = ['Scoring_Reasoning', 'Scoring.Reasoning', 'scoring_reasoning', 'Scoring_reasoning']
        score_reasoning = None
        for col in reasoning_column_options:
            if col in row:
                score_reasoning = row[col]
                break
        if score_reasoning is None:
            raise KeyError(f"Could not find scoring reasoning column. Available columns: {list(row.index)}")

        score = row["Score"]

        # Generate HTML for this row
        html_content = process_single_report(text, score_reasoning, score)

        # Save using the first column value as filename in the output folder
        output_path = os.path.join(output_folder, f"report_{filename}.html")
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(html_content)

        # Store just the filename for the index (not the full path)
        report_files.append(os.path.basename(output_path))
        print(f"Report saved to {output_path}")

    # Generate and save index page in the same folder
    index_html = generate_score_index_page(report_files)
    index_filename = os.path.join(output_folder, f"{os.path.basename(index_name)}.html")
    with open(index_filename, "w", encoding="utf-8") as f:
        f.write(index_html)
    print(f"Index page saved to {index_filename}")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate HTML reports for cell type evaluation results')
    parser.add_argument('--dir', type=str, required=True, 
                        help='Root directory containing evaluation CSV files')
    parser.add_argument('--overwrite', action='store_true', 
                        help='Overwrite existing HTML reports')
    parser.add_argument('--index', action='store_true',
                        help='Generate an index.html file with links to all reports')
    
    args = parser.parse_args()
    
    # Find all evaluation CSV files
    csv_files = find_evaluation_csvs(args.dir)
    
    if not csv_files:
        print(f"No evaluation CSV files found in {args.dir}")
        return
    
    print(f"Found {len(csv_files)} evaluation CSV files")
    
    # Process each CSV file
    for csv_file in csv_files:
        process_evaluation_csv(csv_file, args.overwrite)
    
    # Create index.html if requested
    if args.index:
        create_index_html(csv_files, args.dir)

if __name__ == "__main__":
    main() 