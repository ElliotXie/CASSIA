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

    # Build table rows with collapsible markers
    rows = []
    for idx, row in df.iterrows():
        cluster = row.get('Result ID', '')
        main_type = row.get('main_cell_type', '')
        sub_type = row.get('sub_cell_type', '')
        key_markers = row.get('key_markers', '')
        reason = row.get('reason', '')
        marker_id = f"marker_{idx+1}"
        rows.append(f'''
        <tr>
            <td class="cluster-col">{cluster}</td>
            <td class="annotation-col">
                <div class="main-type">{main_type}</div>
                <div class="sub-type">{sub_type}</div>
            </td>
            <td class="marker-col">
                <button class="marker-toggle" onclick="toggleMarkers('{marker_id}')">Show Markers</button>
                <div id="{marker_id}" class="marker-content" style="display:none;">{key_markers}</div>
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
            .marker-toggle {{ background: #e3f2fd; color: #1976d2; border: none; border-radius: 5px; padding: 3px 10px; font-size: 13px; cursor: pointer; margin-bottom: 4px; }}
            .marker-toggle:hover {{ background: #bbdefb; }}
            .marker-content {{ margin-top: 4px; font-size: 13px; color: #555; background: #f8fafc; border-radius: 5px; padding: 6px 8px; box-shadow: 0 1px 2px rgba(0,0,0,0.03); }}
            .reasoning-col {{ width: 55%; }}
            .reasoning-box {{ background: #fffde7; border-left: 5px solid #ffe082; border-radius: 7px; padding: 12px 16px; font-size: 15px; color: #444; box-shadow: 0 1px 4px rgba(255,193,7,0.07); }}
        </style>
        <script>
        function toggleMarkers(id) {{
            var el = document.getElementById(id);
            if (el.style.display === 'none') {{
                el.style.display = 'block';
            }} else {{
                el.style.display = 'none';
            }}
        }}
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
    </body>
    </html>
    '''
    with open(html_report_path, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f"Subclustering HTML report saved to {html_report_path}")

def process_evaluation_csv(csv_path: str, overwrite: bool = False) -> None:
    try:
        if not os.path.exists(csv_path):
            print(f"File not found: {csv_path}")
            return
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
        # Determine column names
        gold_col = "True Cell Type" if "True Cell Type" in df.columns else "gold_standard"
        pred_col = "Predicted Sub Cell Types" if "Predicted Sub Cell Types" in df.columns else "predicted_celltype"
        
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