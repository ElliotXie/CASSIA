import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob
from typing import List, Dict, Any, Tuple
import re

def find_evaluation_csvs(root_dir: str) -> List[str]:
    """Find all evaluation CSV files in the given directory and its subdirectories."""
    pattern = os.path.join(root_dir, "**", "*_evaluation.csv")
    return glob(pattern, recursive=True)

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

def load_and_summarize_csv(csv_path: str) -> Dict[str, Any]:
    """Load a CSV file and calculate summary metrics."""
    try:
        df = pd.read_csv(csv_path)
        
        # Determine score column
        score_col = None
        for col in ["score", "evaluation_score", "similarity_score"]:
            if col in df.columns:
                score_col = col
                break
        
        if not score_col:
            print(f"No score column found in {csv_path}. Skipping.")
            return None
        
        # Calculate basic statistics
        model_name = extract_model_name(csv_path)
        
        # Overall statistics
        stats = {
            'model': model_name,
            'file': csv_path,
            'count': len(df),
            'mean_score': df[score_col].mean(),
            'median_score': df[score_col].median(),
            'min_score': df[score_col].min(),
            'max_score': df[score_col].max(),
            'std_score': df[score_col].std()
        }
        
        # Group by Tissue and Species if those columns exist
        if 'Tissue' in df.columns and 'Species' in df.columns:
            grouped = df.groupby(['Tissue', 'Species'])[score_col].agg(['mean', 'count']).reset_index()
            grouped_stats = grouped.to_dict('records')
            stats['grouped_stats'] = grouped_stats
        
        return stats
    
    except Exception as e:
        print(f"Error processing {csv_path}: {str(e)}")
        return None

def generate_comparison_visualization(model_stats: List[Dict[str, Any]], output_dir: str) -> None:
    """Generate visualizations comparing models."""
    # Create a DataFrame from model stats
    df = pd.DataFrame(model_stats)
    
    # Sort by mean score
    df = df.sort_values('mean_score', ascending=False)
    
    # Set up the matplotlib figure
    plt.figure(figsize=(12, 8))
    
    # Create a bar chart
    ax = sns.barplot(x='model', y='mean_score', data=df, palette='viridis')
    
    # Customize the plot
    plt.title('Average Score by Model', fontsize=18)
    plt.xlabel('Model', fontsize=14)
    plt.ylabel('Average Score', fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.grid(axis='y', alpha=0.3)
    
    # Add value labels on top of bars
    for i, v in enumerate(df['mean_score']):
        ax.text(i, v + 0.5, f"{v:.2f}", ha='center', fontweight='bold')
    
    # Set y-axis limits
    plt.ylim(0, 100 if df['mean_score'].max() > 10 else 5)
    
    # Save the figure
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'model_comparison.png'), dpi=300)
    plt.close()
    
    # If we have grouped stats, create a grouped comparison
    if 'grouped_stats' in df.iloc[0] and df.iloc[0]['grouped_stats'] is not None:
        create_grouped_comparison(df, output_dir)

def create_grouped_comparison(df: pd.DataFrame, output_dir: str) -> None:
    """Create visualizations for tissue/species breakdown."""
    # Extract and combine all grouped stats
    all_groups = []
    for _, row in df.iterrows():
        model = row['model']
        if 'grouped_stats' in row and row['grouped_stats'] is not None:
            for group_stat in row['grouped_stats']:
                group_stat['model'] = model
                all_groups.append(group_stat)
    
    if not all_groups:
        return
        
    grouped_df = pd.DataFrame(all_groups)
    
    # Create a pivot table for visualization
    pivot_df = grouped_df.pivot_table(
        index=['Tissue', 'Species'], 
        columns='model', 
        values='mean'
    ).reset_index()
    
    # Create heatmap for each tissue+species combination
    tissues = pivot_df['Tissue'].unique()
    for tissue in tissues:
        tissue_df = pivot_df[pivot_df['Tissue'] == tissue]
        
        # If only one species, it's simpler
        if len(tissue_df['Species'].unique()) == 1:
            # Extract the scores part (drop Tissue and Species columns)
            scores_df = tissue_df.drop(['Tissue', 'Species'], axis=1)
            
            # Reshape for better visualization
            plot_df = scores_df.T.reset_index()
            plot_df.columns = ['model', 'score']
            
            # Sort by score
            plot_df = plot_df.sort_values('score', ascending=False)
            
            # Create bar chart
            plt.figure(figsize=(10, 6))
            ax = sns.barplot(x='model', y='score', data=plot_df, palette='viridis')
            plt.title(f'Model Comparison: {tissue} ({tissue_df["Species"].iloc[0]})', fontsize=16)
            plt.xlabel('Model', fontsize=14)
            plt.ylabel('Average Score', fontsize=14)
            plt.xticks(rotation=45, ha='right')
            
            # Add value labels
            for i, v in enumerate(plot_df['score']):
                ax.text(i, v + 0.5, f"{v:.2f}", ha='center', fontweight='bold')
            
            # Set y-axis limits
            plt.ylim(0, 100 if plot_df['score'].max() > 10 else 5)
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'comparison_{tissue.replace(" ", "_")}.png'), dpi=300)
            plt.close()
        else:
            # Multiple species - create a grouped bar chart
            # Melt the dataframe for easier plotting
            melted_df = pd.melt(
                tissue_df, 
                id_vars=['Tissue', 'Species'], 
                var_name='model', 
                value_name='score'
            )
            
            plt.figure(figsize=(12, 8))
            ax = sns.barplot(x='model', y='score', hue='Species', data=melted_df, palette='Set2')
            plt.title(f'Model Comparison: {tissue}', fontsize=16)
            plt.xlabel('Model', fontsize=14)
            plt.ylabel('Average Score', fontsize=14)
            plt.xticks(rotation=45, ha='right')
            plt.legend(title='Species')
            
            # Set y-axis limits
            plt.ylim(0, 100 if melted_df['score'].max() > 10 else 5)
            
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f'comparison_{tissue.replace(" ", "_")}_by_species.png'), dpi=300)
            plt.close()

def generate_markdown_table(model_stats: List[Dict[str, Any]], output_dir: str) -> None:
    """Generate a markdown table with model comparison results."""
    # Create a DataFrame from model stats
    df = pd.DataFrame(model_stats)
    
    # Sort by mean score
    df = df.sort_values('mean_score', ascending=False)
    
    # Create markdown content
    md_content = "# Model Comparison Results\n\n"
    md_content += f"Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
    
    # Overall results table
    md_content += "## Overall Results\n\n"
    md_content += "| Rank | Model | Average Score | Median Score | Min | Max | Std Dev | Count |\n"
    md_content += "|------|-------|--------------|--------------|-----|-----|---------|-------|\n"
    
    for i, (_, row) in enumerate(df.iterrows()):
        md_content += f"| {i+1} | {row['model']} | {row['mean_score']:.2f} | {row['median_score']:.2f} | "
        md_content += f"{row['min_score']:.2f} | {row['max_score']:.2f} | {row['std_score']:.2f} | {row['count']} |\n"
    
    # Add breakdown by tissue and species if available
    if 'grouped_stats' in df.iloc[0] and df.iloc[0]['grouped_stats'] is not None:
        md_content += "\n## Results by Tissue and Species\n\n"
        
        # Extract and combine all grouped stats
        all_groups = []
        for _, row in df.iterrows():
            model = row['model']
            if 'grouped_stats' in row and row['grouped_stats'] is not None:
                for group_stat in row['grouped_stats']:
                    group_stat['model'] = model
                    all_groups.append(group_stat)
        
        grouped_df = pd.DataFrame(all_groups)
        
        # Create pivoted tables for each tissue
        tissues = grouped_df['Tissue'].unique()
        for tissue in tissues:
            tissue_df = grouped_df[grouped_df['Tissue'] == tissue]
            species_list = tissue_df['Species'].unique()
            
            md_content += f"### {tissue}\n\n"
            
            for species in species_list:
                species_df = tissue_df[tissue_df['Species'] == species]
                
                # Sort by mean
                species_df = species_df.sort_values('mean', ascending=False)
                
                md_content += f"**{species}**:\n\n"
                md_content += "| Rank | Model | Average Score | Sample Count |\n"
                md_content += "|------|-------|--------------|---------------|\n"
                
                for i, (_, row) in enumerate(species_df.iterrows()):
                    md_content += f"| {i+1} | {row['model']} | {row['mean']:.2f} | {row['count']} |\n"
                
                md_content += "\n"
    
    # Write to file
    with open(os.path.join(output_dir, 'comparison_results.md'), 'w', encoding='utf-8') as f:
        f.write(md_content)

def generate_html_comparison(model_stats: List[Dict[str, Any]], output_dir: str) -> None:
    """Generate an HTML comparison report."""
    # Create a DataFrame from model stats
    df = pd.DataFrame(model_stats)
    
    # Sort by mean score
    df = df.sort_values('mean_score', ascending=False)
    
    # Create HTML content
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Model Comparison Results</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }}
            h1, h2, h3 {{ color: #2c3e50; }}
            table {{ border-collapse: collapse; width: 100%; margin-bottom: 30px; }}
            th, td {{ text-align: left; padding: 12px; }}
            th {{ background-color: #3498db; color: white; }}
            tr:nth-child(even) {{ background-color: #f2f2f2; }}
            tr:hover {{ background-color: #ddd; }}
            .chart-container {{ margin-top: 30px; margin-bottom: 40px; }}
        </style>
    </head>
    <body>
        <h1>Model Comparison Results</h1>
        <p>Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        
        <h2>Overall Results</h2>
        <div class="chart-container">
            <img src="model_comparison.png" alt="Model Comparison Chart" style="max-width:100%; height:auto;" />
        </div>
        
        <table>
            <tr>
                <th>Rank</th>
                <th>Model</th>
                <th>Average Score</th>
                <th>Median Score</th>
                <th>Min</th>
                <th>Max</th>
                <th>Std Dev</th>
                <th>Count</th>
            </tr>
    """
    
    for i, (_, row) in enumerate(df.iterrows()):
        html_content += f"""
            <tr>
                <td>{i+1}</td>
                <td>{row['model']}</td>
                <td>{row['mean_score']:.2f}</td>
                <td>{row['median_score']:.2f}</td>
                <td>{row['min_score']:.2f}</td>
                <td>{row['max_score']:.2f}</td>
                <td>{row['std_score']:.2f}</td>
                <td>{row['count']}</td>
            </tr>
        """
    
    html_content += """
        </table>
    """
    
    # Add breakdown by tissue and species if available
    if 'grouped_stats' in df.iloc[0] and df.iloc[0]['grouped_stats'] is not None:
        html_content += """
        <h2>Results by Tissue and Species</h2>
        """
        
        # Extract and combine all grouped stats
        all_groups = []
        for _, row in df.iterrows():
            model = row['model']
            if 'grouped_stats' in row and row['grouped_stats'] is not None:
                for group_stat in row['grouped_stats']:
                    group_stat['model'] = model
                    all_groups.append(group_stat)
        
        grouped_df = pd.DataFrame(all_groups)
        
        # Create pivoted tables for each tissue
        tissues = grouped_df['Tissue'].unique()
        for tissue in tissues:
            tissue_df = grouped_df[grouped_df['Tissue'] == tissue]
            species_list = tissue_df['Species'].unique()
            
            html_content += f"""
            <h3>{tissue}</h3>
            <div class="chart-container">
                <img src="comparison_{tissue.replace(' ', '_')}{'' if len(species_list) == 1 else '_by_species'}.png" 
                     alt="Comparison for {tissue}" style="max-width:100%; height:auto;" />
            </div>
            """
            
            for species in species_list:
                species_df = tissue_df[tissue_df['Species'] == species]
                
                # Sort by mean
                species_df = species_df.sort_values('mean', ascending=False)
                
                html_content += f"""
                <h4>{species}</h4>
                <table>
                    <tr>
                        <th>Rank</th>
                        <th>Model</th>
                        <th>Average Score</th>
                        <th>Sample Count</th>
                    </tr>
                """
                
                for i, (_, row) in enumerate(species_df.iterrows()):
                    html_content += f"""
                    <tr>
                        <td>{i+1}</td>
                        <td>{row['model']}</td>
                        <td>{row['mean']:.2f}</td>
                        <td>{row['count']}</td>
                    </tr>
                    """
                
                html_content += """
                </table>
                """
    
    html_content += """
    </body>
    </html>
    """
    
    # Write to file
    with open(os.path.join(output_dir, 'comparison_results.html'), 'w', encoding='utf-8') as f:
        f.write(html_content)

def generate_text_summary(model_stats: List[Dict[str, Any]], output_dir: str) -> None:
    """Generate a simple text summary of model scores by tissue."""
    # Check if we have grouped stats
    if not model_stats or 'grouped_stats' not in model_stats[0] or not model_stats[0]['grouped_stats']:
        print("No grouped statistics available for text summary.")
        return
    
    # Extract and combine all grouped stats
    all_groups = []
    for stats in model_stats:
        model = stats['model']
        if 'grouped_stats' in stats and stats['grouped_stats']:
            for group_stat in stats['grouped_stats']:
                group_stat['model'] = model
                all_groups.append(group_stat)
    
    if not all_groups:
        return
    
    grouped_df = pd.DataFrame(all_groups)
    
    # Create a summary text
    summary_text = "Model Performance Summary by Tissue\n"
    summary_text += "===================================\n\n"
    summary_text += f"Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n"
    
    # Overall average scores by model
    model_df = pd.DataFrame(model_stats)[['model', 'mean_score']]
    model_df = model_df.sort_values('mean_score', ascending=False)
    
    summary_text += "Overall Average Scores:\n"
    summary_text += "----------------------\n"
    for _, row in model_df.iterrows():
        summary_text += f"{row['model']}: {row['mean_score']:.2f}\n"
    
    summary_text += "\n\nScores by Tissue and Species:\n"
    summary_text += "----------------------------\n"
    
    # Group by Tissue and Species
    tissues = grouped_df['Tissue'].unique()
    
    for tissue in tissues:
        summary_text += f"\n{tissue}:\n"
        tissue_df = grouped_df[grouped_df['Tissue'] == tissue]
        species_list = tissue_df['Species'].unique()
        
        for species in species_list:
            summary_text += f"  {species}:\n"
            species_df = tissue_df[tissue_df['Species'] == species]
            
            # Sort by mean score
            species_df = species_df.sort_values('mean', ascending=False)
            
            for _, row in species_df.iterrows():
                summary_text += f"    {row['model']}: {row['mean']:.2f} ({int(row['count'])} samples)\n"
    
    # Write to file
    with open(os.path.join(output_dir, 'model_summary.txt'), 'w', encoding='utf-8') as f:
        f.write(summary_text)
    
    print(f"Text summary saved to {os.path.join(output_dir, 'model_summary.txt')}")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate comparison reports for cell type evaluation results')
    parser.add_argument('--dir', type=str, required=True, 
                        help='Root directory containing evaluation CSV files')
    parser.add_argument('--output', type=str, default=None,
                        help='Output directory for comparison results (default: same as input directory)')
    
    args = parser.parse_args()
    
    # Set output directory
    output_dir = args.output if args.output else args.dir
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all evaluation CSV files
    csv_files = find_evaluation_csvs(args.dir)
    
    if not csv_files:
        print(f"No evaluation CSV files found in {args.dir}")
        return
    
    print(f"Found {len(csv_files)} evaluation CSV files")
    
    # Process each CSV file and gather statistics
    model_stats = []
    for csv_file in csv_files:
        stats = load_and_summarize_csv(csv_file)
        if stats:
            model_stats.append(stats)
            print(f"Processed {stats['model']} - Mean score: {stats['mean_score']:.2f}")
    
    if not model_stats:
        print("No valid statistics gathered. Exiting.")
        return
        
    # Generate visualizations
    generate_comparison_visualization(model_stats, output_dir)
    
    # Generate markdown table
    generate_markdown_table(model_stats, output_dir)
    
    # Generate HTML report
    generate_html_comparison(model_stats, output_dir)
    
    # Generate text summary
    generate_text_summary(model_stats, output_dir)
    
    print(f"Comparison reports generated in {output_dir}")
    print(f"  - comparison_results.html")
    print(f"  - comparison_results.md")
    print(f"  - model_comparison.png")
    print(f"  - model_summary.txt")

if __name__ == "__main__":
    main() 