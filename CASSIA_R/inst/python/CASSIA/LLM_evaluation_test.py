import os
import pandas as pd
from LLM_evaluation import (
    LLMEvaluator,
    generate_simulated_data,
    generate_multiple_celltype_samples,
    calculate_evaluation_metrics
)
import matplotlib.pyplot as plt
import seaborn as sns
import concurrent.futures
import time
import base64
from io import BytesIO
import io

def test_simulated_data():
    print("\n--- Simulated Data Generation Test ---")
    df = generate_simulated_data(5)
    print(df)
    print("\nSimulated data generated successfully.")

def test_single_evaluation():
    print("\n--- Single Celltype Evaluation Test ---")
    api_key = os.environ.get("OPENROUTER_API_KEY")
    if not api_key:
        print("[SKIP] No OpenRouter API key found. Skipping real LLM call.")
        return
    evaluator = LLMEvaluator(api_key=api_key)
    # Use a realistic example
    result = evaluator.evaluate_single_celltype(
        predicted_celltype="CD8+ T cell",
        gold_standard="CD8+ cytotoxic T cell",
        tissue="blood",
        species="human"
    )
    print("LLM Score:", result.get('score'))
    print("Explanation:", result.get('explanation'))

def test_batch_evaluation():
    print("\n--- Batch Evaluation Test ---")
    api_key = os.environ.get("OPENROUTER_API_KEY")
    if not api_key:
        print("[SKIP] No OpenRouter API key found. Skipping real LLM call.")
        return
    evaluator = LLMEvaluator(api_key=api_key)
    df = generate_simulated_data(3)
    print("Input Data:")
    print(df)
    result_df = evaluator.batch_evaluate_from_dataframe(
        df=df,
        predicted_col="predicted_celltype",
        gold_col="gold_standard",
        tissue_col="tissue",
        species_col="species"
    )
    print("\nBatch Evaluation Results:")
    print(result_df[["gold_standard", "predicted_celltype", "evaluation_score", "evaluation_explanation"]])

def test_metrics():
    print("\n--- Metrics Calculation Test ---")
    df = generate_simulated_data(20)
    metrics = calculate_evaluation_metrics(df, score_col="true_accuracy")
    print("Metrics on Simulated Data:")
    for k, v in metrics.items():
        print(f"{k}: {v}")

def visualize_results(result_df, score_col="evaluation_score"):
    print("\n--- Visualization of Evaluation Results ---")
    if score_col not in result_df.columns:
        print(f"[ERROR] Column '{score_col}' not found in results.")
        return
        
    # Get score range info to determine visualization approach
    min_score = result_df[score_col].min()
    max_score = result_df[score_col].max()
    
    # Determine if we're using rule-based (0-5) or similarity (0-100) scoring
    is_similarity_scale = max_score > 10  # Simple heuristic to detect similarity scale
    
    plt.figure(figsize=(10, 5))
    
    if is_similarity_scale:
        # For continuous 0-100 similarity scale
        bins = [0, 20, 40, 60, 80, 100]
        bin_labels = ['0-20', '20-40', '40-60', '60-80', '80-100']
        
        # Create histogram with continuous bins
        ax = sns.histplot(result_df[score_col], bins=bins, kde=False)
        plt.title("Distribution of Similarity Scores")
        plt.xlabel("Similarity Score")
        plt.ylabel("Count")
        plt.xticks([10, 30, 50, 70, 90], bin_labels)
        
        # Add count labels on bars
        for i, patch in enumerate(ax.patches):
            ax.annotate(f'{int(patch.get_height())}', 
                        (patch.get_x() + patch.get_width() / 2, patch.get_height()), 
                        ha='center', va='bottom')
    else:
        # For discrete 0-5 scale
        sns.histplot(result_df[score_col], bins=[-0.5,0.5,1.5,2.5,3.5,4.5,5.5], kde=False, discrete=True)
        plt.title("Distribution of Evaluation Scores")
        plt.xlabel("Score")
        plt.ylabel("Count")
        plt.xticks([0,1,2,3,4,5])
    
    plt.grid(axis='y')
    plt.tight_layout()
    plt.savefig("score_distribution_hist.png")
    plt.close()
    
    # Bar plot for proportions - adjust for similarity scores
    plt.figure(figsize=(8, 4))
    
    if is_similarity_scale:
        # For similarity scale, create bins and count
        bins = [0, 20, 40, 60, 80, 100]
        bin_labels = ['0-20', '20-40', '40-60', '60-80', '80-100']
        binned_scores = pd.cut(result_df[score_col], bins=bins, labels=bin_labels, include_lowest=True)
        score_counts = binned_scores.value_counts().sort_index()
        
        # Use modern seaborn syntax to avoid warning
        sns.barplot(x=score_counts.index, y=score_counts.values, hue=score_counts.index, legend=False)
    else:
        # For rule-based, use value counts directly
        score_counts = result_df[score_col].value_counts().sort_index()
        # Use modern seaborn syntax to avoid warning
        sns.barplot(x=score_counts.index, y=score_counts.values, hue=score_counts.index, legend=False)
        
    plt.title("Score Distribution (Bar Plot)")
    plt.xlabel("Score Range" if is_similarity_scale else "Score")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig("score_distribution_bar.png")
    plt.close()

def generate_html_report(result_df, gold_col, pred_col, score_col="evaluation_score", metrics=None, html_report_path="report.html"):
    """Generate an HTML report for the evaluation results."""
    from datetime import datetime
    
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
        plt.title("Distribution of Similarity Scores", fontsize=22, fontweight='bold')
        plt.xlabel("Similarity Score", fontsize=18)
        plt.ylabel("Count", fontsize=18)
        plt.xticks([10, 30, 50, 70, 90], bin_labels, fontsize=16)
    else:
        # For discrete 0-5 scale
        ax = sns.histplot(result_df[score_col], bins=[-0.5,0.5,1.5,2.5,3.5,4.5,5.5], kde=False, discrete=True, color=palette[3], edgecolor='black')
        plt.title("Distribution of Evaluation Scores", fontsize=22, fontweight='bold')
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

    # Metrics summary
    if metrics is None:
        metrics = calculate_evaluation_metrics(result_df, score_col=score_col)
    
    # Calculate score ratio - adjust for similarity scale if needed
    score_sum = result_df[score_col].sum()
    score_count = result_df[score_col].count()
    max_possible = 100 if is_similarity_scale else 5
    score_ratio = (score_sum / (score_count * max_possible)) if score_count > 0 else 0
    # Only show detailed ratios for discrete scale
    if is_similarity_scale:
        metrics_html = "<ul>" + "".join([
            f"<li><b>{k.replace('_',' ').capitalize()}:</b> {v:.3f}</li>" for k, v in metrics.items() if k in ['mean_score','median_score','min_score','max_score','std_score','count']
        ]) + f"<li><b>Score Ratio:</b> {score_ratio*100:.1f}%</li></ul>"
    else:
        metrics_html = "<ul>" + "".join([
            f"<li><b>{k.replace('_',' ').capitalize()}:</b> {v:.3f}</li>" for k, v in metrics.items()
        ]) + f"<li><b>Score Ratio:</b> {score_ratio*100:.1f}%</li></ul>"

    # Sample results: handle different for similarity vs discrete scale
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
                    reason_col = 'reasoning' if 'reasoning' in row else ('similarity_reasoning' if 'similarity_reasoning' in row else 'evaluation_explanation')
                    expl = row[reason_col] if reason_col in row else ''
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
                    reason_col = 'reasoning' if 'reasoning' in row else ('similarity_reasoning' if 'similarity_reasoning' in row else 'evaluation_explanation')
                    expl = row[reason_col] if reason_col in row else ''
                    sample_rows.append(f'<tr><td>{gold}</td><td>{pred}</td><td>{scr}</td><td>{expl[:200]}...</td></tr>')

    # HTML content
    html = f"""
    <html>
    <head>
        <title>LLM Celltype Annotation Evaluation Report</title>
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
        <h1>LLM Celltype Annotation Evaluation Report</h1>
        <div class="section">
            <b>Generated:</b> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}<br>
            <b>Total Samples:</b> {len(result_df)}<br>
            <b>Score Type:</b> {"Similarity (0-100)" if is_similarity_scale else "Discrete (0-5)"}<br>
        </div>
        <div class="section metrics">
            <h2>Summary Metrics</h2>
            {metrics_html}
        </div>
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

def test_external_dataset(file_path, gold_col, pred_col, tissue_col=None, species_col=None, save_path=None, visualize=False, n=1, max_workers=4, html_report_path=None, retry_zero_score=1, model="deepseek/deepseek-chat-v3-0324", cassia_format=False, scoring_mode="rule", max_rows=None):
    """
    Evaluate an external dataset using the LLM evaluator.
    Supports different scoring modes: 'rule' (0-5 scale) or 'similarity' (0-100 scale).
    After the main evaluation, any rows with score 0 will be retried up to 'retry_zero_score' times (default 1). If a retry yields a score > 0, the score and explanation are updated.
    The LLM model can be selected with the 'model' parameter.
    If cassia_format is True, automatically create a predicted celltype column by splitting the pred_col entries on commas and use that as the prediction column.
    Args:
        scoring_mode (str): Scoring mode to use: 'rule' or 'similarity'. Defaults to 'rule'.
        max_rows (int, optional): Maximum number of rows to process. If None, process all rows.
    """
    print(f"\n--- External Dataset Batch Evaluation (Mode: {scoring_mode}) ---")
    if not os.path.isfile(file_path):
        print(f"[ERROR] File not found: {file_path}")
        return

    api_key = os.environ.get("OPENROUTER_API_KEY")
    if not api_key:
        print("[SKIP] No OpenRouter API key found. Skipping real LLM call.")
        return

    if file_path.endswith('.xlsx'):
        df = pd.read_excel(file_path)
    else:
        df = pd.read_csv(file_path)
        
    # Limit rows if max_rows is specified
    if max_rows is not None and max_rows > 0:
        print(f"[INFO] Limiting to first {max_rows} rows (from {len(df)} total)")
        df = df.head(max_rows)

    # If cassia_format, create a new prediction column from pred_col
    if cassia_format:
        new_pred_col = pred_col + "_processed"
        df[new_pred_col] = df[pred_col].astype(str).apply(lambda x: x.split(",")[0].strip() if "," in x else x.strip())
        pred_col = new_pred_col
        print(f"[INFO] Cassia format detected. Using '{pred_col}' as the prediction column.")

    evaluator = LLMEvaluator(api_key=api_key, model=model)

    results = []
    start_time = time.time()

    # Define unified column names
    if scoring_mode == "similarity":
        score_col_name = "score"
        reasoning_col_name = "reasoning"
    else: # Default to rule-based
        score_col_name = "score"
        reasoning_col_name = "reasoning"

    if n == 1:
        def eval_row(row):
            tissue = row[tissue_col] if tissue_col else "unknown"
            species = row[species_col] if species_col else "human"
            
            if scoring_mode == "similarity":
                llm_result = evaluator.evaluate_single_similarity(
                    predicted_celltype=row[pred_col],
                    gold_standard=row[gold_col],
                    tissue=tissue,
                    species=species
                )
                score = llm_result.get("similarity_score")
                reasoning = llm_result.get("similarity_reasoning")
                print(f"[DEBUG] Single similarity result: score={score}, reasoning={reasoning[:50]}...")

                return {
                    gold_col: row[gold_col],
                    pred_col: row[pred_col],
                    score_col_name: score,
                    reasoning_col_name: reasoning,
                    **({tissue_col: row[tissue_col]} if tissue_col else {}),
                    **({species_col: row[species_col]} if species_col else {})
                }
            else: # Default to rule-based
                llm_result = evaluator.evaluate_single_celltype(
                    predicted_celltype=row[pred_col],
                    gold_standard=row[gold_col],
                    tissue=tissue,
                    species=species
                )
                score = llm_result.get("score")
                reasoning = llm_result.get("explanation")

            return {
                gold_col: row[gold_col],
                pred_col: row[pred_col],
                score_col_name: score,
                reasoning_col_name: reasoning,
                **({tissue_col: row[tissue_col]} if tissue_col else {}),
                **({species_col: row[species_col]} if species_col else {})
            }
        rows = [row for _, row in df.iterrows()]
        total = len(rows)
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_idx = {executor.submit(eval_row, row): idx for idx, row in enumerate(rows)}
            results = [None] * total
            for i, future in enumerate(concurrent.futures.as_completed(future_to_idx)):
                idx = future_to_idx[future]
                try:
                    results[idx] = future.result()
                except Exception as exc:
                    print(f"[ERROR] Row {idx} generated an exception: {exc}")
                if (i+1) % 10 == 0 or (i+1) == total:
                    print(f"Processed {i+1}/{total} rows...")
        result_df = pd.DataFrame(results)
    else:
        # n > 1: use the multiple prompt for each batch
        def eval_batch(batch_df, batch_idx):
            print(f"[INFO] Processing batch {batch_idx+1} ({len(batch_df)} rows) using {scoring_mode} scoring...")
            batch_start = time.time()
            predicted_celltypes = batch_df[pred_col].tolist()
            gold_standards = batch_df[gold_col].tolist()
            tissues = batch_df[tissue_col].tolist() if tissue_col else ["unknown"] * len(batch_df)
            species_list = batch_df[species_col].tolist() if species_col else ["human"] * len(batch_df)
            # Use the first tissue/species in the batch for context
            tissue = tissues[0] if tissues else "unknown"
            specie = species_list[0] if species_list else "human"
            
            scores = []
            explanations = []
            if scoring_mode == "similarity":
                eval_result = evaluator.evaluate_multiple_similarity(
                    predicted_celltypes=predicted_celltypes,
                    gold_standards=gold_standards,
                    tissue=tissue,
                    species=specie
                )
                scores = eval_result.get('similarity_scores', [])
                explanations = eval_result.get('similarity_reasonings', [])
            else: # Default to rule-based
                eval_result = evaluator.evaluate_multiple_celltypes(
                    predicted_celltypes=predicted_celltypes,
                    gold_standards=gold_standards,
                    tissue=tissue,
                    species=specie
                )
                scores = eval_result.get('individual_scores', [])
                explanations = eval_result.get('explanations', [])

            # Expand results for each row in the batch
            batch_results_list = []
            for i, (row_idx, row) in enumerate(batch_df.iterrows()):
                batch_results_list.append({
                    gold_col: row[gold_col],
                    pred_col: row[pred_col],
                    score_col_name: scores[i] if i < len(scores) else None,
                    reasoning_col_name: explanations[i] if i < len(explanations) else None,
                    **({tissue_col: row[tissue_col]} if tissue_col else {}),
                    **({species_col: row[species_col]} if species_col else {})
                })
            print(f"[INFO] Finished batch {batch_idx+1} in {time.time()-batch_start:.1f}s.")
            return batch_results_list
            
        batches = [df.iloc[i:i+n] for i in range(0, len(df), n)]
        total_batches = len(batches)
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(eval_batch, batch, idx) for idx, batch in enumerate(batches)]
            batch_results = []
            for i, future in enumerate(concurrent.futures.as_completed(futures)):
                try:
                    batch_results.extend(future.result())
                except Exception as exc:
                    print(f"[ERROR] Batch {i+1} generated an exception: {exc}")
                print(f"[PROGRESS] {i+1}/{total_batches} batches completed.")
        result_df = pd.DataFrame(batch_results)

    print(f"\nBatch Evaluation Results (Total time: {time.time()-start_time:.1f}s):")
    # Use unified column names for printing
    print(result_df[[gold_col, pred_col, score_col_name, reasoning_col_name]])

    # Retry rows with score 0 up to retry_zero_score times (using the selected mode)
    for retry in range(retry_zero_score):
        # Print some samples for debugging
        print(f"[DEBUG] DataFrame before retry (first 5 rows):\n{result_df.head(5)}")
        
        # Check for None scores as well, especially for similarity mode
        zero_mask = (result_df[score_col_name] == 0) | (result_df[score_col_name].isnull())
        print(f"[DEBUG] Found {zero_mask.sum()} rows with score 0 or None.")
        print(f"[DEBUG] Score column stats: min={result_df[score_col_name].min()}, max={result_df[score_col_name].max()}, null={result_df[score_col_name].isnull().sum()}")
        
        zero_rows = result_df[zero_mask]
        if zero_rows.empty:
            break
        print(f"[RETRY] Attempt {retry+1}: Retrying {len(zero_rows)} rows with score 0 or None...")
        updated = 0

        def retry_eval_row(idx_row_tuple):
            idx, row = idx_row_tuple
            tissue = row[tissue_col] if tissue_col else "unknown"
            species = row[species_col] if species_col else "human"
            score = None
            reasoning = None
            if scoring_mode == "similarity":
                eval_result = evaluator.evaluate_single_similarity(
                    predicted_celltype=row[pred_col],
                    gold_standard=row[gold_col],
                    tissue=tissue,
                    species=species
                )
                score = eval_result.get('similarity_score')
                reasoning = eval_result.get('similarity_reasoning')
            else:
                eval_result = evaluator.evaluate_single_celltype(
                    predicted_celltype=row[pred_col],
                    gold_standard=row[gold_col],
                    tissue=tissue,
                    species=species
                )
                score = eval_result.get('score')
                reasoning = eval_result.get('explanation')
            return idx, score, reasoning

        # Use ThreadPoolExecutor for parallel retries
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_idx = {executor.submit(retry_eval_row, (idx, row)): idx for idx, row in zero_rows.iterrows()}
            for future in concurrent.futures.as_completed(future_to_idx):
                idx = future_to_idx[future]
                try:
                    idx, score, reasoning = future.result()
                    if score is not None and score > 0:
                        result_df.at[idx, score_col_name] = score
                        result_df.at[idx, reasoning_col_name] = reasoning
                        updated += 1
                except Exception as exc:
                    print(f"[RETRY][ERROR] Row {idx} generated an exception: {exc}")
        print(f"[RETRY] Updated {updated} out of {len(zero_rows)} rows.")

    if save_path:
        result_df.to_csv(save_path, index=False)
        print(f"Results saved to {save_path}")
    if visualize:
        # Use unified score column name for visualization
        visualize_results(result_df, score_col=score_col_name)
    if html_report_path:
        # Use unified column names for report generation
        metrics = calculate_evaluation_metrics(result_df, score_col=score_col_name)
        generate_html_report(
            result_df=result_df, 
            gold_col=gold_col, 
            pred_col=pred_col, 
            score_col=score_col_name, 
            metrics=metrics, 
            html_report_path=html_report_path
        )


def main():
    #test_simulated_data()
    #test_single_evaluation()
    #test_batch_evaluation()
    #test_metrics()
    # Example usage in main (commented out):
    test_external_dataset(
        file_path="C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/allresults/QWEN3-235b/combined_csv_QWEN3-235b.csv",
        gold_col="True Cell Type",
        pred_col="Predicted Sub Cell Types",
        tissue_col="Tissue",
        species_col="Species",
        save_path="C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/allresults/QWEN3-235b/combined_scores_QWEN3-235b_evaluation.csv",
        visualize=True,
        n=1,                # or n=5, n=10, etc.
        max_workers=10,     # number of parallel workers
        html_report_path="C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/allresults/QWEN3-235b/combined_scores_QWEN3-235b_evaluation.html",
        model="deepseek/deepseek-chat-v3-0324",
        cassia_format=True,
        scoring_mode="similarity",
        max_rows=100        # Limit to 5 rows for testing
    )

if __name__ == "__main__":
    main()
