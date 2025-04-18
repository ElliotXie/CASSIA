import os
import pandas as pd
from LLM_evalution import (
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


os.environ["OPENROUTER_API_KEY"] = "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"
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
    plt.figure(figsize=(10, 5))
    sns.histplot(result_df[score_col], bins=[-0.5,0.5,1.5,2.5,3.5,4.5,5.5], kde=False, discrete=True)
    plt.title("Distribution of Evaluation Scores")
    plt.xlabel("Score")
    plt.ylabel("Count")
    plt.xticks([0,1,2,3,4,5])
    plt.grid(axis='y')
    plt.close()
    # Bar plot for proportions
    score_counts = result_df[score_col].value_counts().sort_index()
    plt.figure(figsize=(8, 4))
    sns.barplot(x=score_counts.index, y=score_counts.values, palette="viridis")
    plt.title("Score Distribution (Bar Plot)")
    plt.xlabel("Score")
    plt.ylabel("Count")
    plt.close()

def generate_html_report(result_df, gold_col, pred_col, score_col="evaluation_score", metrics=None, html_report_path="report.html"):
    from datetime import datetime
    import base64
    from io import BytesIO
    # Generate a larger, more beautiful histogram and encode as base64 (do not show on screen)
    buf1 = BytesIO()
    plt.figure(figsize=(14, 7))
    sns.set_theme(style="whitegrid")
    palette = sns.color_palette("crest", 6)
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
    # Calculate score ratio
    score_sum = result_df[score_col].sum()
    score_count = result_df[score_col].count()
    score_ratio = (score_sum / (score_count * 5)) if score_count > 0 else 0
    metrics_html = "<ul>" + "".join([f"<li><b>{k}:</b> {v:.3f}</li>" for k, v in metrics.items()]) + f"<li><b>Score Ratio:</b> {score_ratio*100:.1f}%</li></ul>"

    # Sample results: up to 5 for each score (0-5)
    sample_rows = []
    for score in range(6):
        score_df = result_df[result_df[score_col] == score]
        n = min(5, len(score_df))
        if n > 0:
            for _, row in score_df.head(n).iterrows():
                gold = row[gold_col] if gold_col in row else ''
                pred = row[pred_col] if pred_col in row else ''
                scr = row[score_col] if score_col in row else ''
                expl = row['evaluation_explanation'] if 'evaluation_explanation' in row else ''
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
            <p><i>Showing up to 5 examples for each score (0-5).</i></p>
        </div>
    </body>
    </html>
    """
    with open(html_report_path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"HTML report saved to {html_report_path}")

def test_external_dataset(file_path, gold_col, pred_col, tissue_col=None, species_col=None, save_path=None, visualize=False, n=1, max_workers=4, html_report_path=None, retry_zero_score=1, model="deepseek/deepseek-chat-v3-0324"):
    """
    Evaluate an external dataset using the LLM evaluator. After the main evaluation, any rows with score 0 will be retried up to 'retry_zero_score' times (default 1). If a retry yields a score > 0, the score and explanation are updated.
    The LLM model can be selected with the 'model' parameter.
    """
    print("\n--- External Dataset Batch Evaluation ---")
    if not os.path.isfile(file_path):
        print(f"[ERROR] File not found: {file_path}")
        return

    api_key = os.environ.get("OPENROUTER_API_KEY")
    if not api_key:
        print("[SKIP] No OpenRouter API key found. Skipping real LLM call.")
        return

    df = pd.read_csv(file_path)
    evaluator = LLMEvaluator(api_key=api_key, model=model)

    results = []
    start_time = time.time()
    if n == 1:
        def eval_row(row):
            return {
                gold_col: row[gold_col],
                pred_col: row[pred_col],
                "evaluation_score": None,
                "evaluation_explanation": None,
                **({tissue_col: row[tissue_col]} if tissue_col else {}),
                **({species_col: row[species_col]} if species_col else {}),
                **evaluator.evaluate_single_celltype(
                    predicted_celltype=row[pred_col],
                    gold_standard=row[gold_col],
                    tissue=row[tissue_col] if tissue_col else "unknown",
                    species=row[species_col] if species_col else "human"
                )
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
        def eval_batch(batch_df, batch_idx):
            print(f"[INFO] Processing batch {batch_idx+1} ({len(batch_df)} rows)...")
            batch_start = time.time()
            result = evaluator.batch_evaluate_from_dataframe(
                df=batch_df,
                predicted_col=pred_col,
                gold_col=gold_col,
                tissue_col=tissue_col,
                species_col=species_col
            )
            print(f"[INFO] Finished batch {batch_idx+1} in {time.time()-batch_start:.1f}s.")
            return result
        batches = [df.iloc[i:i+n] for i in range(0, len(df), n)]
        total_batches = len(batches)
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(eval_batch, batch, idx) for idx, batch in enumerate(batches)]
            batch_results = []
            for i, future in enumerate(concurrent.futures.as_completed(futures)):
                try:
                    batch_results.append(future.result())
                except Exception as exc:
                    print(f"[ERROR] Batch {i+1} generated an exception: {exc}")
                print(f"[PROGRESS] {i+1}/{total_batches} batches completed.")
        result_df = pd.concat(batch_results, ignore_index=True)

    print(f"\nBatch Evaluation Results (Total time: {time.time()-start_time:.1f}s):")
    print(result_df[[gold_col, pred_col, "evaluation_score", "evaluation_explanation"]])

    # Retry rows with score 0 up to retry_zero_score times
    for retry in range(retry_zero_score):
        zero_mask = result_df["evaluation_score"] == 0
        zero_rows = result_df[zero_mask]
        if zero_rows.empty:
            break
        print(f"[RETRY] Attempt {retry+1}: Retrying {len(zero_rows)} rows with score 0...")
        updated = 0
        for idx, row in zero_rows.iterrows():
            eval_result = evaluator.evaluate_single_celltype(
                predicted_celltype=row[pred_col],
                gold_standard=row[gold_col],
                tissue=row[tissue_col] if tissue_col else "unknown",
                species=row[species_col] if species_col else "human"
            )
            if eval_result.get('score', 0) > 0:
                result_df.at[idx, 'evaluation_score'] = eval_result.get('score', 0)
                result_df.at[idx, 'evaluation_explanation'] = eval_result.get('explanation', '')
                updated += 1
        print(f"[RETRY] Updated {updated} out of {len(zero_rows)} rows.")

    if save_path:
        result_df.to_csv(save_path, index=False)
        print(f"Results saved to {save_path}")
    if visualize:
        visualize_results(result_df, score_col="evaluation_score")
    if html_report_path:
        metrics = calculate_evaluation_metrics(result_df, score_col="evaluation_score")
        generate_html_report(result_df, gold_col, pred_col, score_col="evaluation_score", metrics=metrics, html_report_path=html_report_path)

def main():
    #test_simulated_data()
    #test_single_evaluation()
    #test_batch_evaluation()
    #test_metrics()
    # Example usage in main (commented out):
    test_external_dataset(
        file_path="C:/Users/ellio/Downloads/mode_evaluation.xlsx - Llama3.2.csv",
        gold_col="True Cell Type",
        pred_col="Predicted",
        tissue_col="Tissue",
        species_col="Species",
        save_path="results_tested14.csv",
        visualize=True,
        n=3,                # or n=5, n=10, etc.
        max_workers=6,       # number of parallel workers
        html_report_path="report14.html",
        model="deepseek/deepseek-chat-v3-0324"
    )

if __name__ == "__main__":
    main()
