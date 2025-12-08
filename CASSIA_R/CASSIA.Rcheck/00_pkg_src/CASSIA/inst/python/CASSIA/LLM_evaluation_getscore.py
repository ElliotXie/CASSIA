import pandas as pd

result_df = pd.read_csv("C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/o4mini/combined_scores_o4mini_evaluation.csv")

    #group based on the Tissue	Species and then calcualte the average score in the score column，save that stat to a new csv file

result_df_stat= result_df.groupby(["Tissue", "Species"])["score"].mean().reset_index()

# Save a title row to the csv file, using the model name
model_name = "O4Mini"
output_path = "C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/o4mini/combined_scores_o4mini_evaluation_average.csv"

with open(output_path, 'w', encoding='utf-8') as f:
    f.write(f"Model: {model_name}\n")
    result_df_stat.to_csv(f, index=False)
    

def process_multiple_files(file_model_pairs, output_txt_path):
    """
    Processes a list of (CSV file path, model name) pairs, computes grouped averages,
    and writes all results to a TXT file. Each section shows the model name, file name,
    grouped average scores, and overall average score for that model.
    """
    with open(output_txt_path, 'w', encoding='utf-8') as out_f:
        for file_path, model_name in file_model_pairs:
            try:
                # Read the CSV file
                df = pd.read_csv(file_path)
                
                # Calculate the overall average for this model (from the full CSV)
                overall_model_avg = df["score"].mean()
                
                # Calculate the grouped averages
                result_df_stat = df.groupby(["Tissue", "Species"])["score"].mean().reset_index()
                
                # Write model and file info
                out_f.write(f"Model: {model_name}\n")
                out_f.write(f"File: {file_path}\n")
                
                # Write the overall model average
                out_f.write(f"Overall average score for {model_name}: {overall_model_avg:.4f}\n\n")
                
                # Write the grouped statistics
                out_f.write("Grouped statistics:\n")
                out_f.write(result_df_stat.to_csv(index=False))
                out_f.write("\n\n")
                
            except Exception as e:
                out_f.write(f"Model: {model_name}\nFile: {file_path}\nError: {e}\n\n")

# Example usage:
file_model_list = [
    ("C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/o4mini/combined_scores_o4mini_evaluation.csv", "O4Mini"),
    ("C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/Gemini2.5pro_2/combined_scores_gemini-2.5pro2_evaluation.csv", "Gemini-2.5pro"),
    ("C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/gemini2.5flash/combined_scores_gemini2.5flash_evaluation.csv", "Gemini2.5flash"),
    ("C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/GPT4.1/combined_scores_GPT4.1_evaluation.csv", "GPT4.1"),
    ("C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/claude3.7/combined_scores_claude3.7_evaluation.csv", "claude3.7")
]
process_multiple_files(file_model_list, "C:/Users/ellio/OneDrive - UW-Madison/CASSIA+/combined_results.txt")





