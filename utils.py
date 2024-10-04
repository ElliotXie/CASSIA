import json
import csv

def create_csv_from_json(json_filename):
    # Extract tissue name and number of genes from the JSON filename
    filename_parts = json_filename.split('_')
    tissue_name = filename_parts[0]
    num_genes = next(part for part in filename_parts if part.endswith('genes')).split('genes')[0]

    # Read the JSON file
    with open(json_filename, 'r') as json_file:
        data = json.load(json_file)

    # Function to write CSV files
    def write_csv(filename, headers, row_data):
        with open(filename, 'w', newline='', encoding='utf-8') as csv_file:
            writer = csv.writer(csv_file)
            writer.writerow(headers)
            writer.writerows(row_data)

    # Prepare data for both CSV files
    full_data = []
    summary_data = []

    for true_cell_type, details in data.items():
        main_cell_type = details['analysis_result']['main_cell_type']
        sub_cell_types = ', '.join(details['analysis_result']['sub_cell_types'])
        marker_number = details['analysis_result']['num_markers']
        conversation_history = ' | '.join([f"{entry[0]}: {entry[1]}" for entry in details['conversation_history']])
        
        full_data.append([true_cell_type, main_cell_type, sub_cell_types, 
                          marker_number, conversation_history])
        summary_data.append([true_cell_type, main_cell_type, sub_cell_types])

    # Write the full data CSV
    full_csv_filename = f'{tissue_name}_cell_type_analysis_results_full_{num_genes}genes.csv'
    write_csv(full_csv_filename, 
              ['True Cell Type', 'Predicted Main Cell Type', 'Predicted Sub Cell Types', 
               'Marker Number', 'Conversation History'],
              full_data)

    # Write the summary data CSV
    summary_csv_filename = f'{tissue_name}_cell_type_analysis_results_summary_{num_genes}genes.csv'
    write_csv(summary_csv_filename,
              ['True Cell Type', 'Predicted Main Cell Type', 'Predicted Sub Cell Types'],
              summary_data)

    print("Two CSV files have been created:")
    print(f"1. {full_csv_filename} (full data)")
    print(f"2. {summary_csv_filename} (summary data)")