import json
import csv

# Read the JSON file
with open('cell_type_analysis_results_10.json', 'r') as json_file:
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
    confidence_score = details['analysis_result']['confidence_score']
    conversation_history = ' | '.join([f"{entry[0]}: {entry[1]}" for entry in details['conversation_history']])
    
    full_data.append([true_cell_type, main_cell_type, sub_cell_types, 
                      marker_number, confidence_score, conversation_history])
    summary_data.append([true_cell_type, main_cell_type, sub_cell_types])

# Write the full data CSV
write_csv('cell_type_analysis_results_full10.csv', 
          ['True Cell Type', 'Predicted Main Cell Type', 'Predicted Sub Cell Types', 
           'Marker Number', 'Confidence Score', 'Conversation History'],
          full_data)

# Write the summary data CSV
write_csv('cell_type_analysis_results_summary10.csv',
          ['True Cell Type', 'Predicted Main Cell Type', 'Predicted Sub Cell Types'],
          summary_data)

print("Two CSV files have been created:")
print("1. cell_type_analysis_results_full.csv (full data)")
print("2. cell_type_analysis_results_summary.csv (summary data)")