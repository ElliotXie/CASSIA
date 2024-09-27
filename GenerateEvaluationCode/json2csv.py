import json
import csv

# Read the JSON file
with open('cell_type_analysis_results_10.json', 'r') as json_file:
    data = json.load(json_file)

# Open a CSV file for writing
with open('cell_type_analysis_results_full3.csv', 'w', newline='', encoding='utf-8') as csv_file:
    writer = csv.writer(csv_file)
    
    # Write the header
    writer.writerow(['True Cell Type', 'Predicted Main Cell Type', 'Predicted Sub Cell Types', 
                     'Marker Number', 'Confidence Score', 'Conversation History'])
    
    # Iterate through the JSON data
    for true_cell_type, details in data.items():
        main_cell_type = details['analysis_result']['main_cell_type']
        sub_cell_types = ', '.join(details['analysis_result']['sub_cell_types'])
        marker_number = details['analysis_result']['num_markers']
        confidence_score = details['analysis_result']['confidence_score']
        
        # Join the conversation history into a single string
        conversation_history = ' | '.join([f"{entry[0]}: {entry[1]}" for entry in details['conversation_history']])
        
        # Write the row
        writer.writerow([true_cell_type, main_cell_type, sub_cell_types, 
                         marker_number, confidence_score, conversation_history])

print("CSV file has been created: cell_type_analysis_results_full.csv")