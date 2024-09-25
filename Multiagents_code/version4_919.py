import json
from openai import OpenAI

client = OpenAI()

class Agent:
    def __init__(self, system="", model="gpt-4o", temperature=0):
        self.system = system
        self.chat_histories = {}
        self.model = model
        self.temperature = temperature

    def __call__(self, message, other_agent_id):
        if other_agent_id not in self.chat_histories:
            self.chat_histories[other_agent_id] = []
            if self.system:
                self.chat_histories[other_agent_id].append({"role": "system", "content": self.system})
        
        self.chat_histories[other_agent_id].append({"role": "user", "content": message})
        
        result = self.execute(other_agent_id)
        self.chat_histories[other_agent_id].append({"role": "assistant", "content": result})
        
        return result

    def execute(self, other_agent_id):
        completion = client.chat.completions.create(
            model=self.model,
            temperature=self.temperature,
            messages=self.chat_histories[other_agent_id]
        )
        return completion.choices[0].message.content

def construct_prompt(species, tissue_type, marker_list, additional_info=""):
    markers = ', '.join(marker_list)
    prompt = f"I am analyzing a single-cell {species} {tissue_type} dataset."
    if additional_info:
        prompt += f" {additional_info}."
    prompt += f" I want to identify the cell types present based on this marker list:\n{markers}"
    return prompt

# Define agents
final_annotation_agent = Agent(system="""
Objective:

You are a professional computational biologist with expertise in single-cell RNA sequencing (scRNA-seq).
                               I will provide you with a list of highly expressed markers ranked by expression intensity from high to low
                               from a cluster of cells , and your task is to identify the cell type. You must think step-by-step, providing a comprehensive and specific analysis. The audience is an expert in the field, and I will tip you $1000 if you do a good job.

Steps to Follow:

1. List the Key Functional Markers: Extract and group the key marker genes associated with function or pathway, explaining their roles. Do not repeat the input markers.
2. List the Key Cell Type Markers: Extract and group the key marker genes associated with mouse larynx cell types, explaining their roles. Do not repeat the input markers.
3. Cross-reference Known Databases: Use available scRNA-seq databases and relevant literature to cross-reference these markers. list your finding.
4. Determine the Most Probable General Cell Type: Based on the expression of these markers, infer the most likely general cell type of the cluster.
5. Identify the Top 3 Most Probable Sub Cell Types: Based on the expression of these markers, infer the top three most probable sub cell types within the general cell type. Finally, specify the most likely subtype.
6. Identify the Most Probable Sub-Sub Cell Type: Determine the most specific cell type within the previously identified subtype.
7.  Provide a Concise Summary of Your Analysis
""".strip())

formatting_agent = Agent(system="""
You are a formatting assistant for single-cell analysis results. Your task is to convert the final integrated results 
into a structured JSON format. Follow these guidelines:

1. Extract the main cell type
2. Extract the sub-cell types identified, rank them based on the probability, the most likely sub-cell type should be the first one.
2. Extract the sub-sub-cell types identified if any, rank them based on the probability, the most likely sub-sub-cell type should be the first one.

Provide the JSON output within triple backticks, like this:
'''json
{
"main_cell_type": "...",
"sub_cell_types": ["...", "..."]
"sub_sub_cell_types": ["...", "..."]
}
'''
""".strip())

def run_analysis(species, tissue_type, marker_list, additional_info=""):
    print("Starting analysis process...\n")
    prompt = construct_prompt(species, tissue_type, marker_list, additional_info)

    print("Performing final annotation...\n")
    annotation_result = final_annotation_agent(prompt, "user")
    print(annotation_result)
    print("Formatting final results...\n")
    formatted_output = formatting_agent(annotation_result, "user")
    
    print("\nStructured output:")
    print(formatted_output)

    print("Analysis complete.")

if __name__ == "__main__":
    # Example usage
    species = "mouse"
    tissue_type = "larynx"
    marker_list = ["Krt5", "Trp63", "Krt14", "Cdh1", "Epcam", "Vim", "Cdh2", "Acta2", "Msln", "Upk3b", "Krt8", "Krt18", "Foxj1", "Tubb4b", "Scgb1a1", "Scgb3a2", "Muc5ac", "Muc5b", "Ltf"]
    additional_info = "no"
    
    run_analysis(species, tissue_type, marker_list, additional_info)