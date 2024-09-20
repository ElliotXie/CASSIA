import json
from openai import OpenAI

client = OpenAI()

class Agent:
    def __init__(self, system="", model="gpt-4", temperature=0):
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
You are a professional computational biologist with expertise in single-cell RNA sequencing (scRNA-seq). Analyze the provided marker list and identify the cell type. Think step-by-step, providing a comprehensive and specific analysis.

Steps to Follow:
1. List and group the key functional markers, explaining their roles.
2. List and group the key cell type markers, explaining their roles.
3. Cross-reference known scRNA-seq databases and relevant literature.
4. Determine the most probable general cell type.
5. Identify the top 3 most probable sub cell types.
6. Identify the most probable sub-sub cell type.
7. Provide a concise summary of your analysis.

Include your step-by-step detailed reasoning.
""")

formatting_agent = Agent(system="""
You are a formatting assistant for single-cell analysis results. Convert the final integrated results 
into a structured JSON format. Follow these guidelines:
1. Extract the main cell type and any sub-cell types identified.
2. Include only information explicitly stated in the input.
3. Ensure the output is valid JSON.

Provide the JSON output within triple backticks.
""")

def run_analysis(species, tissue_type, marker_list, additional_info=""):
    print("Starting analysis process...\n")
    prompt = construct_prompt(species, tissue_type, marker_list, additional_info)

    print("Performing final annotation...\n")
    annotation_result = final_annotation_agent(prompt, "user")

    print("Formatting final results...\n")
    formatted_output = formatting_agent(annotation_result, "user")
    
    print("\nStructured output:")
    print(formatted_output)

    print("Analysis complete.")

if __name__ == "__main__":
    # Example usage
    species = "mouse"
    tissue_type = "larynx"
    marker_list = ["Krt5", "Trp63", "Krt14", "Cdh1", "Epcam", "Vim", "Cdh2", "Acta2", "Msln", "Upk3b", "Krt8", "Krt18", "Foxj1", "Tubb4b", "Scgb1a1", "Scgb3a2", "Muc5ac", "Muc5b", "Ltf", "Bpifb1"]
    additional_info = "Focus on epithelial cell types"
    
    run_analysis(species, tissue_type, marker_list, additional_info)