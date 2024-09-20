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

def final_annotation(agent, prompt):
    response = agent(prompt, "user")
    return response

def coupling_validation(agent, annotation_result, species, tissue_type, marker_list, additional_info):
    validation_message = f"""Please validate the following annotation result:

Annotation Result:
{annotation_result}

Context from input:
Species: {species}
Tissue Type: {tissue_type}
Marker List: {', '.join(marker_list)}
Additional Info: {additional_info if additional_info else 'None'}

Validate the annotation based on this context.
"""
    response = agent(validation_message, "final_annotation")
    return response

def format_results(agent, final_annotations):
    return agent(final_annotations, "user")

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

coupling_validator_agent = Agent(system="""
You are a coupling validator for single-cell analysis results. Validate the final annotation results by checking:
1. Consistency between the identified cell type and the provided markers.
2. Alignment of the sub-cell types with the main cell type.
3. Appropriateness of the annotation given the species and tissue type.
4. Proper consideration of the additional information provided during onboarding.

Provide your validation result, highlighting any inconsistencies or areas of concern.
If everything looks good, say "VALIDATION PASSED". Otherwise, say "VALIDATION FAILED" and explain why.
""")

formatting_agent = Agent(system="""
You are a formatting assistant for single-cell analysis results. Convert the final integrated results 
into a structured JSON format. Follow these guidelines:
1. Extract the main cell type and any sub-cell types identified.
2. Include only information explicitly stated in the input.
3. Ensure the output is valid JSON.

Provide the JSON output within triple backticks.
""")

# Main execution
def run_analysis(species, tissue_type, marker_list, additional_info="", max_iterations=3):
    print("Starting analysis process...\n")
    prompt = construct_prompt(species, tissue_type, marker_list, additional_info)

    validation_passed = False
    iteration = 0

    while not validation_passed and iteration < max_iterations:
        iteration += 1
        print(f"\nStarting final annotation (Iteration {iteration})...\n")
        annotation_result = final_annotation(final_annotation_agent, prompt)

        print("Validating annotation...\n")
        validation_result = coupling_validation(coupling_validator_agent, annotation_result, species, tissue_type, marker_list, additional_info)

        if "VALIDATION PASSED" in validation_result:
            validation_passed = True
        else:
            print(f"Validation failed (Attempt {iteration}/{max_iterations}). Refining annotation...\n")
            prompt = f"Previous annotation attempt failed validation. Please address the following feedback and provide an updated annotation:\n\n{validation_result}\n\nOriginal prompt: {prompt}"

    if validation_passed:
        print("Formatting final results...\n")
        formatted_output = format_results(formatting_agent, annotation_result)
        
        print("\nStructured output:")
        print(formatted_output)
    else:
        print(f"Validation failed after {max_iterations} attempts. Please review the annotation results and validation feedback.")
        print(f"Last annotation result:\n{annotation_result}")
        print(f"Last validation feedback:\n{validation_result}")

    print("Analysis complete.")

if __name__ == "__main__":
    # Example usage
    species = "mouse"
    tissue_type = "larynx"
    marker_list = ["Krt5", "Trp63", "Krt14", "Cdh1", "Epcam", "Vim", "Cdh2", "Acta2", "Msln", "Upk3b", "Krt8", "Krt18", "Foxj1", "Tubb4b", "Scgb1a1", "Scgb3a2", "Muc5ac", "Muc5b", "Ltf", "Bpifb1"]
    additional_info = "Focus on epithelial cell types"
    
    run_analysis(species, tissue_type, marker_list, additional_info)
    
    # Uncomment and modify the line below to change the maximum iterations
    # run_analysis(species, tissue_type, marker_list, additional_info, max_iterations=5)