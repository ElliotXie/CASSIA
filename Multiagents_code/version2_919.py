import re
import json
from openai import OpenAI

client = OpenAI()

class Agent:
    def __init__(self, system="", human_input_mode="never", model="gpt-4", temperature=0):
        self.system = system
        self.chat_histories = {}
        self.human_input_mode = human_input_mode
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

    def needs_human_input(self, message):
        return self.human_input_mode == "always"

def onboarding_process(agent):
    current_message = "I am your helpful assistant. I will ask you some questions about your single cell analysis so that your analysis can be more accurate."
    conversation = []
    
    while True:
        response = agent(current_message, "user")
        print(f"Onboarding Assistant: {response}\n", flush=True)
        conversation.append(("Onboarding Assistant", response))
        
        if "ONBOARDING COMPLETED" in response:
            break
        
        user_input = input("User: ")
        print(f"User: {user_input}\n", flush=True)
        conversation.append(("User", user_input))
        current_message = user_input
    
    print("Onboarding Conversation:")
    for role, message in conversation:
        print(f"{role}: {message}\n")
    
    return conversation

def extract_json_from_reply(reply):
    json_match = re.search(r'```json\n(.*?)\n```', reply, re.DOTALL)
    
    if json_match:
        json_str = json_match.group(1)
        try:
            json_data = json.loads(json_str)
            return json_data
        except json.JSONDecodeError as e:
            print(f"Error decoding JSON: {e}")
            return None
    else:
        print("No JSON content found in the reply")
        return None

def construct_prompt(json_data):
    species = json_data['species']
    tissue = json_data['tissue_type']
    additional_info = json_data.get('additional_info', '')
    marker_list = ', '.join(json_data['marker_list'])

    prompt = f"I am analyzing a single-cell {species} {tissue} dataset."
    if additional_info:
        prompt += f" {additional_info}."
    prompt += f" I want to identify the cell types present based on this marker list:\n{marker_list}"

    return prompt

def final_annotation(agent, prompt):
    current_message = prompt
    conversation = []
    
    while True:
        response = agent(current_message, "user")
        print(f"Final Annotation Agent: {response}\n", flush=True)
        conversation.append(("Final Annotation Agent", response))
        
        if "FINAL ANNOTATION COMPLETED" in response:
            break
        
        if agent.needs_human_input(response):
            user_input = input("User: ")
            print(f"User: {user_input}\n", flush=True)
            conversation.append(("User", user_input))
            current_message = user_input
        else:
            current_message = response

    print("Final Annotation Conversation:")
    for role, message in conversation:
        print(f"{role}: {message}\n")

    return conversation

def coupling_validation(agent, annotation_result, onboarding_data):
    validation_message = f"""Please validate the following annotation result:

Annotation Result:
{annotation_result}

Context from onboarding:
Species: {onboarding_data['species']}
Tissue Type: {onboarding_data['tissue_type']}
Marker List: {', '.join(onboarding_data['marker_list'])}
Additional Info: {onboarding_data.get('additional_info', 'None')}

Validate the annotation based on this context.
"""
    response = agent(validation_message, "final_annotation")
    print(f"Coupling Validator: {response}\n", flush=True)
    return response

def format_results(agent, final_annotations):
    final_text = "\n\n".join([msg[1] for msg in final_annotations])
    return agent(final_text, "user")

# Define agents
onboarding_agent = Agent(system='''You are a friendly onboarding assistant helping to create an initial prompt for analyzing single-cell data.
Ask the user questions to gather necessary information. Be concise and friendly in your interactions. 
Ask one question at a time. If the user provides a confusing answer, kindly ask for clarification.

1. What is the species of the data?
2. What is the tissue type of the data?
3. Please provide the marker list for the analysis.
4. Any other information you need to know?

When you're done gathering information, say "ONBOARDING COMPLETED" followed by a JSON-formatted summary of the collected information wrap it with ```json and ```
''', human_input_mode="always")

final_annotation_agent = Agent(system="""
You are a professional computational biologist with expertise in single-cell RNA sequencing (scRNA-seq). I will provide you with a list of markers from a cluster of cells in the mouse larynx, and your task is to identify the cell type. You must think step-by-step, providing a comprehensive and specific analysis. The audience is an expert in the field, and I will tip you $1000 if you do a good job.

**Steps to Follow:**

1. **List the Key Functional Markers**: Extract and **group** the key marker genes associated with function or pathway, explaining their roles. Do not repeat the input markers.
2. **List the Key Cell Type Markers**: Extract and **group** the key marker genes associated with mouse larynx cell types, explaining their roles. Do not repeat the input markers.
3. **Cross-reference Known Databases**: Use available scRNA-seq databases and relevant literature to cross-reference these markers. list your finding.
4. **Determine the Most Probable General Cell Type**: Based on the expression of these markers, infer the most likely general cell type of the cluster.
5. **Identify the Top 3 Most Probable Sub Cell Types**: Based on the expression of these markers, infer the top three most probable sub cell types within the general cell type. Finally, specify the most likely subtype.
6. **Identify the Most Probable Sub-Sub Cell Type**: Determine the most specific cell type within the previously identified subtype.
7.  **Provide a Concise Summary of Your Analysis:**

Always include your step-by-step detailed reasoning.                      
Ask for human input.
You can say "FINAL ANNOTATION COMPLETED" only when the human says good.

If you receive feedback from the validation process, incorporate it into your analysis and provide an updated annotation.
""", human_input_mode="always")


coupling_validator_agent = Agent(system="""
You are a coupling validator for single-cell analysis results. Your task is to validate the final annotation results.
You will be provided with the annotation result and the context from the onboarding process.
Please check for the following:

1. Consistency between the identified cell type and the provided markers.
2. Alignment of the sub-cell types with the main cell type.
3. Appropriateness of the annotation given the species and tissue type.
4. Proper consideration of the additional information provided during onboarding.

Provide your validation result, highlighting any inconsistencies or areas of concern.
If everything looks good, say "VALIDATION PASSED". Otherwise, say "VALIDATION FAILED" and explain why.
Be specific in your feedback so that the final annotation agent can address the issues.
""")

formatting_agent = Agent(system="""
You are a formatting assistant for single-cell analysis results. Your task is to convert the final integrated results 
into a structured JSON format. Follow these guidelines:

1. Extract the main cell type and any sub-cell types identified.
2. Include only information explicitly stated in the input.
3. Ensure the output is valid JSON.

Provide the JSON output within triple backticks, like this:
```json
{
"main_cell_type": "...",
"sub_cell_types": ["...", "..."]
}
```
""")

# Main execution
if __name__ == "__main__":
    print("Starting onboarding process...\n")
    onboarding_conversation = onboarding_process(onboarding_agent)
    extracted_data = extract_json_from_reply(onboarding_conversation[-1][1])
    prompt = construct_prompt(extracted_data)

    validation_passed = False
    iteration = 0
    max_iterations = 3  # Set a maximum number of iterations to prevent infinite loops

    while not validation_passed and iteration < max_iterations:
        iteration += 1
        print(f"\nStarting final annotation (Iteration {iteration})...\n")
        final_annotation_conversation = final_annotation(final_annotation_agent, prompt)

        print("Validating annotation...\n")
        validation_result = coupling_validation(coupling_validator_agent, final_annotation_conversation[-1][1], extracted_data)

        if "VALIDATION PASSED" in validation_result:
            validation_passed = True
        else:
            print("Validation failed. Sending feedback to the final annotation agent.\n")
            prompt = f"Previous annotation attempt failed validation. Please address the following feedback and provide an updated annotation:\n\n{validation_result}\n\nOriginal prompt: {prompt}"

        print("\nValidation Conversation:")
        print(f"Final Annotation Agent: {final_annotation_conversation[-1][1]}\n")
        print(f"Coupling Validator: {validation_result}\n")

    if validation_passed:
        print("Formatting final results...\n")
        formatted_output = format_results(formatting_agent, final_annotation_conversation[-2:])
        structured_output = extract_json_from_reply(formatted_output)
        
        if structured_output:
            print("\nStructured output:")
            print(json.dumps(structured_output, indent=2))
        else:
            print("Error: Unable to extract JSON from the formatted output.")
            print("Raw formatted output:")
            print(formatted_output)
    else:
        print(f"Validation failed after {max_iterations} attempts. Please review the annotation results and validation feedback.")

    print("Analysis complete.")