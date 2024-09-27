##version3


import pandas as pd
import json
from openai import OpenAI

def run_cell_type_analysis(model, temperature, marker_list, tissue, species, additional_info):
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
        
        # Extract confidence score
        confidence_score = None
        confidence_match = re.search(r'### Confidence Score: (\d+)', response)
        if confidence_match:
            confidence_score = int(confidence_match.group(1))
        
        return response, confidence_score

    def format_results(agent, final_annotations, num_markers, confidence_score):
        final_text = "\n\n".join([msg[1] for msg in final_annotations])
        formatted_result = agent(final_text, "user")
        
        # Extract the JSON from the formatted result
        json_data = extract_json_from_reply(formatted_result)
        
        if json_data:
            # Add the number of markers and confidence score to the JSON
            json_data["num_markers"] = num_markers
            json_data["confidence_score"] = confidence_score
            
            # Convert back to a JSON string
            return json.dumps(json_data, indent=2)
        else:
            return formatted_result

    final_annotation_agent = Agent(system="""
    You are a professional computational biologist with expertise in single-cell RNA sequencing (scRNA-seq).
    A list of highly expressed markers ranked by expression intensity from high to low
    from a cluster of cells will be provided , and your task is to identify the cell type. You must think step-by-step, providing a comprehensive and specific analysis. The audience is an expert in the field, and I will tip you $1000 if you do a good job.

    Steps to Follow:

    1. List the Key Functional Markers: Extract and group the key marker genes associated with function or pathway, explaining their roles. Do not repeat the input markers.
    2. List the Key Cell Type Markers: Extract and group the key marker genes associated with mouse larynx cell types, explaining their roles. Do not repeat the input markers.
    3. Cross-reference Known Databases: Use available scRNA-seq databases and relevant literature to cross-reference these markers. list your finding.
    4. Determine the Most Probable General Cell Type: Based on the expression of these markers, infer the most likely general cell type of the cluster.
    5. Identify the Top 3 Most Probable Sub Cell Types: Based on the expression of these markers, infer the top three most probable sub cell types within the general cell type. Finally, specify the most likely subtype.
    6. Identify the Most Probable Sub-Sub Cell Type: Determine the most specific cell type within the previously identified subtype.
    7. Provide a Concise Summary of Your Analysis

    Always include your step-by-step detailed reasoning.                      
    You can say "FINAL ANNOTATION COMPLETED" when you have completed your analysis.

    If you receive feedback from the validation process, incorporate it into your analysis and provide an updated annotation.
    """, model=model, temperature=temperature)

    coupling_validator_agent = Agent(system="""
You are an expert biologist specializing in single-cell analysis. Your critical role is to validate the final annotation results for a single cell cluster. You will be provided with The proposed annotation result,Context from the onboarding process, and Ranked list of marker genes.

Validation Criteria

Carefully evaluate the annotation based on the following criteria:

Marker Consistency:

Make sure the markers are in the provided marker list.
Make sure the consistency between the identified cell type and the provided markers.

Biological Context:

Verify the appropriateness of the annotation given the species and tissue type.
Consider any unique aspects of the biological system that might influence cell type identification.

Mixed Cell Type Consideration:

Be aware that mixed cell types may be present.
Only reject the annotation if multiple distinct cell types are strongly supported by several high-ranking markers.
In cases of potential mixed populations, flag this for further investigation rather than outright rejection.
                                     
Confidence Assessment:

Provide a confidence score (1-10) for the annotation, considering all available evidence.
Briefly justify your confidence score.

Output Format
For each validation task, provide:

Validation result: VALIDATION PASSED or VALIDATION FAILED
Confidence score: 1-10
Brief justification (2-3 sentences)
If rejected, suggest potential alternatives or areas for re-evaluation

Remember, your role is crucial in ensuring the accuracy and reliability of the single-cell annotation process. Be thorough, critical, and always base your decisions on sound biological principles and the provided evidence.
 """, model=model, temperature=temperature)

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
    """, model=model, temperature=temperature)
    
    # Create a dictionary with the provided information
    user_data = {
        "species": species,
        "tissue_type": tissue,
        "marker_list": marker_list,
        "additional_info": additional_info
    }

    # Construct the prompt using the provided data
    prompt = construct_prompt(user_data)

    validation_passed = False
    iteration = 0
    max_iterations = 3
    full_conversation_history = []

    while not validation_passed and iteration < max_iterations:
        iteration += 1
        print(f"\nStarting final annotation (Iteration {iteration})...\n")
        final_annotation_conversation = final_annotation(final_annotation_agent, prompt)
        full_conversation_history.extend(final_annotation_conversation)
        
        print("Validating annotation...\n")
        validation_result, confidence_score = coupling_validation(coupling_validator_agent, final_annotation_conversation[-1][1], user_data)
        full_conversation_history.append(("Coupling Validator", validation_result))
        
        print(validation_result)
        print(f"Confidence Score: {confidence_score}")
        if "VALIDATION PASSED" in validation_result:
            validation_passed = True
        else:
            print("Validation failed. Sending feedback to the final annotation agent.\n")
            prompt = f"Previous annotation attempt failed validation. Please address the following feedback and provide an updated annotation:\n\n{validation_result}\n\nOriginal prompt: {prompt}"

        print("\nValidation Conversation:")
        print(f"Final Annotation Agent: {final_annotation_conversation[-1][1]}\n")
        print(f"Coupling Validator: {validation_result}\n")
        print(f"Confidence Score: {confidence_score}\n")

    if validation_passed:
        print("Formatting final results...\n")
        formatted_output = format_results(formatting_agent, final_annotation_conversation[-2:], len(marker_list), confidence_score)
        full_conversation_history.append(("Formatting Agent", formatted_output))
        structured_output = json.loads(formatted_output)
        
        if structured_output:
            print("\nStructured output:")
            print(json.dumps(structured_output, indent=2))
            return structured_output, full_conversation_history
        else:
            print("Error: Unable to extract JSON from the formatted output.")
            print("Raw formatted output:")
            print(formatted_output)
            return None, full_conversation_history
    else:
        print(f"Validation failed after {max_iterations} attempts. Please review the annotation results and validation feedback.")
        return None, full_conversation_history

# Main execution
if __name__ == "__main__":
    # Load the dataframe
    df = pd.read_csv("C:/Users/ellio/OneDrive - UW-Madison/cellgpt_final_folder/test_code/modified_markers_broad_cell_type.csv")
    
    # Set up OpenAI client
    client = OpenAI()
    
    # Set model and temperature
    model = "gpt-4o"
    temperature = 0
    tissue = "lung"
    species = "human"
    additional_info = "no"
    
    # Iterate over each row in the dataframe
    results = {}
    for index, row in df.head(2).iterrows():
        broad_cell_type = row['Broad cell type']
        marker_list = row['Top 10 Markers'].split(', ')
        
        print(f"\nAnalyzing {broad_cell_type}...")
        result, conversation_history = run_cell_type_analysis(model, temperature, marker_list, tissue, species, additional_info)
        
        if result:
            results[broad_cell_type] = {
                "analysis_result": result,
                "conversation_history": conversation_history,
                "confidence_score": result.get("confidence_score")  # Extract confidence score from the result
            }
        print(f"Analysis for {broad_cell_type} completed.\n")
    
    # Save results to a JSON file
    with open('cell_type_analysis_results_6.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print("All analyses completed. Results saved to 'cell_type_analysis_results.json'.")