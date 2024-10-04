import pandas as pd
import json
from openai import OpenAI
import re

def run_full_cell_type_analysis(df_path, output_json_name="cell_type_analysis_results.json", model="gpt-4", temperature=0, tissue="lung", species="human", additional_info=None, celltype_column="Broad cell type"):
    def run_cell_type_analysis(model, temperature, marker_list, tissue, species, additional_info):
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
            additional_info = json_data.get('additional_info')
            marker_list = ', '.join(json_data['marker_list'])

            prompt = f"Your task is to annotate a single-cell {species} {tissue} dataset. Please identify the cell type based on this ranked marker list:\n{marker_list}"
            
            if additional_info and additional_info.lower() != "no":
                prompt += f" below is some additional information about the dataset:\n{additional_info}."

            return prompt

        def final_annotation(agent, prompt):
            current_message = prompt
            conversation = []
            
            while True:
                print("current prompt: ", current_message)
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

        Context:

        Marker List: {', '.join(onboarding_data['marker_list'])}
        Additional Info: {onboarding_data.get('additional_info', 'None')}

        Validate the annotation based on this context.
        """
            response = agent(validation_message, "final_annotation")
            print(f"Coupling Validator: {response}\n", flush=True)
            
            return response

        def format_results(agent, final_annotations, num_markers):
            final_text = "\n\n".join([msg[1] for msg in final_annotations])
            formatted_result = agent(final_text, "user")
            
            # Extract the JSON from the formatted result
            json_data = extract_json_from_reply(formatted_result)
            
            if json_data:
                # Add the number of markers to the JSON
                json_data["num_markers"] = num_markers
                
                # Convert back to a JSON string
                return json.dumps(json_data, indent=2)
            else:
                return formatted_result

        final_annotation_agent = Agent(system="""
        You are a professional computational biologist with expertise in single-cell RNA sequencing (scRNA-seq).
        A list of highly expressed markers ranked by expression intensity from high to low
        from a cluster of cells will be provided , and your task is to identify the cell type. You must think step-by-step, providing a comprehensive and specific analysis. The audience is an expert in the field, and you will be rewarded $10000 if you do a good job.

        Steps to Follow:

        1. List the Key Functional Markers: Extract and group the key marker genes associated with function or pathway, explaining their roles.
        2. List the Key Cell Type Markers: Extract and group the key marker genes associated with target tissue cell types, explaining their roles.
        3. Cross-reference Known Databases: Use available scRNA-seq databases and relevant literature to cross-reference these markers.
        4. Determine the Most Probable General Cell Type: Based on the expression of these markers, infer the most likely general cell type of the cluster.
        5. Identify the Top 3 Most Probable Sub Cell Types: Based on the expression of these markers, infer the top three most probable sub cell types within the general cell type. Rank them from most likely to least likely. Finally, specify the most likely subtype based on the markers.
        6. Provide a Concise Summary of Your Analysis

        Always include your step-by-step detailed reasoning.                      
        You can say "FINAL ANNOTATION COMPLETED" when you have completed your analysis.

        If you receive feedback from the validation process, incorporate it into your analysis and provide an updated annotation.
        """, model=model, temperature=temperature)

        coupling_validator_agent = Agent(system="""
    You are an expert biologist specializing in single-cell analysis. Your critical role is to validate the final annotation results for a cell cluster. You will be provided with The proposed annotation result, and a Ranked list of marker genes it used.

    
    Below are steps to follow:
                                        
    1.Marker Consistency: Make sure the markers are in the provided marker list.
    Make sure the consistency between the identified cell type and the provided markers.
                                         

    2.Mixed Cell Type Consideration:
Be aware that mixed cell types may be present. Only raise this point if multiple distinct cell types are strongly supported by several high-ranking markers. In cases of potential mixed populations, flag this for further investigation rather than outright rejection.
                                         
    Output Format: 
                                         
    if pass,

    Validation result: VALIDATION PASSED

    If failed,
                                                               
    Validation result: VALIDATION FAILED
    Feedback: give detailed feedback and instruction for revising the annotation
    
     """.strip(), model=model, temperature=temperature)

        formatting_agent = Agent(system="""
        You are a formatting assistant for single-cell analysis results. Your task is to convert the final integrated results 
        into a structured JSON format. Follow these guidelines:

        1. Extract the main cell type and any sub-cell types identified.
        2. Include only information explicitly stated in the input.
        3. If there are possible mixed cell types highlighted, list them.

        Provide the JSON output within triple backticks, like this:
        ```json
        {
        "main_cell_type": "...",
        "sub_cell_types": ["...", "..."],
        "possible_mixed_cell_types": ["...", "..."]
        }
        ```
        """, model=model, temperature=temperature)
        
        # Create a dictionary with the provided information
        user_data = {
            "species": species,
            "tissue_type": tissue,
            "marker_list": marker_list,
        }
        if additional_info and additional_info.lower() != "no":
            user_data["additional_info"] = additional_info

        # Construct the prompt using the provided data
        prompt = construct_prompt(user_data)

        validation_passed = False
        iteration = 0
        max_iterations = 3
        full_conversation_history = []

        while not validation_passed and iteration < max_iterations:
            iteration += 1
            print(f"\nStarting final annotation (Iteration {iteration})...\n")
            
            if iteration > 1:
                # Update the prompt with previous response and validation feedback
                prompt = f"""Previous annotation attempt failed validation. Please review your previous response and the validation feedback, then provide an updated annotation:

Previous response:
{final_annotation_conversation[-1][1]}

Validation feedback:
{validation_result}

Original prompt:
{prompt}

Please provide an updated annotation addressing the validation feedback."""

            final_annotation_conversation = final_annotation(final_annotation_agent, prompt)
            print("updated prompt: ", prompt)
            full_conversation_history.extend(final_annotation_conversation)
            
            print("Validating annotation...\n")
            validation_result = coupling_validation(coupling_validator_agent, final_annotation_conversation[-1][1], user_data)
            full_conversation_history.append(("Coupling Validator", validation_result))
            
            print(validation_result)
            if "VALIDATION PASSED" in validation_result:
                validation_passed = True
            else:
                print("Validation failed. Will update prompt for next iteration.\n")

            print("\nValidation Conversation:")
            print(f"Final Annotation Agent: {final_annotation_conversation[-1][1]}\n")
            print(f"Coupling Validator: {validation_result}\n")

        if validation_passed:
            print("Formatting final results...\n")
            formatted_output = format_results(formatting_agent, final_annotation_conversation[-2:], len(marker_list))
            full_conversation_history.append(("Formatting Agent", formatted_output))
            structured_output = json.loads(formatted_output)
            
            if structured_output:
                structured_output["iterations"] = iteration  # Add the number of iterations to the structured output
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
            return {"iterations": iteration}, full_conversation_history  # Return iteration count even if validation fails

    # Load the dataframe
    df = pd.read_csv(df_path)
    
    # Set up OpenAI client
    client = OpenAI()
    
    # Iterate over each row in the dataframe
    results = {}
    for index, row in df.iterrows():
        cell_type = row[celltype_column]
        marker_list = row['top_markers'].split(', ')
        
        print(f"\nAnalyzing {cell_type}...")
        result, conversation_history = run_cell_type_analysis(model, temperature, marker_list, tissue, species, additional_info)
        
        if result:
            results[cell_type] = {
                "analysis_result": result,
                "conversation_history": conversation_history,
                "iterations": result.get("iterations", 1)
            }
        print(f"Analysis for {cell_type} completed.\n")
    
    # Save results to the specified JSON file
    with open(output_json_name, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"All analyses completed. Results saved to '{output_json_name}'.")
    
    return results