import json
import re
from openai import OpenAI
import os
import anthropic
import requests

def run_cell_type_analysis(model, temperature, marker_list, tissue, species, additional_info):
    client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))

    class Agent:
        def __init__(self, system="", human_input_mode="never", model="gpt-4o", temperature=0):
            self.system = system
            self.chat_histories = {}
            self.model = model
            self.temperature = temperature

        def __call__(self, message, other_agent_id):
            if other_agent_id not in self.chat_histories:
                self.chat_histories[other_agent_id] = [{"role": "system", "content": self.system}] if self.system else []
            
            self.chat_histories[other_agent_id].append({"role": "user", "content": message})
            
            completion = client.chat.completions.create(
                model=self.model,
                temperature=self.temperature,
                messages=self.chat_histories[other_agent_id]
            )
            result = completion.choices[0].message.content
            self.chat_histories[other_agent_id].append({"role": "assistant", "content": result})
            
            return result

    def extract_json_from_reply(reply):
        json_match = re.search(r'```json\n(.*?)\n```', reply, re.DOTALL)
        if json_match:
            try:
                return json.loads(json_match.group(1))
            except json.JSONDecodeError as e:
                print(f"Error decoding JSON: {e}")
        else:
            print("No JSON content found in the reply")
        return None

    def construct_prompt(json_data):
        marker_list = ', '.join(json_data['marker_list'])
        prompt = f"Your task is to annotate a single-cell {json_data['species']} dataset"
        
        if json_data['tissue_type'] and json_data['tissue_type'].lower() not in ['none', 'tissue blind']:
            prompt += f" from {json_data['tissue_type']} tissue"
        
        prompt += f". Please identify the cell type based on this ranked marker list:\n{marker_list}"
        
        if json_data.get('additional_info') and json_data['additional_info'].lower() != "no":
            prompt += f" Below is some additional information about the dataset:\n{json_data['additional_info']}."
        return prompt

    def final_annotation(agent, prompt):
        conversation = []
        while True:
            response = agent(prompt, "user")
            conversation.append(("Final Annotation Agent", response))
            
            if "FINAL ANNOTATION COMPLETED" in response:
                break
            
            prompt = response

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
        return agent(validation_message, "final_annotation")

    def format_results(agent, final_annotations, num_markers):
        final_text = "\n\n".join([msg[1] for msg in final_annotations])
        formatted_result = agent(final_text, "user")
        
        json_data = extract_json_from_reply(formatted_result)
        if json_data:
            json_data["num_markers"] = num_markers
            return json.dumps(json_data, indent=2)
        return formatted_result

    final_annotation_system_v1 = """
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
    """

    final_annotation_system_v2 = """
    You are a professional computational biologist with expertise in single-cell RNA sequencing (scRNA-seq).
    A list of highly expressed markers ranked by expression intensity from high to low
    from a cluster of cells will be provided, and your task is to identify the cell type. The tissue of origin is not specified, so you must consider multiple possibilities. You must think step-by-step, providing a comprehensive and specific analysis. The audience is an expert in the field, and you will be rewarded $10000 if you do a good job.

    Steps to Follow:

    1. List the Key Functional Markers: Extract and group the key marker genes associated with function or pathway, explaining their roles.
    2. List the Key Cell Type Markers: Extract and group the key marker genes associated with various cell types, explaining their roles.
    3. Cross-reference Known Databases: Use available scRNA-seq databases and relevant literature to cross-reference these markers.
    4. Determine the possible tissue type: Determine the possible tissue type based on the marker list, and provide a detailed explanation for your reasoning.
    5. Determine the Most Probable General Cell Type: Based on the expression of these markers, infer the most likely general cell type of the cluster.
    6. Identify the Top 3 Most Probable Sub Cell Types: Based on the expression of these markers, infer the top three most probable sub cell types. Rank them from most likely to least likely. Finally, specify the most likely subtype based on the markers.
    7. Provide a Concise Summary of Your Analysis

    Always include your step-by-step detailed reasoning.                      
    You can say "FINAL ANNOTATION COMPLETED" when you have completed your analysis.

    If you receive feedback from the validation process, incorporate it into your analysis and provide an updated annotation.
    """

    coupling_validator_system_v1 = """
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
    """

    coupling_validator_system_v2 = """
    You are an expert biologist specializing in single-cell analysis. Your critical role is to validate the final annotation results for a cell cluster where the tissue of origin is not specified. You will be provided with the proposed annotation result and a ranked list of marker genes it used.

    Below are steps to follow:
                                    
    1. Marker Consistency: Make sure the markers are in the provided marker list.
       Ensure consistency between the identified cell type and the provided markers.

    2. Tissue-Agnostic Validation: 
       Ensure that the suggested possible tissues of origin are consistent with the marker expression.

    3. Mixed Cell Type Consideration:
       Be aware that mixed cell types may be present. Only raise this point if multiple distinct cell types are strongly supported by several high-ranking markers. In cases of potential mixed populations, flag this for further investigation rather than outright rejection.
                                        
    Output Format: 
                                        
    If pass:
    Validation result: VALIDATION PASSED

    If failed:
    Validation result: VALIDATION FAILED
    Feedback: give detailed feedback and instruction for revising the annotation
    """

    is_tissue_blind = tissue.lower() in ['none', 'tissue blind'] if tissue else True

    formatting_system_tissue_blind = """
    You are a formatting assistant for single-cell analysis results. Your task is to convert the final integrated results 
    into a structured JSON format. Follow these guidelines:

    1. Extract the main cell type and any sub-cell types identified.
    2. Include only information explicitly stated in the input.
    3. If there are possible mixed cell types highlighted, list them.
    4. Include possible tissues.

    Provide the JSON output within triple backticks, like this:
    ```json
    {
    "main_cell_type": "...",
    "sub_cell_types": ["...", "..."],
    "possible_mixed_cell_types": ["...", "..."],
    "possible_tissues": ["...", "..."]
    }
    ```
    """

    formatting_system_non_tissue_blind = """
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
    """

    formatting_system_failed = """
    You are a formatting assistant for single-cell analysis results. Your task is to convert the final integrated results 
    into a structured JSON format, with special consideration for uncertain or conflicting annotations. Follow these guidelines:

    1. The analysis failed after multiple attempts. Please try to extract as much information as possible. Summarize what has gone wrong and what has been tried.
    2. Provide a detailed feedback on why the analysis failed, and what has been tried and why it did not work.
    3. Finally, provide a detailed step-by-step reasoning of how to fix the analysis.

    Provide the JSON output within triple backticks, like this:
    ```json
    {
    "main_cell_type": "if any",
    "sub_cell_types": "if any",
    "possible_cell_types": "if any",
    "feedback": "...",
    "next_steps": "..."
    }
    ```
    """

    final_annotation_agent = Agent(
        system=final_annotation_system_v2 if is_tissue_blind else final_annotation_system_v1, 
        model=model, 
        temperature=temperature
    )

    coupling_validator_agent = Agent(
        system=coupling_validator_system_v2.strip() if is_tissue_blind else coupling_validator_system_v1.strip(), 
        model=model, 
        temperature=temperature
    )

    formatting_agent = Agent(
        system=formatting_system_tissue_blind if is_tissue_blind else formatting_system_non_tissue_blind,
        model=model,
        temperature=temperature
    )
    
    user_data = {
        "species": species,
        "tissue_type": tissue,
        "marker_list": marker_list,
    }
    if additional_info and additional_info.lower() != "no":
        user_data["additional_info"] = additional_info

    prompt = construct_prompt(user_data)

    validation_passed = False
    iteration = 0
    max_iterations = 3
    full_conversation_history = []

    while not validation_passed and iteration < max_iterations:
        iteration += 1
        
        if iteration > 1:
            prompt = f"""Previous annotation attempt failed validation. Please review your previous response and the validation feedback, then provide an updated annotation:

Previous response:
{final_annotation_conversation[-1][1]}

Validation feedback:
{validation_result}

Original prompt:
{prompt}

Please provide an updated annotation addressing the validation feedback."""

        final_annotation_conversation = final_annotation(final_annotation_agent, prompt)
        full_conversation_history.extend(final_annotation_conversation)
        
        validation_result = coupling_validation(coupling_validator_agent, final_annotation_conversation[-1][1], user_data)
        full_conversation_history.append(("Coupling Validator", validation_result))
        
        if "VALIDATION PASSED" in validation_result:
            validation_passed = True

    formatting_system = formatting_system_tissue_blind if is_tissue_blind else formatting_system_non_tissue_blind if validation_passed else formatting_system_failed

    formatting_agent.system = formatting_system
    formatted_output = format_results(formatting_agent, final_annotation_conversation[-2:], len(marker_list))
    full_conversation_history.append(("Formatting Agent", formatted_output))
    structured_output = json.loads(formatted_output)
    
    if structured_output:
        structured_output["iterations"] = iteration
        return structured_output, full_conversation_history
    else:
        print("Error: Unable to extract JSON from the formatted output.")
        print("Raw formatted output:")
        print(formatted_output)
        return None, full_conversation_history
    



def run_cell_type_analysis_claude(model, temperature, marker_list, tissue, species, additional_info):
    client = anthropic.Anthropic(api_key=os.environ.get("ANTHROPIC_API_KEY"))

    class Agent:
        def __init__(self, system="", human_input_mode="never", model=model, temperature=0):
            self.system = system
            self.chat_histories = {}
            self.model = model
            self.temperature = temperature

        def __call__(self, message, other_agent_id):
            if other_agent_id not in self.chat_histories:
                self.chat_histories[other_agent_id] = []
            
            # Combine system prompt with user message if it's the first message
            if not self.chat_histories[other_agent_id] and self.system:
                full_message = f"{self.system}\n\n{message}"
            else:
                full_message = message

            response = client.messages.create(
                model=self.model,
                temperature=self.temperature,
                max_tokens=7000,
                system=self.system,  # System message goes here as a separate parameter
                messages=[
                    *[{"role": msg["role"], "content": msg["content"]} 
                    for msg in self.chat_histories[other_agent_id]],
                    {"role": "user", "content": full_message}
                ]
            )
            
            result = response.content[0].text
            self.chat_histories[other_agent_id].append({"role": "user", "content": full_message})
            self.chat_histories[other_agent_id].append({"role": "assistant", "content": result})
            
            return result

    def extract_json_from_reply(reply):
        json_match = re.search(r'```json\n(.*?)\n```', reply, re.DOTALL)
        if json_match:
            try:
                return json.loads(json_match.group(1))
            except json.JSONDecodeError as e:
                print(f"Error decoding JSON: {e}")
        else:
            print("No JSON content found in the reply")
        return None

    def construct_prompt(json_data):
        marker_list = ', '.join(json_data['marker_list'])
        prompt = f"Your task is to annotate a single-cell {json_data['species']} dataset"
        
        if json_data['tissue_type'] and json_data['tissue_type'].lower() not in ['none', 'tissue blind']:
            prompt += f" from {json_data['tissue_type']} tissue"
        
        prompt += f". Please identify the cell type based on this ranked marker list:\n{marker_list}"
        
        if json_data.get('additional_info') and json_data['additional_info'].lower() != "no":
            prompt += f" Below is some additional information about the dataset:\n{json_data['additional_info']}."
        return prompt

    def final_annotation(agent, prompt):
        conversation = []
        while True:
            response = agent(prompt, "user")
            conversation.append(("Final Annotation Agent", response))
            
            if "FINAL ANNOTATION COMPLETED" in response:
                break
            
            prompt = response

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
        return agent(validation_message, "final_annotation")

    def format_results(agent, final_annotations, num_markers):
        final_text = "\n\n".join([msg[1] for msg in final_annotations])
        formatted_result = agent(final_text, "user")
        
        json_data = extract_json_from_reply(formatted_result)
        if json_data:
            json_data["num_markers"] = num_markers
            return json.dumps(json_data, indent=2)
        return formatted_result

    final_annotation_system_v1 = """
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
    """

    final_annotation_system_v2 = """
    You are a professional computational biologist with expertise in single-cell RNA sequencing (scRNA-seq).
    A list of highly expressed markers ranked by expression intensity from high to low
    from a cluster of cells will be provided, and your task is to identify the cell type. The tissue of origin is not specified, so you must consider multiple possibilities. You must think step-by-step, providing a comprehensive and specific analysis. The audience is an expert in the field, and you will be rewarded $10000 if you do a good job.

    Steps to Follow:

    1. List the Key Functional Markers: Extract and group the key marker genes associated with function or pathway, explaining their roles.
    2. List the Key Cell Type Markers: Extract and group the key marker genes associated with various cell types, explaining their roles.
    3. Cross-reference Known Databases: Use available scRNA-seq databases and relevant literature to cross-reference these markers.
    4. Determine the possible tissue type: Determine the possible tissue type based on the marker list, and provide a detailed explanation for your reasoning.
    5. Determine the Most Probable General Cell Type: Based on the expression of these markers, infer the most likely general cell type of the cluster.
    6. Identify the Top 3 Most Probable Sub Cell Types: Based on the expression of these markers, infer the top three most probable sub cell types. Rank them from most likely to least likely. Finally, specify the most likely subtype based on the markers.
    7. Provide a Concise Summary of Your Analysis

    Always include your step-by-step detailed reasoning.                      
    You can say "FINAL ANNOTATION COMPLETED" when you have completed your analysis.

    If you receive feedback from the validation process, incorporate it into your analysis and provide an updated annotation.
    """

    coupling_validator_system_v1 = """
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
    """

    coupling_validator_system_v2 = """
    You are an expert biologist specializing in single-cell analysis. Your critical role is to validate the final annotation results for a cell cluster where the tissue of origin is not specified. You will be provided with the proposed annotation result and a ranked list of marker genes it used.

    Below are steps to follow:
                                    
    1. Marker Consistency: Make sure the markers are in the provided marker list.
       Ensure consistency between the identified cell type and the provided markers.

    2. Tissue-Agnostic Validation: 
       Ensure that the suggested possible tissues of origin are consistent with the marker expression.

    3. Mixed Cell Type Consideration:
       Be aware that mixed cell types may be present. Only raise this point if multiple distinct cell types are strongly supported by several high-ranking markers. In cases of potential mixed populations, flag this for further investigation rather than outright rejection.
                                        
    Output Format: 
                                        
    If pass:
    Validation result: VALIDATION PASSED

    If failed:
    Validation result: VALIDATION FAILED
    Feedback: give detailed feedback and instruction for revising the annotation
    """

    is_tissue_blind = tissue.lower() in ['none', 'tissue blind'] if tissue else True

    formatting_system_tissue_blind = """
    You are a formatting assistant for single-cell analysis results. Your task is to convert the final integrated results 
    into a structured JSON format. Follow these guidelines:

    1. Extract the main cell type and any sub-cell types identified.
    2. Include only information explicitly stated in the input.
    3. If there are possible mixed cell types highlighted, list them.
    4. Include possible tissues.

    Provide the JSON output within triple backticks, like this:
    ```json
    {
    "main_cell_type": "...",
    "sub_cell_types": ["...", "..."],
    "possible_mixed_cell_types": ["...", "..."],
    "possible_tissues": ["...", "..."]
    }
    ```
    """

    formatting_system_non_tissue_blind = """
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
    """

    formatting_system_failed = """
    You are a formatting assistant for single-cell analysis results. Your task is to convert the final integrated results 
    into a structured JSON format, with special consideration for uncertain or conflicting annotations. Follow these guidelines:

    1. The analysis failed after multiple attempts. Please try to extract as much information as possible. Summarize what has gone wrong and what has been tried.
    2. Provide a detailed feedback on why the analysis failed, and what has been tried and why it did not work.
    3. Finally, provide a detailed step-by-step reasoning of how to fix the analysis.

    Provide the JSON output within triple backticks, like this:
    ```json
    {
    "main_cell_type": "if any",
    "sub_cell_types": "if any",
    "possible_cell_types": "if any",
    "feedback": "...",
    "next_steps": "..."
    }
    ```
    """

    final_annotation_agent = Agent(
        system=final_annotation_system_v2 if is_tissue_blind else final_annotation_system_v1, 
        temperature=temperature
    )

    coupling_validator_agent = Agent(
        system=coupling_validator_system_v2.strip() if is_tissue_blind else coupling_validator_system_v1.strip(), 
        temperature=temperature
    )

    formatting_agent = Agent(
        system=formatting_system_tissue_blind if is_tissue_blind else formatting_system_non_tissue_blind,
        temperature=temperature
    )
    
    user_data = {
        "species": species,
        "tissue_type": tissue,
        "marker_list": marker_list,
    }
    if additional_info and additional_info.lower() != "no":
        user_data["additional_info"] = additional_info

    prompt = construct_prompt(user_data)

    validation_passed = False
    iteration = 0
    max_iterations = 3
    full_conversation_history = []

    while not validation_passed and iteration < max_iterations:
        iteration += 1
        
        if iteration > 1:
            prompt = f"""Previous annotation attempt failed validation. Please review your previous response and the validation feedback, then provide an updated annotation:

Previous response:
{final_annotation_conversation[-1][1]}

Validation feedback:
{validation_result}

Original prompt:
{prompt}

Please provide an updated annotation addressing the validation feedback."""

        final_annotation_conversation = final_annotation(final_annotation_agent, prompt)
        full_conversation_history.extend(final_annotation_conversation)
        
        validation_result = coupling_validation(coupling_validator_agent, final_annotation_conversation[-1][1], user_data)
        full_conversation_history.append(("Coupling Validator", validation_result))
        
        if "VALIDATION PASSED" in validation_result:
            validation_passed = True

    formatting_system = formatting_system_tissue_blind if is_tissue_blind else formatting_system_non_tissue_blind if validation_passed else formatting_system_failed

    formatting_agent.system = formatting_system
    formatted_output = format_results(formatting_agent, final_annotation_conversation[-2:], len(marker_list))
    full_conversation_history.append(("Formatting Agent", formatted_output))
    structured_output = json.loads(formatted_output)
    
    if structured_output:
        structured_output["iterations"] = iteration
        return structured_output, full_conversation_history
    else:
        print("Error: Unable to extract JSON from the formatted output.")
        print("Raw formatted output:")
        print(formatted_output)
        return None, full_conversation_history
    





def run_cell_type_analysis_openrouter(model, temperature, marker_list, tissue, species, additional_info):

    class Agent:
        def __init__(self, system="", human_input_mode="never", model="gpt-4o", temperature=0):
            self.system = system
            self.chat_histories = {}
            self.model = model
            self.temperature = temperature

        def __call__(self, message, other_agent_id):
            if other_agent_id not in self.chat_histories:
                self.chat_histories[other_agent_id] = [{"role": "system", "content": self.system}] if self.system else []
            
            self.chat_histories[other_agent_id].append({"role": "user", "content": message})
            
            response = requests.post(
                url="https://openrouter.ai/api/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {os.environ.get('OPENROUTER_API_KEY')}",
                    "HTTP-Referer": "https://localhost:5000",
                    "Content-Type": "application/json"
                },
                json={
                    "model": self.model,
                    "temperature": self.temperature,
                    "messages": self.chat_histories[other_agent_id]
                }
            )
            
            if response.status_code == 200:
                response_data = response.json()
                result = response_data['choices'][0]['message']['content']
                self.chat_histories[other_agent_id].append({"role": "assistant", "content": result})
                return result
            else:
                raise Exception(f"OpenRouter API error: {response.status_code}")

    def extract_json_from_reply(reply):
        json_match = re.search(r'```json\n(.*?)\n```', reply, re.DOTALL)
        if json_match:
            try:
                return json.loads(json_match.group(1))
            except json.JSONDecodeError as e:
                print(f"Error decoding JSON: {e}")
        else:
            print("No JSON content found in the reply")
        return None

    def construct_prompt(json_data):
        marker_list = ', '.join(json_data['marker_list'])
        prompt = f"Your task is to annotate a single-cell {json_data['species']} dataset"
        
        if json_data['tissue_type'] and json_data['tissue_type'].lower() not in ['none', 'tissue blind']:
            prompt += f" from {json_data['tissue_type']} tissue"
        
        prompt += f". Please identify the cell type based on this ranked marker list:\n{marker_list}"
        
        if json_data.get('additional_info') and json_data['additional_info'].lower() != "no":
            prompt += f" Below is some additional information about the dataset:\n{json_data['additional_info']}."
        return prompt

    def final_annotation(agent, prompt):
        conversation = []
        while True:
            response = agent(prompt, "user")
            conversation.append(("Final Annotation Agent", response))
            
            if "FINAL ANNOTATION COMPLETED" in response:
                break
            
            prompt = response

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
        return agent(validation_message, "final_annotation")

    def format_results(agent, final_annotations, num_markers):
        final_text = "\n\n".join([msg[1] for msg in final_annotations])
        formatted_result = agent(final_text, "user")
        
        json_data = extract_json_from_reply(formatted_result)
        if json_data:
            json_data["num_markers"] = num_markers
            return json.dumps(json_data, indent=2)
        return formatted_result

    final_annotation_system_v1 = """
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
    """

    final_annotation_system_v2 = """
    You are a professional computational biologist with expertise in single-cell RNA sequencing (scRNA-seq).
    A list of highly expressed markers ranked by expression intensity from high to low
    from a cluster of cells will be provided, and your task is to identify the cell type. The tissue of origin is not specified, so you must consider multiple possibilities. You must think step-by-step, providing a comprehensive and specific analysis. The audience is an expert in the field, and you will be rewarded $10000 if you do a good job.

    Steps to Follow:

    1. List the Key Functional Markers: Extract and group the key marker genes associated with function or pathway, explaining their roles.
    2. List the Key Cell Type Markers: Extract and group the key marker genes associated with various cell types, explaining their roles.
    3. Cross-reference Known Databases: Use available scRNA-seq databases and relevant literature to cross-reference these markers.
    4. Determine the possible tissue type: Determine the possible tissue type based on the marker list, and provide a detailed explanation for your reasoning.
    5. Determine the Most Probable General Cell Type: Based on the expression of these markers, infer the most likely general cell type of the cluster.
    6. Identify the Top 3 Most Probable Sub Cell Types: Based on the expression of these markers, infer the top three most probable sub cell types. Rank them from most likely to least likely. Finally, specify the most likely subtype based on the markers.
    7. Provide a Concise Summary of Your Analysis

    Always include your step-by-step detailed reasoning.                      
    You can say "FINAL ANNOTATION COMPLETED" when you have completed your analysis.

    If you receive feedback from the validation process, incorporate it into your analysis and provide an updated annotation.
    """

    coupling_validator_system_v1 = """
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
    """

    coupling_validator_system_v2 = """
    You are an expert biologist specializing in single-cell analysis. Your critical role is to validate the final annotation results for a cell cluster where the tissue of origin is not specified. You will be provided with the proposed annotation result and a ranked list of marker genes it used.

    Below are steps to follow:
                                    
    1. Marker Consistency: Make sure the markers are in the provided marker list.
       Ensure consistency between the identified cell type and the provided markers.

    2. Tissue-Agnostic Validation: 
       Ensure that the suggested possible tissues of origin are consistent with the marker expression.

    3. Mixed Cell Type Consideration:
       Be aware that mixed cell types may be present. Only raise this point if multiple distinct cell types are strongly supported by several high-ranking markers. In cases of potential mixed populations, flag this for further investigation rather than outright rejection.
                                        
    Output Format: 
                                        
    If pass:
    Validation result: VALIDATION PASSED

    If failed:
    Validation result: VALIDATION FAILED
    Feedback: give detailed feedback and instruction for revising the annotation
    """

    is_tissue_blind = tissue.lower() in ['none', 'tissue blind'] if tissue else True

    formatting_system_tissue_blind = """
    You are a formatting assistant for single-cell analysis results. Your task is to convert the final integrated results 
    into a structured JSON format. Follow these guidelines:

    1. Extract the main cell type and any sub-cell types identified.
    2. Include only information explicitly stated in the input.
    3. If there are possible mixed cell types highlighted, list them.
    4. Include possible tissues.

    Provide the JSON output within triple backticks, like this:
    ```json
    {
    "main_cell_type": "...",
    "sub_cell_types": ["...", "..."],
    "possible_mixed_cell_types": ["...", "..."],
    "possible_tissues": ["...", "..."]
    }
    ```
    """

    formatting_system_non_tissue_blind = """
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
    """

    formatting_system_failed = """
    You are a formatting assistant for single-cell analysis results. Your task is to convert the final integrated results 
    into a structured JSON format, with special consideration for uncertain or conflicting annotations. Follow these guidelines:

    1. The analysis failed after multiple attempts. Please try to extract as much information as possible. Summarize what has gone wrong and what has been tried.
    2. Provide a detailed feedback on why the analysis failed, and what has been tried and why it did not work.
    3. Finally, provide a detailed step-by-step reasoning of how to fix the analysis.

    Provide the JSON output within triple backticks, like this:
    ```json
    {
    "main_cell_type": "if any",
    "sub_cell_types": "if any",
    "possible_cell_types": "if any",
    "feedback": "...",
    "next_steps": "..."
    }
    ```
    """

    final_annotation_agent = Agent(
        system=final_annotation_system_v2 if is_tissue_blind else final_annotation_system_v1, 
        model=model, 
        temperature=temperature
    )

    coupling_validator_agent = Agent(
        system=coupling_validator_system_v2.strip() if is_tissue_blind else coupling_validator_system_v1.strip(), 
        model=model, 
        temperature=temperature
    )

    formatting_agent = Agent(
        system=formatting_system_tissue_blind if is_tissue_blind else formatting_system_non_tissue_blind,
        model=model,
        temperature=temperature
    )
    
    user_data = {
        "species": species,
        "tissue_type": tissue,
        "marker_list": marker_list,
    }
    if additional_info and additional_info.lower() != "no":
        user_data["additional_info"] = additional_info

    prompt = construct_prompt(user_data)

    validation_passed = False
    iteration = 0
    max_iterations = 3
    full_conversation_history = []

    while not validation_passed and iteration < max_iterations:
        iteration += 1
        
        if iteration > 1:
            prompt = f"""Previous annotation attempt failed validation. Please review your previous response and the validation feedback, then provide an updated annotation:

Previous response:
{final_annotation_conversation[-1][1]}

Validation feedback:
{validation_result}

Original prompt:
{prompt}

Please provide an updated annotation addressing the validation feedback."""

        final_annotation_conversation = final_annotation(final_annotation_agent, prompt)
        full_conversation_history.extend(final_annotation_conversation)
        
        validation_result = coupling_validation(coupling_validator_agent, final_annotation_conversation[-1][1], user_data)
        full_conversation_history.append(("Coupling Validator", validation_result))
        
        if "VALIDATION PASSED" in validation_result:
            validation_passed = True

    formatting_system = formatting_system_tissue_blind if is_tissue_blind else formatting_system_non_tissue_blind if validation_passed else formatting_system_failed

    formatting_agent.system = formatting_system
    formatted_output = format_results(formatting_agent, final_annotation_conversation[-2:], len(marker_list))
    full_conversation_history.append(("Formatting Agent", formatted_output))
    structured_output = json.loads(formatted_output)
    
    if structured_output:
        structured_output["iterations"] = iteration
        return structured_output, full_conversation_history
    else:
        print("Error: Unable to extract JSON from the formatted output.")
        print("Raw formatted output:")
        print(formatted_output)
        return None, full_conversation_history