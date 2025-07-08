import json
import re
import os
from llm_utils import call_llm

def run_cell_type_analysis_unified(model, temperature, marker_list, tissue, species, additional_info, provider="openai"):
    """
    Unified function for cell type analysis using any provider via call_llm.
    """
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
        marker_list_str = ', '.join(json_data['marker_list'])
        prompt = f"Your task is to annotate a single-cell {json_data['species']} dataset"
        if json_data['tissue_type'] and json_data['tissue_type'].lower() not in ['none', 'tissue blind']:
            prompt += f" from {json_data['tissue_type']} tissue"
        prompt += f". Please identify the cell type based on this ranked marker list:\n{marker_list_str}"
        if json_data.get('additional_info') and json_data['additional_info'].lower() != "no":
            prompt += f" Below is some additional information about the dataset:\n{json_data['additional_info']}."
        return prompt

    def final_annotation(prompt, system_prompt, model, temperature, provider):
        conversation = []
        current_prompt = prompt
        for _ in range(10):  # safeguard max iterations
            response = call_llm(
                prompt=current_prompt,
                provider=provider,
                model=model,
                temperature=temperature,
                system_prompt=system_prompt
            )
            conversation.append(("Final Annotation Agent", response))
            if "FINAL ANNOTATION COMPLETED" in response:
                break
            current_prompt = response
        return conversation

    def coupling_validation(annotation_result, onboarding_data, system_prompt, model, temperature, provider):
        # Handle different marker_list formats
        marker_list = onboarding_data.get('marker_list', [])
        
        # If marker_list is a list with one element that contains commas, it's likely a single string with all genes
        if (isinstance(marker_list, list) and 
            len(marker_list) == 1 and 
            isinstance(marker_list[0], str) and 
            ',' in marker_list[0]):
            
            # Split the single string into individual genes
            import re
            markers = re.split(r',\s*', marker_list[0])
            markers = [m.strip() for m in markers if m.strip()]
            marker_list_str = ', '.join(markers)
        elif isinstance(marker_list, list):
            # Normal case: list of individual genes
            marker_list_str = ', '.join(marker_list)
        elif isinstance(marker_list, str):
            # Single string case: use as-is
            marker_list_str = marker_list
        else:
            # Fallback
            marker_list_str = str(marker_list)
        
        validation_message = f"""Please validate the following annotation result:\n\nAnnotation Result:\n{annotation_result}\n\nContext:\n\nMarker List: {marker_list_str}\nAdditional Info: {onboarding_data.get('additional_info', 'None')}\n\nValidate the annotation based on this context.\n"""
        response = call_llm(
            prompt=validation_message,
            provider=provider,
            model=model,
            temperature=temperature,
            system_prompt=system_prompt
        )
        return response

    def format_results(final_annotations, num_markers, system_prompt, model, temperature, provider):
        final_text = "\n\n".join([msg[1] for msg in final_annotations])
        formatted_result = call_llm(
            prompt=final_text,
            provider=provider,
            model=model,
            temperature=temperature,
            system_prompt=system_prompt
        )
        json_data = extract_json_from_reply(formatted_result)
        if json_data:
            json_data["num_markers"] = num_markers
            return json.dumps(json_data, indent=2)
        return formatted_result

    # System prompts (unchanged)
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

    1. Extract the main cell type and the three most likely sub-cell types identified from step 4 and step 5 of the Final Annotation Agent response. Even the main cell type is the same as the sub-cell types, you still need to list it as a sub-cell type. Strictly follow the order of the sub-cell types.
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

    is_tissue_blind = tissue.lower() in ['none', 'tissue blind'] if tissue else True
    user_data = {
        "species": species,
        "tissue_type": tissue,
        "marker_list": marker_list,
    }
    if additional_info and additional_info.lower() != "no":
        user_data["additional_info"] = additional_info
    prompt = construct_prompt(user_data)

    # Select system prompts
    final_annotation_system = final_annotation_system_v2 if is_tissue_blind else final_annotation_system_v1
    coupling_validator_system = coupling_validator_system_v2.strip() if is_tissue_blind else coupling_validator_system_v1.strip()
    formatting_system = formatting_system_tissue_blind if is_tissue_blind else formatting_system_non_tissue_blind

    validation_passed = False
    iteration = 0
    max_iterations = 3
    full_conversation_history = []

    while not validation_passed and iteration < max_iterations:
        iteration += 1
        if iteration > 1:
            prompt = f"""Previous annotation attempt failed validation. Please review your previous response and the validation feedback, then provide an updated annotation:\n\nPrevious response:\n{final_annotation_conversation[-1][1]}\n\nValidation feedback:\n{validation_result}\n\nOriginal prompt:\n{prompt}\n\nPlease provide an updated annotation addressing the validation feedback."""
        final_annotation_conversation = final_annotation(prompt, final_annotation_system, model, temperature, provider)
        full_conversation_history.extend(final_annotation_conversation)
        validation_result = coupling_validation(final_annotation_conversation[-1][1], user_data, coupling_validator_system, model, temperature, provider)
        full_conversation_history.append(("Coupling Validator", validation_result))
        if "VALIDATION PASSED" in validation_result:
            validation_passed = True

    formatting_system_final = formatting_system if validation_passed else formatting_system_failed
    formatted_output = format_results(final_annotation_conversation[-2:], len(marker_list), formatting_system_final, model, temperature, provider)
    full_conversation_history.append(("Formatting Agent", formatted_output))
    try:
        structured_output = json.loads(formatted_output)
    except Exception:
        structured_output = None
    if structured_output:
        structured_output["iterations"] = iteration
        return structured_output, full_conversation_history
    else:
        print("Error: Unable to extract JSON from the formatted output.")
        print("Raw formatted output:")
        print(formatted_output)
        return None, full_conversation_history

# Backward-compatible wrappers

def run_cell_type_analysis(model, temperature, marker_list, tissue, species, additional_info):
    return run_cell_type_analysis_unified(model, temperature, marker_list, tissue, species, additional_info, provider="openai")

def run_cell_type_analysis_claude(model, temperature, marker_list, tissue, species, additional_info):
    return run_cell_type_analysis_unified(model, temperature, marker_list, tissue, species, additional_info, provider="anthropic")

def run_cell_type_analysis_openrouter(model, temperature, marker_list, tissue, species, additional_info):
    return run_cell_type_analysis_unified(model, temperature, marker_list, tissue, species, additional_info, provider="openrouter")