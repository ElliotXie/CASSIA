import json
import re
from openai import OpenAI
import os
import anthropic
import requests
from typing import Optional, Dict, Any

# Import call_llm for the shared Agent class
try:
    from ..core.llm_utils import call_llm
except ImportError:
    from core.llm_utils import call_llm


# ----------------- Shared Agent Class -----------------

class Agent:
    """
    Unified Agent class that works with all LLM providers via call_llm.

    Maintains conversation history and supports reasoning tokens for models that support it.
    """

    def __init__(
        self,
        system: str = "",
        model: str = None,
        temperature: float = 0.7,
        provider: str = "openrouter",
        reasoning: Optional[Dict[str, Any]] = None,
        api_key: Optional[str] = None,
        base_url: Optional[str] = None  # For custom OpenAI-compatible endpoints
    ):
        self.system = system
        self.model = model
        self.temperature = temperature
        self.provider = provider
        self.reasoning = reasoning
        self.api_key = api_key
        self.base_url = base_url
        self.chat_histories: Dict[str, list] = {}

    def __call__(self, message: str, other_agent_id: str) -> str:
        """
        Send a message and get a response, maintaining conversation history.

        Args:
            message: The message to send
            other_agent_id: Identifier for the conversation (allows multiple parallel conversations)

        Returns:
            The assistant's response text
        """
        # Initialize history for this conversation
        if other_agent_id not in self.chat_histories:
            self.chat_histories[other_agent_id] = []
            if self.system:
                self.chat_histories[other_agent_id].append(
                    {"role": "system", "content": self.system}
                )

        # Add user message to history
        self.chat_histories[other_agent_id].append(
            {"role": "user", "content": message}
        )

        # Determine the provider string for call_llm
        # For custom endpoints, use the base_url as provider
        effective_provider = self.base_url if self.base_url else self.provider

        # Call LLM using call_llm with conversation history
        result = call_llm(
            prompt=message,
            provider=effective_provider,
            model=self.model,
            temperature=self.temperature,
            system_prompt=self.system,
            reasoning=self.reasoning,
            api_key=self.api_key,
            additional_params={"messages": self.chat_histories[other_agent_id]}
        )

        # Add assistant response to history
        self.chat_histories[other_agent_id].append(
            {"role": "assistant", "content": result}
        )

        return result


# ----------------- Helper Functions and Prompts (Defined Once) -----------------

def extract_json_from_reply(reply):
    """Extracts a JSON object from a string, typically the LLM's reply."""
    json_match = re.search(r'```json\n(.*?)\n```', reply, re.DOTALL)
    if json_match:
        json_str = json_match.group(1)
        try:
            # Use strict=False to allow control characters (e.g., newlines) within strings
            return json.loads(json_str, strict=False)
        except json.JSONDecodeError as e:
            print(f"Error decoding JSON even with non-strict parsing: {e}")
    else:
        print("No JSON content found in the reply")
    return None

def construct_prompt(json_data):
    """Constructs the initial prompt for the cell type annotation task."""
    marker_list = ', '.join(json_data['marker_list'])
    prompt = f"Your task is to annotate a single-cell {json_data['species']} dataset"
    
    if json_data.get('tissue_type') and json_data['tissue_type'].lower() not in ['none', 'tissue blind']:
        prompt += f" from {json_data['tissue_type']} tissue"
    
    prompt += f". Please identify the cell type based on this ranked marker list:\n{marker_list}"
    
    if json_data.get('additional_info') and json_data['additional_info'].lower() != "no":
        prompt += f" Below is some additional information about the dataset:\n{json_data['additional_info']}."
    return prompt

def final_annotation(agent, prompt):
    """Manages the conversation with the final annotation agent."""
    conversation = []
    while True:
        response = agent(prompt, "user")
        conversation.append(("Final Annotation Agent", response))
        if "FINAL ANNOTATION COMPLETED" in response:
            break
        prompt = response
    return conversation

def coupling_validation(agent, annotation_result, onboarding_data):
    """Constructs and sends the validation prompt to the coupling validator agent."""
    
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
    
    validation_message = f"""Please validate the following annotation result:

Annotation Result:
{annotation_result}

Context:

Marker List: {marker_list_str}
Additional Info: {onboarding_data.get('additional_info', 'None')}

Validate the annotation based on this context.
"""
    return agent(validation_message, "final_annotation")

def format_results(agent, final_annotations):
    """Calls the formatting agent and returns its raw response."""
    final_text = "\n\n".join([msg[1] for msg in final_annotations])
    return agent(final_text, "user")

# ----------------- System Prompts (Defined Once) -----------------

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

coupling_validator_system_v0 = """
You are an expert biologist specializing in single-cell analysis. Your critical role is to
validate the final annotation results for a cell cluster. You will be provided with the proposed
annotation result, and a ranked list of marker genes it used. Below are steps to follow:
1. Marker Consistency: Make sure the markers are in the provided marker list. Make sure
there is consistency between the identified cell type and the provided markers.
2. Mixed Cell Type Consideration: Be aware that mixed cell types may be present. Only raise
this point if multiple distinct cell types are strongly supported by several high-ranking
markers. In cases of potential mixed populations, flag this for further investigation rather
than outright rejection.
Output Format: if passed, Validation result: VALIDATION PASSED. If failed, Validation
result: VALIDATION FAILED. Feedback: give detailed feedback and instructions for revising
the annotation.
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
5. IMPORTANT: Ensure that all string values in the JSON are properly escaped. For example, any newline characters inside a string must be represented as `\\\\n`.

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
4. IMPORTANT: Ensure that all string values in the JSON are properly escaped. For example, any newline characters inside a string must be represented as `\\\\n`.

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

# ----------------- Unified Analysis Logic -----------------

def _run_analysis_logic(final_annotation_agent, coupling_validator_agent, formatting_agent, user_data, is_tissue_blind):
    """
    Core analysis workflow. This function is called by the provider-specific wrappers.
    """
    prompt = construct_prompt(user_data)

    validation_passed = False
    iteration = 0
    max_iterations = 3
    full_conversation_history = []
    final_annotation_conversation = []

    while not validation_passed and iteration < max_iterations:
        iteration += 1
        
        current_prompt = prompt
        if iteration > 1:
            current_prompt = f"""Previous annotation attempt failed validation. Please review your previous response and the validation feedback, then provide an updated annotation:

Previous response:
{final_annotation_conversation[-1][1]}

Validation feedback:
{validation_result}

Original prompt:
{prompt}

Please provide an updated annotation addressing the validation feedback."""

        final_annotation_conversation = final_annotation(final_annotation_agent, current_prompt)
        full_conversation_history.extend(final_annotation_conversation)
        
        validation_result = coupling_validation(coupling_validator_agent, final_annotation_conversation[-1][1], user_data)
        full_conversation_history.append(("Coupling Validator", validation_result))
        
        if "VALIDATION PASSED" in validation_result:
            validation_passed = True

    # Determine which formatting system to use based on validation success and tissue context
    if validation_passed:
        formatting_agent.system = formatting_system_tissue_blind if is_tissue_blind else formatting_system_non_tissue_blind
    else:
        formatting_agent.system = formatting_system_failed
        
    # Format and return results
    raw_formatted_output = format_results(formatting_agent, final_annotation_conversation[-2:])
    full_conversation_history.append(("Formatting Agent", raw_formatted_output))
    
    structured_output = extract_json_from_reply(raw_formatted_output)
    
    if structured_output:
        structured_output["iterations"] = iteration
        structured_output["num_markers"] = len(user_data['marker_list'])
        return structured_output, full_conversation_history
    else:
        print("Error: Unable to extract JSON from the formatted output.")
        print("Raw formatted output:")
        print(raw_formatted_output)

    # Fallback if JSON parsing fails
    return None, full_conversation_history

# ----------------- Public API Functions (Wrappers) -----------------

def run_cell_type_analysis(model, temperature, marker_list, tissue, species, additional_info, validator_involvement="v1", reasoning=None):
    """OpenAI provider for cell type analysis."""
    is_tissue_blind = tissue.lower() in ['none', 'tissue blind'] if tissue else True
    final_annotation_agent = Agent(
        system=final_annotation_system_v2 if is_tissue_blind else final_annotation_system_v1,
        model=model,
        temperature=temperature,
        provider="openai",
        reasoning=reasoning
    )
    
    # Select validator system based on involvement level
    if validator_involvement == "v0":
        validator_system = coupling_validator_system_v0.strip()
    elif validator_involvement == "v1":
        validator_system = coupling_validator_system_v2.strip() if is_tissue_blind else coupling_validator_system_v1.strip()
    else:
        validator_system = coupling_validator_system_v2.strip() if is_tissue_blind else coupling_validator_system_v1.strip()
    
    coupling_validator_agent = Agent(
        system=validator_system,
        model=model,
        temperature=temperature,
        provider="openai",
        reasoning=reasoning
    )
    formatting_agent = Agent(
        system="",  # System prompt is set inside the logic function
        model=model,
        temperature=temperature,
        provider="openai",
        reasoning=reasoning
    )

    user_data = {"species": species, "tissue_type": tissue, "marker_list": marker_list}
    if additional_info and additional_info.lower() != "no":
        user_data["additional_info"] = additional_info

    return _run_analysis_logic(final_annotation_agent, coupling_validator_agent, formatting_agent, user_data, is_tissue_blind)


def run_cell_type_analysis_claude(model, temperature, marker_list, tissue, species, additional_info, validator_involvement="v1", reasoning=None):
    """Anthropic Claude provider for cell type analysis."""
    is_tissue_blind = tissue.lower() in ['none', 'tissue blind'] if tissue else True

    final_annotation_agent = Agent(
        system=final_annotation_system_v2 if is_tissue_blind else final_annotation_system_v1,
        model=model,
        temperature=temperature,
        provider="anthropic",
        reasoning=reasoning
    )

    # Select validator system based on involvement level
    if validator_involvement == "v0":
        validator_system = coupling_validator_system_v0.strip()
    elif validator_involvement == "v1":
        validator_system = coupling_validator_system_v2.strip() if is_tissue_blind else coupling_validator_system_v1.strip()
    else:
        validator_system = coupling_validator_system_v2.strip() if is_tissue_blind else coupling_validator_system_v1.strip()

    coupling_validator_agent = Agent(
        system=validator_system,
        model=model,
        temperature=temperature,
        provider="anthropic",
        reasoning=reasoning
    )
    formatting_agent = Agent(
        system="",
        model=model,
        temperature=temperature,
        provider="anthropic",
        reasoning=reasoning
    )

    user_data = {"species": species, "tissue_type": tissue, "marker_list": marker_list}
    if additional_info and additional_info.lower() != "no":
        user_data["additional_info"] = additional_info

    return _run_analysis_logic(final_annotation_agent, coupling_validator_agent, formatting_agent, user_data, is_tissue_blind)

def run_cell_type_analysis_openrouter(model, temperature, marker_list, tissue, species, additional_info, validator_involvement="v1", reasoning=None):
    """OpenRouter provider for cell type analysis. Supports reasoning tokens for compatible models."""
    is_tissue_blind = tissue.lower() in ['none', 'tissue blind'] if tissue else True

    final_annotation_agent = Agent(
        system=final_annotation_system_v2 if is_tissue_blind else final_annotation_system_v1,
        model=model,
        temperature=temperature,
        provider="openrouter",
        reasoning=reasoning
    )

    # Select validator system based on involvement level
    if validator_involvement == "v0":
        validator_system = coupling_validator_system_v0.strip()
    elif validator_involvement == "v1":
        validator_system = coupling_validator_system_v2.strip() if is_tissue_blind else coupling_validator_system_v1.strip()
    else:
        validator_system = coupling_validator_system_v2.strip() if is_tissue_blind else coupling_validator_system_v1.strip()

    coupling_validator_agent = Agent(
        system=validator_system,
        model=model,
        temperature=temperature,
        provider="openrouter",
        reasoning=reasoning
    )
    formatting_agent = Agent(
        system="",
        model=model,
        temperature=temperature,
        provider="openrouter",
        reasoning=reasoning
    )

    user_data = {"species": species, "tissue_type": tissue, "marker_list": marker_list}
    if additional_info and additional_info.lower() != "no":
        user_data["additional_info"] = additional_info

    return _run_analysis_logic(final_annotation_agent, coupling_validator_agent, formatting_agent, user_data, is_tissue_blind)

def run_cell_type_analysis_custom(base_url, api_key, model, temperature, marker_list, tissue, species, additional_info, validator_involvement="v1", reasoning=None):
    """
    Custom API provider for OpenAI-compatible endpoints.
    """
    is_tissue_blind = tissue.lower() in ['none', 'tissue blind'] if tissue else True

    final_annotation_agent = Agent(
        system=final_annotation_system_v2 if is_tissue_blind else final_annotation_system_v1,
        model=model,
        temperature=temperature,
        base_url=base_url,  # Use base_url for custom endpoints
        api_key=api_key,
        reasoning=reasoning
    )

    # Select validator system based on involvement level
    if validator_involvement == "v0":
        validator_system = coupling_validator_system_v0.strip()
    elif validator_involvement == "v1":
        validator_system = coupling_validator_system_v2.strip() if is_tissue_blind else coupling_validator_system_v1.strip()
    else:
        validator_system = coupling_validator_system_v2.strip() if is_tissue_blind else coupling_validator_system_v1.strip()

    coupling_validator_agent = Agent(
        system=validator_system,
        model=model,
        temperature=temperature,
        base_url=base_url,
        api_key=api_key,
        reasoning=reasoning
    )
    formatting_agent = Agent(
        system="",
        model=model,
        temperature=temperature,
        base_url=base_url,
        api_key=api_key,
        reasoning=reasoning
    )

    user_data = {"species": species, "tissue_type": tissue, "marker_list": marker_list}
    if additional_info and additional_info.lower() != "no":
        user_data["additional_info"] = additional_info

    return _run_analysis_logic(final_annotation_agent, coupling_validator_agent, formatting_agent, user_data, is_tissue_blind)