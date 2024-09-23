import re
from dotenv import load_dotenv
import openai
import re
import httpx
import os
import json
from dotenv import load_dotenv
import anthropic

_ = load_dotenv()
from openai import OpenAI
client = OpenAI()
_ = load_dotenv()
from openai import OpenAI
client = OpenAI()

class Agent:
    def __init__(self, system="", human_input_mode="never", human_input_keyword="need human input", 
                 ai_provider="openai", model="gpt-4", temperature=0, max_tokens=1000):
        self.system = system
        self.chat_histories = {}
        self.human_input_mode = human_input_mode
        self.human_input_keyword = human_input_keyword.lower()
        self.ai_provider = ai_provider.lower()
        self.model = model
        self.temperature = temperature
        self.max_tokens = max_tokens

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
        if self.ai_provider == "anthropic":
            # Use Claude AI
            client = anthropic.Anthropic()
            
            # Prepare messages for Claude (excluding system message)
            claude_messages = []
            for msg in self.chat_histories[other_agent_id]:
                if msg["role"] != "system":
                    claude_messages.append({
                        "role": msg["role"],
                        "content": [{"type": "text", "text": msg["content"]}]
                    })
            
            response = client.messages.create(
                model=self.model,
                max_tokens=self.max_tokens,
                temperature=self.temperature,
                system=self.system,  # System message as a separate parameter
                messages=claude_messages
            )
            return response.content[0].text
        else:
            # Use OpenAI
            # For OpenAI, we keep the system message in the chat history
            completion = client.chat.completions.create(
                model=self.model,
                temperature=self.temperature,
                max_tokens=self.max_tokens,
                messages=self.chat_histories[other_agent_id]
            )
            return completion.choices[0].message.content

    def needs_human_input(self, message):
        if self.human_input_mode == "always":
            return True
        elif self.human_input_mode == "keyword" and self.human_input_keyword in message.lower():
            return True
        return False

    def wants_to_terminate(self, message):
        return "TERMINAL" in message.upper()

def two_agent_conversation_with_validation_celltype(agent1, agent2, initial_message, marker_list):
    conversation_history = []
    current_message = f"{initial_message}\nMarker list: {marker_list}"
    conversation_history.append(("Initial", current_message))
    print(f"Initial: {current_message}\n")
    
    agent1_responses = []
    while True:
        # Agent 1's analysis loop
        while True:
            response1 = agent1(current_message, "user")
            conversation_history.append(("Agent 1", response1))
            agent1_responses.append(response1)
            print(f"Agent 1: {response1}\n", flush=True)

            if "CELL TYPE ANALYSIS COMPLETED" in response1:
                break
            
            if agent1.needs_human_input(response1):
                human_input = input("User input: ")
                conversation_history.append(("User", human_input))
                print(f"User: {human_input}\n", flush=True)
                current_message = human_input
            else:
                current_message = response1

        # Agent 2's validation
        last_n_responses = "\n".join(agent1_responses[-3:])  # Send the last 5 responses from Agent 1
        validation_message = f"Marker list: {marker_list}\n\nAgent 1's analysis:\n{last_n_responses}"
        response2 = agent2(validation_message, "agent1")
        conversation_history.append(("Agent 2", response2))
        print(f"Agent 2: {response2}\n")
        
        if "VALIDATION COMPLETED" in response2:
            break
        elif "NEED MORE WORK" in response2:
            current_message = f"Agent 2 feedback: {response2}\nPlease revise your analysis based on this feedback. Marker list: {marker_list}"
            # Immediately get and print Agent 1's response to the feedback
            response1 = agent1(current_message, "agent2")
            conversation_history.append(("Agent 1", response1))
            print(f"Agent 1: {response1}\n", flush=True)
            current_message = response1
        else:
            print("Unexpected response from Agent 2. Ending conversation.")
            break
    
    return conversation_history, agent1_responses

def two_agent_conversation_with_validation_functional(agent1, agent2, initial_message, marker_list):
    conversation_history = []
    current_message = f"{initial_message}\nMarker list: {marker_list}"
    conversation_history.append(("Initial", current_message))
    print(f"Initial: {current_message}\n")
    
    agent1_responses = []
    while True:
        # Agent 1's analysis loop
        while True:
            response1 = agent1(current_message, "user")
            conversation_history.append(("Agent 1", response1))
            agent1_responses.append(response1)
            print(f"Agent 1: {response1}\n", flush=True)

            if "FUNCTIONAL ANALYSIS COMPLETED" in response1:
                break
            
            if agent1.needs_human_input(response1):
                human_input = input("User input: ")
                conversation_history.append(("User", human_input))
                print(f"User: {human_input}\n", flush=True)
                current_message = human_input
            else:
                current_message = response1

        # Agent 2's validation
        last_n_responses = "\n".join(agent1_responses[-3:])  # Send the last 5 responses from Agent 1
        validation_message = f"Marker list: {marker_list}\n\nAgent 1's analysis:\n{last_n_responses}"
        response2 = agent2(validation_message, "agent1")
        conversation_history.append(("Agent 2", response2))
        print(f"Agent 2: {response2}\n")
        
        if "VALIDATION COMPLETED" in response2:
            break
        elif "NEED MORE WORK" in response2:
            current_message = f"Agent 2 feedback: {response2}\nPlease revise your analysis based on this feedback. Marker list: {marker_list}"
            # Immediately get and print Agent 1's response to the feedback
            response1 = agent1(current_message, "agent2")
            conversation_history.append(("Agent 1", response1))
            print(f"Agent 1: {response1}\n", flush=True)
            current_message = response1
        else:
            print("Unexpected response from Agent 2. Ending conversation.")
            break
    
    return conversation_history, agent1_responses

def generate_instruction_prompt(prompt, marker_list, functional_analysis, celltype_analysis):
    instruction = f"""
    {prompt}
    
    Below is the Functional Analysis results:
    {functional_analysis}

    Below is the Cell Type Analysis results:
    {celltype_analysis}
    """.strip()
    return instruction



def integrate_and_annotate(agent, instruction):
    print(f"Instruction: {instruction}\n")
    current_message = instruction
    print(f"Instruction: {current_message}\n")
    conversation_history = []
    
    while True:
        response = agent(current_message, "user")
        print(f"Agent: {response}\n", flush=True)
        conversation_history.append(("Agent", response))
        
        if "FINAL ANNOTATION COMPLETED" in response:
            break
        
        if agent.needs_human_input(response):
            user_input = input("User input: ")
            print(f"User: {user_input}\n", flush=True)
            conversation_history.append(("User", user_input))
            current_message = user_input
        else:
            current_message = response

    return conversation_history

def extract_agent1_analysis(conversation_history):
    analysis_results = []
    current_analysis = ""
    analyzing = False

    for role, message in conversation_history:
        if role == "Agent 1":
            if "ANALYSIS BEGIN" in message:
                analyzing = True
                current_analysis = message
            elif analyzing:
                current_analysis += "\n" + message
            if "ANALYSIS COMPLETED" in message:
                analyzing = False
                analysis_results.append(current_analysis)
                current_analysis = ""

    return analysis_results

def extract_group_content(text):
    pattern = r'### Group \d+:.*?(?=### Group \d+:|### Feedback Request|$)'
    matches = re.findall(pattern, text, re.DOTALL)
    cleaned_matches = []
    for match in matches:
        cleaned = '\n'.join(line.strip() for line in match.split('\n') if line.strip())
        cleaned_matches.append(cleaned)
    return cleaned_matches

def onboarding_process(agent3):
    current_message = "I am your helpful assistant. I will ask you some questions about your single cell analysis so that your analysis can be more accurate."
    conversation = []
    user_data = {}
    
    while True:
        response = agent3(current_message, "user")
        print(f"Onboarding Assistant: {response}\n", flush=True)
        conversation.append(("Onboarding Assistant", response))
        
        if "ONBOARDING COMPLETED" in response:
            break
        
        user_input = input("User: ")
        print(f"User: {user_input}\n", flush=True)
        conversation.append(("User", user_input))
        current_message = user_input
        
    return conversation


def extract_json_from_reply(reply):
    # Remove the extra backslashes, quotes, and backticks
    cleaned_reply = reply.strip("'").replace("\\n", "\n").replace('\\"', '"').replace("```", "")
    
    # Use regex to find the JSON content
    json_match = re.search(r'json\s*(\{.*?\})', cleaned_reply, re.DOTALL)
    
    if json_match:
        json_str = json_match.group(1)
        try:
            json_data = json.loads(json_str)
            return json_data
        except json.JSONDecodeError as e:
            print(f"Error decoding JSON: {e}")
            print(f"JSON string: {json_str}")
            return None
    else:
        print("No JSON content found in the reply")
        print(f"Cleaned reply: {cleaned_reply}")
        return None

def list_to_comma_separated_string(input_list):
    return ', '.join(input_list)

def construct_prompt(json_data):
    species = json_data['species']
    tissue = json_data['tissue_type']
    additional_info = json_data.get('additional_info', '')  # Use .get() with a default value
    marker_list = json_data['marker_list']

    initial_message = f"I am analyzing a single-cell {species} {tissue} dataset."
    if additional_info:
        initial_message += f" {additional_info}."
    initial_message += " I want to identify the cell types present based on a marker list:"

    marker_string = ', '.join(marker_list)

    full_prompt = f"{initial_message}\nBelow is the marker list shows highly expressed genes ranked by average expression from highest to lowest:\n{marker_string}"

    return full_prompt

agent1_functional = Agent(system=
'''You are a functional gene professional. Given the marker list.
Please organize these markers into comprehensive functional gene groups. Include as many group as you can. One extra group will give you 1000$ bonus. 
Start your analysis saying "ANLAYSIS BEGIN round n":
For each group:

1. Provide a name for the functional group
2. List the genes from the input that belong to this group.(if you found series of gene have similar name pattern in one group 
        such as Krt1,Krt2,Krt3.. or Cldn1,Cldn2,Cldn3 you are not allowed to list all the name separately, just summerize
        them as Krt family or Cldn family, thats it)
3. Briefly explain the relevance of this group to cell type identification
4. If there are any genes in the list that you're unsure about or that don't clearly fit into a functional group, 
    please list them separately and suggest possible groupings or areas for further investigation.
5. say "### Feedback Request" and then Ask for human input for feedback on your analysis.

Only after you see human input, and human says good. say "FUNCTIONAL ANALYSIS COMPLETED".
 
'''.strip(),
    human_input_mode="always",
    ai_provider="anthropic",
    model="claude-3-5-sonnet-20240620")

agent2_functional = Agent(system='''You are a very careful biologist. 
               Your task is to validate the functional gene enriched groups using the marker list.
                Please perform the following tasks:

        1. Make sure to get ALL listed groups from previous chat.
               
        2. Check if each gene exists in the provided marker list. Identify any genes that were incorrectly included (not present in the marker list).

        4. Say "NEED MORE WORK" if any genes are not in the marker list.
           Say "VALIDATION COMPLETED" otherwise.

'''.strip(),     human_input_mode="never",
    ai_provider="anthropic",
    model="claude-3-5-sonnet-20240620")

agent3 = Agent(system='''You are a friendly onboarding assistant helping to create an initial prompt for analyzing single-cell mouse larynx data.
Ask the user questions to gather necessary information. Be concise and friendly in your interactions. 
               Ask one question at a time. If the user provide confusing answer, friendly ask user to provide clear answer.

1. What is the species of the data?
2. What is the tissue type of the data?
3. Please provide the marker list for the analysis.
4. Any other information you need to know?

When you're done gathering information, say "ONBOARDING COMPLETED" followed by a JSON-formatted summary of the collected information.Name each item by species, tissue_type, marker_list, additional_information
wrap it with ```json and ```
'''.strip(), human_input_mode="always",
    ai_provider="anthropic",
    model="claude-3-5-sonnet-20240620")

agent1_celltype = Agent(system=
'''
        You are a cell type identification professional. Given a list of gene markers, 
        please organize these markers into celltype groups. Try to include as many groups as you can. A huge 10000$ bonus will be provided if you do a good job.
        Start your analysis saying "ANALYSIS BEGIN round n":
        For each group:

        1. Provide a name for the cell type or cell type-related group
        2. List the genes from the input that belong to this group (if you found series of gene have similar name pattern in one group 
        such as Krt1,Krt2,Krt3.. or Cldn1,Cldn2,Cldn3 you are not allowed to list all the name separately, just summerize
        them as Krt family or Cldn family, thats it)
        3. Briefly explain the relevance of this group to cell type identification
        4. If there are any genes in the list that you're unsure about or that don't clearly fit into a cell type group, 
           please list them separately and suggest possible groupings or areas for further investigation.
        5. say "### Feedback Request" and then Ask for human input for feedback on your analysis.

        Only after you see human input, and human says good, say "CELL TYPE ANALYSIS COMPLETED".
 
''',
    human_input_mode="always",
    ai_provider="anthropic",
    model="claude-3-5-sonnet-20240620")

agent2_celltype = Agent(system=''' You are a very careful cell biologist specializing in cell type identification. 
        Your task is to validate all cell type-specific gene groups mentioned using the provided marker list.
        Please perform the following tasks:
        1. Make sure to get ALL listed groups from agent 1.
               
        2. Check if each gene exists in the provided marker list.

        3. For each cell type group:
           - Identify any genes that were incorrectly included (not present in the marker list).

        4. Say "NEED MORE WORK" if any genes are not in the marker list.
           Say "VALIDATION COMPLETED" otherwise.
               
''', 
    human_input_mode="never",
    ai_provider="anthropic",
    model="claude-3-5-sonnet-20240620")




integrative_agent = Agent(system="""
You are an expert cell biologist specializing in specializing in cell type annotation. 
Your task is to make use a list of highly expressed marker, functional analysis results
and cell type analyses results from expert to provide a final cell type annotation step by step. You will received
99999$ bonus if you do a good job. Remember to respect the tissue type and species, only make reasonable annotation.
                          
Please do:
1. Anlysis the functional gene group:include your own reasoning and analyze the provided analysis results.
2. Anlysis the cell type group:include your own reasoning and analyze the provided analysis results.
3. **Cross-reference Known Databases**: Use available scRNA-seq databases and relevant literature to cross-reference these markers. list your finding.
4. **Determine the Most Probable General Cell Type**: Based on the expression of these markers, infer the most likely general cell type of the cluster.
5. **Identify the Top 3 Most Probable Sub Cell Types**: Based on the expression of these markers, infer the top three most probable sub cell types within the general cell type. Finally, specify the most likely subtype.
6. **Identify the Most Probable Sub-Sub Cell Type**: Determine the most specific cell type within the previously identified subtype.
7.  **Provide a Concise Summary of Your Analysis

The provided analysis results are only for reference. Always include your step by step detailed reasoning.                      
Ask for human input.
You can say "FINAL ANNOTATION COMPLETED" only when human say good.
                          

""".strip(), human_input_mode="always",
    ai_provider="anthropic",
    model="claude-3-5-sonnet-20240620")


formatting_agent = Agent(system="""
You are a formatting assistant for single-cell analysis results. Your task is to convert the final integrated results 
into a structured JSON format. Follow these guidelines:

1. Extract the main cell type and any sub-cell types identified.
2. Include only information explicitly stated in the input.
3. Ensure the output is valid JSON.

Provide the JSON output within triple backticks, like this:
'''json
{
"main_cell_type": "...",
"sub_cell_types": ["...", "..."]
}
'''

""", human_input_mode="never",
    ai_provider="anthropic",
    model="claude-3-5-sonnet-20240620")

def format_results(agent, final_annotations):
    # Convert the final annotations to a single string
    final_text = "\n\n".join([msg[1] for msg in final_annotations])
    return agent(final_text, "user")


def extract_analysis_content(text):
    start_marker = "ANALYSIS BEGIN"
    end_marker = "### Feedback Request"
    
    start_index = text.find(start_marker)
    end_index = text.find(end_marker)
    
    if start_index != -1 and end_index != -1:
        extracted_content = text[start_index:end_index].strip()
        return extracted_content
    else:
        return "Analysis content not found or improperly formatted."


# Main execution
if __name__ == "__main__":
    print("Starting onboarding process...\n")
    onboarding_conversation = onboarding_process(agent3)
    extracted_data = extract_json_from_reply(onboarding_conversation[-1][1])
    marker_list = extracted_data['marker_list']
    marker_list = list_to_comma_separated_string(extracted_data['marker_list'])
    prompt = construct_prompt(extracted_data)

    print("Starting functional analysis...\n")
    conversation_history_functional, agent1_responses_functional = two_agent_conversation_with_validation_functional(agent1_functional, agent2_functional, prompt, marker_list)
    analysis_results_functional = extract_agent1_analysis(conversation_history_functional)

    print("Starting cell type analysis...\n")
    conversation_history_celltype, agent1_responses_celltype = two_agent_conversation_with_validation_celltype(agent1_celltype, agent2_celltype, prompt, marker_list)
    analysis_results_celltype = extract_agent1_analysis(conversation_history_celltype)

    print("Starting integrative analysis...\n")
    functional_content = extract_analysis_content(analysis_results_functional[-1])
    celltype_content = extract_analysis_content(analysis_results_celltype[-1])
    instruction_prompt = generate_instruction_prompt(prompt, marker_list, functional_content, celltype_content)
    
    final_annotation_conversation = integrate_and_annotate(integrative_agent, instruction_prompt)

    print("Formatting final results...")
    print(final_annotation_conversation[-1:])

    # Get the last two replies from the integrative agent
    final_annotations = final_annotation_conversation[-3:]
    formatted_output = format_results(formatting_agent, final_annotations)
    
    structured_output = extract_json_from_reply(formatted_output)
    if structured_output:
        print("\nStructured output:")
        print(json.dumps(structured_output, indent=2))
    else:
        print("Error: Unable to extract JSON from the formatted output.")
        print("Raw formatted output:")
        print(formatted_output)

    print("Analysis complete.")