

# Update Note

Next Step:

1. Marker Grouping

~~2. Incorporating more LLM (CLAUDE)~~

3. Incorporating RAG

~~4. Build Web UI~~

5. Build Python package, R package

6. Incorporating subcluster and cluster funcation



# Version Overview

This project has multiple versions with varying levels of complexity and automation:

Version 1 (Full Feature Set)

The most comprehensive version, including:

- Onboarding Agent

- Functional Gene Expert Agent

- Cell Type Gene Expert Agent

- Integration Agent

- Validator Agent

- Human Input Component

- Final Formatting Agent

Version 2 (Reduced Expert Agents)

Similar to Version 1, but excludes:

- Functional Gene Expert Agent

- Cell Type Gene Expert Agent

Version 3 (Automated Process)

Fully automatic version, excluding:

- Onboarding Agent

- Human Input Agent

Version 4 (Minimal Automatic)

The most streamlined version, including only:

- Integration Agent















# Table of Contents

1. Introduction

2. Prerequisites

3. Script Overview

4. Key Components

   4.1 Agent Class

   4.2 Conversation Functions

   4.3 Data Processing Functions

   4.4 Specialized Agents

5. Main Execution Flow

6. Output

7. Customization and Extension

8. Troubleshooting

9. Conclusion

1. Introduction

This script implements a multi-agent system for analyzing single-cell data. It uses a series of AI agents to perform functional analysis, cell type identification, and integrative annotation of gene markers.

2. Prerequisites

- Python 3.x

- Required libraries: `openai`, `dotenv`, `httpx`, `re`, `os`, `json`

- OpenAI API key (set in environment variables)

3. Script Overview

The script performs the following high-level steps:

1. Onboarding process to gather initial information

2. Functional analysis of gene markers

3. Cell type analysis of gene markers

4. Integrative analysis combining results from steps 2 and 3

5. Formatting and structuring the final output

4. Key Components

4.1 Agent Class

The `Agent` class is the core component of the script. It represents an AI agent capable of engaging in conversations and performing specific tasks.

Key features:

- Initialization with system prompt, model, and interaction settings

- Conversation history management

- Execution of AI model queries

- Human input handling

4.2 Conversation Functions

- `two_agent_conversation_with_validation_celltype`: Manages the conversation for cell type analysis

- `two_agent_conversation_with_validation_functional`: Manages the conversation for functional analysis

- `integrate_and_annotate`: Handles the integrative analysis conversation

- `onboarding_process`: Gathers initial information from the user

4.3 Data Processing Functions

- `extract_agent1_analysis`: Extracts analysis results from conversation history

- `extract_json_from_reply`: Parses JSON data from agent responses

- `list_to_comma_separated_string`: Converts list to comma-separated string

- `construct_prompt`: Builds the initial prompt based on user input

4.4 Specialized Agents

The script defines several specialized agents:

- `agent1_functional` and `agent2_functional`: For functional analysis

- `agent1_celltype` and `agent2_celltype`: For cell type analysis

- `agent3`: For onboarding

- `integrative_agent`: For integrative analysis

- `formatting_agent`: For formatting final results

5. Main Execution Flow

1. Start onboarding process

2. Conduct functional analysis

3. Perform cell type analysis

4. Execute integrative analysis

5. Format and structure final results

6. Output

The script produces a structured JSON output containing:

- Main cell type identified

- Sub-cell types (if applicable)

7. Customization and Extension

- Modify agent system prompts to adjust behavior

- Add new specialized agents for additional analyses

- Extend the `Agent` class for more complex interactions

8. Troubleshooting

- Ensure all required libraries are installed

- Verify that the OpenAI API key is correctly set in environment variables

- Check for any rate limiting or API usage issues with OpenAI

9. Conclusion

This script provides a flexible framework for single-cell data analysis using AI agents. It can be adapted and extended for various types of biological data analysis tasks.
