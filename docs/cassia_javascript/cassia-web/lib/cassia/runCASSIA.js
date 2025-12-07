import { callLLM } from './llm_utils.js';

// ----------------- System Prompts -----------------
// These prompts are faithfully ported from Python main_function_code.py

const finalAnnotationSystemV1 = `
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
`.trim();

const finalAnnotationSystemV2 = `
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
`.trim();

const couplingValidatorSystemV0 = `
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
`.trim();

const couplingValidatorSystemV1 = `
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
`.trim();

const couplingValidatorSystemV2 = `
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
`.trim();

// Formatting system for tissue-blind mode (when tissue is not specified)
const formattingSystemTissueBlind = `
You are a formatting assistant for single-cell analysis results. Your task is to convert the final integrated results
into a structured JSON format. Follow these guidelines:

1. Extract the main cell type and any sub-cell types identified.
2. Include only information explicitly stated in the input.
3. If there are possible mixed cell types highlighted, list them.
4. Include possible tissues.
5. IMPORTANT: Ensure that all string values in the JSON are properly escaped. For example, any newline characters inside a string must be represented as \\\\n.

Provide the JSON output within triple backticks, like this:
\`\`\`json
{
"main_cell_type": "...",
"sub_cell_types": ["...", "..."],
"possible_mixed_cell_types": ["...", "..."],
"possible_tissues": ["...", "..."]
}
\`\`\`
`.trim();

// Formatting system for non-tissue-blind mode (when tissue is specified)
const formattingSystemNonTissueBlind = `
You are a formatting assistant for single-cell analysis results. Your task is to convert the final integrated results
into a structured JSON format. Follow these guidelines:

1. Extract the main cell type and the three most likely sub-cell types identified from step 4 and step 5 of the Final Annotation Agent response. Even the main cell type is the same as the sub-cell types, you still need to list it as a sub-cell type. Strictly follow the order of the sub-cell types.
2. Include only information explicitly stated in the input.
3. If there are possible mixed cell types highlighted, list them.
4. IMPORTANT: Ensure that all string values in the JSON are properly escaped. For example, any newline characters inside a string must be represented as \\\\n.

Provide the JSON output within triple backticks, like this:
\`\`\`json
{
"main_cell_type": "...",
"sub_cell_types": ["...", "..."],
"possible_mixed_cell_types": ["...", "..."]
}
\`\`\`
`.trim();

// Formatting system for failed analyses (when validation fails after max iterations)
const formattingSystemFailed = `
You are a formatting assistant for single-cell analysis results. Your task is to convert the final integrated results
into a structured JSON format, with special consideration for uncertain or conflicting annotations. Follow these guidelines:

1. The analysis failed after multiple attempts. Please try to extract as much information as possible. Summarize what has gone wrong and what has been tried.
2. Provide a detailed feedback on why the analysis failed, and what has been tried and why it did not work.
3. Finally, provide a detailed step-by-step reasoning of how to fix the analysis.

Provide the JSON output within triple backticks, like this:
\`\`\`json
{
"main_cell_type": "if any",
"sub_cell_types": "if any",
"possible_cell_types": "if any",
"feedback": "...",
"next_steps": "..."
}
\`\`\`
`.trim();

// ----------------- Helper Functions -----------------

function extractJsonFromReply(reply) {
    const jsonMatch = reply.match(/```json\n(.*?)\n```/s);
    if (jsonMatch) {
        const jsonStr = jsonMatch[1];
        try {
            return JSON.parse(jsonStr);
        } catch (e) {
            console.error(`Error decoding JSON: ${e}`);
        }
    } else {
        console.log("No JSON content found in the reply");
    }
    return null;
}

function constructPrompt(jsonData) {
    const markerList = Array.isArray(jsonData.marker_list) 
        ? jsonData.marker_list.join(', ')
        : jsonData.marker_list;
    
    let prompt = `Your task is to annotate a single-cell ${jsonData.species} dataset`;
    
    if (jsonData.tissue_type && !['none', 'tissue blind'].includes(jsonData.tissue_type.toLowerCase())) {
        prompt += ` from ${jsonData.tissue_type} tissue`;
    }
    
    prompt += `. Please identify the cell type based on this ranked marker list:\n${markerList}`;
    
    if (jsonData.additional_info && jsonData.additional_info.toLowerCase() !== "no") {
        prompt += ` Below is some additional information about the dataset:\n${jsonData.additional_info}.`;
    }
    
    return prompt;
}

async function finalAnnotation(agent, prompt) {
    const conversation = [];
    let currentPrompt = prompt;
    
    while (true) {
        const response = await agent.call(currentPrompt, "user");
        conversation.push(["Final Annotation Agent", response]);
        
        if (response.includes("FINAL ANNOTATION COMPLETED")) {
            break;
        }
        
        currentPrompt = response;
    }
    
    return conversation;
}

async function couplingValidation(agent, annotationResult, onboardingData) {
    let markerList = onboardingData.marker_list || [];
    
    if (Array.isArray(markerList) && 
        markerList.length === 1 && 
        typeof markerList[0] === 'string' && 
        markerList[0].includes(',')) {
        
        const markers = markerList[0].split(/,\s*/).map(m => m.trim()).filter(m => m);
        markerList = markers.join(', ');
    } else if (Array.isArray(markerList)) {
        markerList = markerList.join(', ');
    } else if (typeof markerList === 'string') {
        markerList = markerList;
    } else {
        markerList = String(markerList);
    }
    
    const validationMessage = `Please validate the following annotation result:

Annotation Result:
${annotationResult}

Context:

Marker List: ${markerList}
Additional Info: ${onboardingData.additional_info || 'None'}

Validate the annotation based on this context.`;

    return await agent.call(validationMessage, "final_annotation");
}

async function formatResults(agent, finalAnnotations) {
    const finalText = finalAnnotations.map(msg => msg[1]).join('\n\n');
    return await agent.call(finalText, "user");
}

// ----------------- Agent Class -----------------

class Agent {
    constructor(system = "", model = "google/gemini-2.5-flash-preview", temperature = 0, provider = "openrouter", apiKey = null, reasoningConfig = null) {
        this.system = system;
        this.model = model;
        this.temperature = temperature;
        this.provider = provider;
        this.apiKey = apiKey;
        this.chatHistories = {};
        this.reasoningConfig = reasoningConfig;
    }
    
    async call(message, otherAgentId) {
        if (!this.chatHistories[otherAgentId]) {
            this.chatHistories[otherAgentId] = [];
        }
        
        // Add system message if this is the first message for this conversation
        if (this.chatHistories[otherAgentId].length === 0 && this.system) {
            this.chatHistories[otherAgentId].push({ role: "system", content: this.system });
        }
        
        this.chatHistories[otherAgentId].push({ role: "user", content: message });
        
        try {
            const response = await callLLM(
                message,
                this.provider,
                this.model,
                this.apiKey,
                this.temperature,
                7000, // max_tokens
                this.system,
                { messages: this.chatHistories[otherAgentId] },
                this.reasoningConfig
            );
            
            this.chatHistories[otherAgentId].push({ role: "assistant", content: response });
            
            return response;
        } catch (error) {
            console.error(`Error calling LLM: ${error.message}`);
            throw error;
        }
    }
}

// ----------------- Core Analysis Logic -----------------
// Faithfully ported from Python _run_analysis_logic in main_function_code.py

async function runAnalysisLogic(finalAnnotationAgent, couplingValidatorAgent, formattingAgent, userData, isTissueBlind) {
    const onboardingData = userData;

    // Construct initial prompt
    const initialPrompt = constructPrompt(onboardingData);

    let validationPassed = false;
    let iteration = 0;
    let finalAnnotationConversation = [];
    let allConversations = [];
    let validationResult = "";
    const maxIterations = 3;

    // Validation loop - matches Python logic
    while (!validationPassed && iteration < maxIterations) {
        iteration++;
        console.log(`Starting iteration ${iteration}...`);

        // Build prompt - Python-faithful retry logic with actual feedback included
        let currentPrompt;
        if (iteration === 1) {
            currentPrompt = initialPrompt;
        } else {
            // Include previous response and validation feedback like Python does
            const previousResponse = finalAnnotationConversation.length > 0
                ? finalAnnotationConversation[finalAnnotationConversation.length - 1][1]
                : "";

            currentPrompt = `Previous annotation attempt failed validation. Please review your previous response and the validation feedback, then provide an updated annotation:

Previous response:
${previousResponse}

Validation feedback:
${validationResult}

Original prompt:
${initialPrompt}

Please provide an updated annotation addressing the validation feedback.`;
        }

        finalAnnotationConversation = await finalAnnotation(finalAnnotationAgent, currentPrompt);
        allConversations.push({
            iteration: iteration,
            annotation: [...finalAnnotationConversation]
        });

        // Run validation on the last annotation response
        validationResult = await couplingValidation(
            couplingValidatorAgent,
            finalAnnotationConversation[finalAnnotationConversation.length - 1][1],
            onboardingData
        );
        allConversations[allConversations.length - 1].validation_result = validationResult;

        if (validationResult.includes("VALIDATION PASSED")) {
            validationPassed = true;
            console.log(`Validation passed on iteration ${iteration}`);
        } else {
            console.log(`Validation failed on iteration ${iteration}, retrying...`);
        }

        allConversations[allConversations.length - 1].validation_passed = validationPassed;
    }

    // Select appropriate formatting system based on validation status and tissue context
    // This matches Python logic in _run_analysis_logic
    if (validationPassed) {
        formattingAgent.system = isTissueBlind ? formattingSystemTissueBlind : formattingSystemNonTissueBlind;
    } else {
        formattingAgent.system = formattingSystemFailed;
    }

    // Format results - pass last 2 annotation messages like Python does
    let formattedResult;
    try {
        const lastAnnotations = finalAnnotationConversation.slice(-2);
        const formattingResponse = await formatResults(formattingAgent, lastAnnotations);
        allConversations.push({
            role: "Formatting Agent",
            content: formattingResponse
        });

        formattedResult = extractJsonFromReply(formattingResponse);

        if (!formattedResult) {
            throw new Error("Failed to extract JSON from formatting response");
        }
    } catch (error) {
        console.warn("Formatting failed, using fallback format");
        console.error("Error: Unable to extract JSON from the formatted output.");

        // Fallback formatting
        const annotationText = finalAnnotationConversation
            .map(msg => msg[1])
            .join(' ');

        formattedResult = {
            main_cell_type: "Unknown",
            sub_cell_types: [],
            possible_mixed_cell_types: [],
            feedback: "Failed to parse formatting response",
            next_steps: annotationText.substring(0, 200) + "..."
        };
    }

    // Add metadata like Python does
    formattedResult.iterations = iteration;
    formattedResult.num_markers = Array.isArray(onboardingData.marker_list)
        ? onboardingData.marker_list.length
        : 0;

    // Prepare conversation history
    const conversationHistory = {
        all_iterations: allConversations,
        final_result: formattedResult,
        input_data: onboardingData
    };

    return [formattedResult, conversationHistory];
}

// ----------------- Provider-Specific Analysis Functions -----------------

async function runCellTypeAnalysis(model, temperature, markerList, tissue, species, additionalInfo, validatorInvolvement = "v1", apiKey = null, reasoningConfig = null) {
    // Determine if tissue-blind
    const isTissueBlind = tissue ? ['none', 'tissue blind'].includes(tissue.toLowerCase()) : true;

    // Reasoning config for OpenAI models that support it (GPT-5 series)
    // Use reasoningConfig passed as parameter

    // Create agents
    const finalAnnotationAgent = new Agent(
        isTissueBlind ? finalAnnotationSystemV2 : finalAnnotationSystemV1,
        model,
        temperature,
        "openai",
        apiKey,
        reasoningConfig
    );

    let validatorSystem;
    if (validatorInvolvement === "v0") {
        validatorSystem = couplingValidatorSystemV0;
    } else if (validatorInvolvement === "v1") {
        validatorSystem = isTissueBlind ? couplingValidatorSystemV2 : couplingValidatorSystemV1;
    } else {
        validatorSystem = isTissueBlind ? couplingValidatorSystemV2 : couplingValidatorSystemV1;
    }

    const couplingValidatorAgent = new Agent(validatorSystem, model, temperature, "openai", apiKey, reasoningConfig);
    // Formatting agent system is set dynamically in runAnalysisLogic based on validation status
    const formattingAgent = new Agent("", model, temperature, "openai", apiKey, reasoningConfig);

    // Prepare user data
    const userData = {
        species: species,
        tissue_type: tissue,
        marker_list: markerList
    };
    
    if (additionalInfo && additionalInfo.toLowerCase() !== "no") {
        userData.additional_info = additionalInfo;
    }
    
    return await runAnalysisLogic(finalAnnotationAgent, couplingValidatorAgent, formattingAgent, userData, isTissueBlind);
}

async function runCellTypeAnalysisClaude(model, temperature, markerList, tissue, species, additionalInfo, validatorInvolvement = "v1", apiKey = null, reasoningConfig = null) {
    const isTissueBlind = tissue ? ['none', 'tissue blind'].includes(tissue.toLowerCase()) : true;

    // Reasoning config for Anthropic models that support it (Claude Opus 4.5)
    // Use reasoningConfig passed as parameter

    const finalAnnotationAgent = new Agent(
        isTissueBlind ? finalAnnotationSystemV2 : finalAnnotationSystemV1,
        model,
        temperature,
        "anthropic",
        apiKey,
        reasoningConfig
    );

    let validatorSystem;
    if (validatorInvolvement === "v0") {
        validatorSystem = couplingValidatorSystemV0;
    } else if (validatorInvolvement === "v1") {
        validatorSystem = isTissueBlind ? couplingValidatorSystemV2 : couplingValidatorSystemV1;
    } else {
        validatorSystem = isTissueBlind ? couplingValidatorSystemV2 : couplingValidatorSystemV1;
    }

    const couplingValidatorAgent = new Agent(validatorSystem, model, temperature, "anthropic", apiKey, reasoningConfig);
    // Formatting agent system is set dynamically in runAnalysisLogic based on validation status
    const formattingAgent = new Agent("", model, temperature, "anthropic", apiKey, reasoningConfig);

    const userData = {
        species: species,
        tissue_type: tissue,
        marker_list: markerList
    };
    
    if (additionalInfo && additionalInfo.toLowerCase() !== "no") {
        userData.additional_info = additionalInfo;
    }
    
    return await runAnalysisLogic(finalAnnotationAgent, couplingValidatorAgent, formattingAgent, userData, isTissueBlind);
}

async function runCellTypeAnalysisOpenRouter(model, temperature, markerList, tissue, species, additionalInfo, validatorInvolvement = "v1", apiKey = null, reasoningConfig = null) {
    const isTissueBlind = tissue ? ['none', 'tissue blind'].includes(tissue.toLowerCase()) : true;
    
    // Reasoning config for OpenRouter models that support it
    // Use reasoningConfig passed as parameter

    const finalAnnotationAgent = new Agent(
        isTissueBlind ? finalAnnotationSystemV2 : finalAnnotationSystemV1,
        model,
        temperature,
        "openrouter",
        apiKey,
        reasoningConfig
    );

    let validatorSystem;
    if (validatorInvolvement === "v0") {
        validatorSystem = couplingValidatorSystemV0;
    } else if (validatorInvolvement === "v1") {
        validatorSystem = isTissueBlind ? couplingValidatorSystemV2 : couplingValidatorSystemV1;
    } else {
        validatorSystem = isTissueBlind ? couplingValidatorSystemV2 : couplingValidatorSystemV1;
    }

    const couplingValidatorAgent = new Agent(validatorSystem, model, temperature, "openrouter", apiKey, reasoningConfig);
    // Formatting agent system is set dynamically in runAnalysisLogic based on validation status
    const formattingAgent = new Agent("", model, temperature, "openrouter", apiKey, reasoningConfig);

    const userData = {
        species: species,
        tissue_type: tissue,
        marker_list: markerList
    };
    
    if (additionalInfo && additionalInfo.toLowerCase() !== "no") {
        userData.additional_info = additionalInfo;
    }
    
    return await runAnalysisLogic(finalAnnotationAgent, couplingValidatorAgent, formattingAgent, userData, isTissueBlind);
}

async function runCellTypeAnalysisCustom(baseUrl, apiKey, model, temperature, markerList, tissue, species, additionalInfo, validatorInvolvement = "v1", reasoningConfig = null) {
    const isTissueBlind = tissue ? ['none', 'tissue blind'].includes(tissue.toLowerCase()) : true;

    // Reasoning config for custom endpoints that support it
    // Use reasoningConfig passed as parameter

    const finalAnnotationAgent = new Agent(
        isTissueBlind ? finalAnnotationSystemV2 : finalAnnotationSystemV1,
        model,
        temperature,
        baseUrl,
        apiKey,
        reasoningConfig
    );

    let validatorSystem;
    if (validatorInvolvement === "v0") {
        validatorSystem = couplingValidatorSystemV0;
    } else if (validatorInvolvement === "v1") {
        validatorSystem = isTissueBlind ? couplingValidatorSystemV2 : couplingValidatorSystemV1;
    } else {
        validatorSystem = isTissueBlind ? couplingValidatorSystemV2 : couplingValidatorSystemV1;
    }

    const couplingValidatorAgent = new Agent(validatorSystem, model, temperature, baseUrl, apiKey, reasoningConfig);
    // Formatting agent system is set dynamically in runAnalysisLogic based on validation status
    const formattingAgent = new Agent("", model, temperature, baseUrl, apiKey, reasoningConfig);

    const userData = {
        species: species,
        tissue_type: tissue,
        marker_list: markerList
    };
    
    if (additionalInfo && additionalInfo.toLowerCase() !== "no") {
        userData.additional_info = additionalInfo;
    }
    
    return await runAnalysisLogic(finalAnnotationAgent, couplingValidatorAgent, formattingAgent, userData, isTissueBlind);
}

// ----------------- Main runCASSIA Function -----------------

/**
 * Wrapper function to run cell type analysis using OpenAI, Anthropic, OpenRouter, or a custom OpenAI-compatible provider.
 * 
 * @param {string} model - Model name to use
 * @param {number} temperature - Temperature parameter for the model
 * @param {Array} markerList - List of markers to analyze
 * @param {string} tissue - Tissue type
 * @param {string} species - Species type
 * @param {string} additionalInfo - Additional information for the analysis
 * @param {string} provider - AI provider to use ('openai', 'anthropic', 'openrouter', or a base URL)
 * @param {string} validatorInvolvement - Validator involvement level ('v0', 'v1', 'v2')
 * @returns {Promise<Array>} [analysis_result, conversation_history]
 */
export async function runCASSIA(
    model = "google/gemini-2.5-flash-preview",
    temperature = 0,
    markerList = null,
    tissue = "lung",
    species = "human",
    additionalInfo = null,
    provider = "openrouter",
    validatorInvolvement = "v1",
    apiKey = null,
    reasoningEffort = null
) {
    // Build reasoning config from effort level (case-insensitive)
    const reasoningConfig = reasoningEffort && reasoningEffort.toLowerCase() !== 'none'
        ? { effort: reasoningEffort.toLowerCase() }
        : null;

    if (provider.toLowerCase() === "openai") {
        return await runCellTypeAnalysis(model, temperature, markerList, tissue, species, additionalInfo, validatorInvolvement, apiKey, reasoningConfig);
    } else if (provider.toLowerCase() === "anthropic") {
        return await runCellTypeAnalysisClaude(model, temperature, markerList, tissue, species, additionalInfo, validatorInvolvement, apiKey, reasoningConfig);
    } else if (provider.toLowerCase() === "openrouter") {
        return await runCellTypeAnalysisOpenRouter(model, temperature, markerList, tissue, species, additionalInfo, validatorInvolvement, apiKey, reasoningConfig);
    } else if (provider.toLowerCase().startsWith("http")) {
        if (!apiKey) {
            throw new Error("API key is required for custom endpoint");
        }
        return await runCellTypeAnalysisCustom(
            provider,
            apiKey,
            model,
            temperature,
            markerList,
            tissue,
            species,
            additionalInfo,
            validatorInvolvement,
            reasoningConfig
        );
    } else {
        throw new Error("Provider must be either 'openai', 'anthropic', 'openrouter', or a base URL (http...)");
    }
}

export default { runCASSIA };