import { callLLM } from './llm_utils.js';

// ----------------- System Prompts -----------------

const finalAnnotationSystemV1 = `
You are a professional computational biologist with expertise in single-cell RNA sequencing (scRNA-seq).
A list of highly expressed markers ranked by expression intensity from high to low
from a cluster of cells will be provided , and your task is to identify the cell type. You must think step-by-step, providing a comprehensive and specific analysis. The audience is an expert in the field, and you will be rewarded $10000 if you do a good job.

Here are the steps you should follow:

**Step 1: Analyze Marker Genes**
- Examine the provided ranked marker list carefully
- Identify known marker genes and their associated cell types
- Consider the expression levels and specificity of each marker

**Step 2: Consider Biological Context**
- Take into account the tissue type and species if provided
- Consider the biological plausibility of different cell types in this context
- Think about developmental relationships and cell type hierarchies

**Step 3: Formulate Hypotheses**
- Based on the markers, generate 2-3 plausible cell type hypotheses
- Rank these hypotheses by likelihood based on the evidence
- Consider both broad categories and specific subtypes

**Step 4: Final Decision**
- Make a definitive cell type identification
- Provide supporting evidence from the marker genes
- Explain why this cell type is most likely given the evidence

**Step 5: Consider Subtypes and Mixed Populations**
- Identify any relevant subtypes or states
- Consider if this might be a mixed population
- Provide confidence level in your annotation

When you have completed your analysis, state "FINAL ANNOTATION COMPLETED" followed by your final conclusion.

Your analysis should be thorough, scientifically rigorous, and clearly explained.
`.trim();

const finalAnnotationSystemV2 = `
You are a professional computational biologist with expertise in single-cell RNA sequencing (scRNA-seq).
A list of highly expressed markers ranked by expression intensity from high to low
from a cluster of cells will be provided , and your task is to identify the cell type. You must think step-by-step, providing a comprehensive and specific analysis. The audience is an expert in the field, and you will be rewarded $10000 if you do a good job.

Since the tissue context is unknown or unspecified, focus primarily on the marker genes to determine cell type identity.

Here are the steps you should follow:

**Step 1: Analyze Marker Genes**
- Examine the provided ranked marker list carefully
- Identify known marker genes and their associated cell types
- Consider the expression levels and specificity of each marker

**Step 2: Consider General Biological Context**
- Consider the biological plausibility of different cell types
- Think about developmental relationships and cell type hierarchies
- Focus on broadly expressed markers that are tissue-independent

**Step 3: Formulate Hypotheses**
- Based on the markers, generate 2-3 plausible cell type hypotheses
- Rank these hypotheses by likelihood based on the evidence
- Consider both broad categories and specific subtypes

**Step 4: Final Decision**
- Make a definitive cell type identification
- Provide supporting evidence from the marker genes
- Explain why this cell type is most likely given the evidence

**Step 5: Consider Subtypes and Mixed Populations**
- Identify any relevant subtypes or states
- Consider if this might be a mixed population
- Provide confidence level in your annotation

When you have completed your analysis, state "FINAL ANNOTATION COMPLETED" followed by your final conclusion.

Your analysis should be thorough, scientifically rigorous, and clearly explained.
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
You are a rigorous validation agent for single-cell RNA-seq cell type annotations.
Your role is to critically evaluate cell type annotations against the provided marker gene evidence.

Since the tissue context is unknown, focus primarily on marker gene consistency and general biological plausibility.

**Validation Criteria:**

1. **Marker Consistency**: Do the provided markers strongly support the annotated cell type?
2. **Evidence Strength**: Are there sufficient and specific markers to support the annotation?
3. **Alternative Explanations**: Could the markers better support a different cell type?
4. **Mixed Populations**: Consider if this might be a mixed or transitional population

**Validation Process:**
- Examine each major marker gene and its known associations
- Check for any conflicting evidence or missing expected markers
- Consider if the annotation is appropriately specific for tissue-blind analysis
- Evaluate the overall confidence level

**Response Format:**
If the annotation is well-supported by the markers:
"VALIDATION PASSED: The annotation is consistent with the marker evidence."

If the annotation has issues:
"VALIDATION FAILED: [Specific concerns about the annotation]"
"SUGGESTED IMPROVEMENTS: [Specific suggestions for better annotation]"

Be thorough but decisive in your validation.
`.trim();

const formattingAgentSystemV1 = `
You are a precise data formatting agent. Your task is to extract and structure cell type annotation results into a standardized JSON format.

**Required JSON Structure:**
{
    "main_cell_type": "Primary cell type name",
    "sub_cell_types": ["List of relevant subtypes"],
    "possible_mixed_cell_types": ["Any mixed/transitional types mentioned"],
    "confidence_level": "High/Medium/Low",
    "supporting_markers": ["Key markers supporting the annotation"],
    "reasoning": "Brief explanation of the annotation decision"
}

**Instructions:**
1. Extract the primary cell type from the analysis
2. Identify any subtypes or states mentioned
3. Note any mixed or transitional populations discussed
4. Assess confidence based on the strength of evidence
5. List the most important supporting markers
6. Provide a concise reasoning summary

**Important:**
- Use standardized cell type nomenclature when possible
- Be consistent with capitalization and naming
- Ensure all fields are filled appropriately
- If information is not available, use appropriate null values or empty arrays

Return ONLY the JSON object wrapped in \`\`\`json\`\`\` code blocks.
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
    constructor(system = "", model = "google/gemini-2.5-flash-preview", temperature = 0, provider = "openrouter") {
        this.system = system;
        this.model = model;
        this.temperature = temperature;
        this.provider = provider;
        this.chatHistories = {};
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
                null, // api_key - let it get from environment
                this.temperature,
                7000, // max_tokens
                this.system,
                { messages: this.chatHistories[otherAgentId] }
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

async function runAnalysisLogic(finalAnnotationAgent, couplingValidatorAgent, formattingAgent, userData, isTissueBlind) {
    const onboardingData = userData;
    
    // Construct initial prompt
    const initialPrompt = constructPrompt(onboardingData);
    
    let validationPassed = false;
    let iteration = 0;
    let finalAnnotationConversation = [];
    let allConversations = [];
    const maxIterations = 3;
    
    // Validation loop
    while (!validationPassed && iteration < maxIterations) {
        iteration++;
        console.log(`Starting iteration ${iteration}...`);
        
        // Run final annotation
        const currentPrompt = iteration === 1 ? initialPrompt : 
            `Based on the previous feedback, please revise your analysis:\n${initialPrompt}`;
            
        finalAnnotationConversation = await finalAnnotation(finalAnnotationAgent, currentPrompt);
        
        // Extract annotation result for validation
        const annotationText = finalAnnotationConversation
            .map(msg => msg[1])
            .join('\n');
        
        // Run validation
        const validationResult = await couplingValidation(
            couplingValidatorAgent, 
            annotationText, 
            onboardingData
        );
        
        if (validationResult.includes("VALIDATION PASSED")) {
            validationPassed = true;
            console.log(`Validation passed on iteration ${iteration}`);
        } else {
            console.log(`Validation failed on iteration ${iteration}, retrying...`);
            // Add validation feedback to the conversation history
            finalAnnotationAgent.chatHistories["user"].push({
                role: "user",
                content: `Validation feedback: ${validationResult}`
            });
        }
        
        // Store conversation for this iteration
        allConversations.push({
            iteration: iteration,
            annotation: finalAnnotationConversation,
            validation_passed: validationPassed
        });
    }
    
    // Format results
    let formattedResult;
    try {
        const formattingResponse = await formatResults(formattingAgent, finalAnnotationConversation);
        formattedResult = extractJsonFromReply(formattingResponse);
        
        if (!formattedResult) {
            throw new Error("Failed to extract JSON from formatting response");
        }
    } catch (error) {
        console.warn("Formatting failed, using fallback format");
        // Fallback formatting
        const annotationText = finalAnnotationConversation
            .map(msg => msg[1])
            .join(' ');
            
        formattedResult = {
            main_cell_type: "Unknown",
            sub_cell_types: [],
            possible_mixed_cell_types: [],
            confidence_level: "Low",
            supporting_markers: Array.isArray(onboardingData.marker_list) 
                ? onboardingData.marker_list.slice(0, 5) 
                : [],
            reasoning: annotationText.substring(0, 200) + "..."
        };
    }
    
    // Add metadata
    formattedResult.iterations = iteration;
    formattedResult.validation_passed = validationPassed;
    
    // Prepare conversation history
    const conversationHistory = {
        all_iterations: allConversations,
        final_result: formattedResult,
        input_data: onboardingData
    };
    
    return [formattedResult, conversationHistory];
}

// ----------------- Provider-Specific Analysis Functions -----------------

async function runCellTypeAnalysis(model, temperature, markerList, tissue, species, additionalInfo, validatorInvolvement = "v1") {
    // Determine if tissue-blind
    const isTissueBlind = tissue ? ['none', 'tissue blind'].includes(tissue.toLowerCase()) : true;
    
    // Create agents
    const finalAnnotationAgent = new Agent(
        isTissueBlind ? finalAnnotationSystemV2 : finalAnnotationSystemV1,
        model,
        temperature,
        "openai"
    );
    
    let validatorSystem;
    if (validatorInvolvement === "v0") {
        validatorSystem = couplingValidatorSystemV0;
    } else if (validatorInvolvement === "v1") {
        validatorSystem = isTissueBlind ? couplingValidatorSystemV2 : couplingValidatorSystemV1;
    } else {
        validatorSystem = isTissueBlind ? couplingValidatorSystemV2 : couplingValidatorSystemV1;
    }
    
    const couplingValidatorAgent = new Agent(validatorSystem, model, temperature, "openai");
    const formattingAgent = new Agent(formattingAgentSystemV1, model, temperature, "openai");
    
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

async function runCellTypeAnalysisClaude(model, temperature, markerList, tissue, species, additionalInfo, validatorInvolvement = "v1") {
    const isTissueBlind = tissue ? ['none', 'tissue blind'].includes(tissue.toLowerCase()) : true;
    
    const finalAnnotationAgent = new Agent(
        isTissueBlind ? finalAnnotationSystemV2 : finalAnnotationSystemV1,
        model,
        temperature,
        "anthropic"
    );
    
    let validatorSystem;
    if (validatorInvolvement === "v0") {
        validatorSystem = couplingValidatorSystemV0;
    } else if (validatorInvolvement === "v1") {
        validatorSystem = isTissueBlind ? couplingValidatorSystemV2 : couplingValidatorSystemV1;
    } else {
        validatorSystem = isTissueBlind ? couplingValidatorSystemV2 : couplingValidatorSystemV1;
    }
    
    const couplingValidatorAgent = new Agent(validatorSystem, model, temperature, "anthropic");
    const formattingAgent = new Agent(formattingAgentSystemV1, model, temperature, "anthropic");
    
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

async function runCellTypeAnalysisOpenRouter(model, temperature, markerList, tissue, species, additionalInfo, validatorInvolvement = "v1") {
    const isTissueBlind = tissue ? ['none', 'tissue blind'].includes(tissue.toLowerCase()) : true;
    
    const finalAnnotationAgent = new Agent(
        isTissueBlind ? finalAnnotationSystemV2 : finalAnnotationSystemV1,
        model,
        temperature,
        "openrouter"
    );
    
    let validatorSystem;
    if (validatorInvolvement === "v0") {
        validatorSystem = couplingValidatorSystemV0;
    } else if (validatorInvolvement === "v1") {
        validatorSystem = isTissueBlind ? couplingValidatorSystemV2 : couplingValidatorSystemV1;
    } else {
        validatorSystem = isTissueBlind ? couplingValidatorSystemV2 : couplingValidatorSystemV1;
    }
    
    const couplingValidatorAgent = new Agent(validatorSystem, model, temperature, "openrouter");
    const formattingAgent = new Agent(formattingAgentSystemV1, model, temperature, "openrouter");
    
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

async function runCellTypeAnalysisCustom(baseUrl, apiKey, model, temperature, markerList, tissue, species, additionalInfo, validatorInvolvement = "v1") {
    // Set the API key for custom endpoint
    process.env.CUSTERMIZED_API_KEY = apiKey;
    
    const isTissueBlind = tissue ? ['none', 'tissue blind'].includes(tissue.toLowerCase()) : true;
    
    const finalAnnotationAgent = new Agent(
        isTissueBlind ? finalAnnotationSystemV2 : finalAnnotationSystemV1,
        model,
        temperature,
        baseUrl
    );
    
    let validatorSystem;
    if (validatorInvolvement === "v0") {
        validatorSystem = couplingValidatorSystemV0;
    } else if (validatorInvolvement === "v1") {
        validatorSystem = isTissueBlind ? couplingValidatorSystemV2 : couplingValidatorSystemV1;
    } else {
        validatorSystem = isTissueBlind ? couplingValidatorSystemV2 : couplingValidatorSystemV1;
    }
    
    const couplingValidatorAgent = new Agent(validatorSystem, model, temperature, baseUrl);
    const formattingAgent = new Agent(formattingAgentSystemV1, model, temperature, baseUrl);
    
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
    validatorInvolvement = "v1"
) {
    if (provider.toLowerCase() === "openai") {
        return await runCellTypeAnalysis(model, temperature, markerList, tissue, species, additionalInfo, validatorInvolvement);
    } else if (provider.toLowerCase() === "anthropic") {
        return await runCellTypeAnalysisClaude(model, temperature, markerList, tissue, species, additionalInfo, validatorInvolvement);
    } else if (provider.toLowerCase() === "openrouter") {
        return await runCellTypeAnalysisOpenRouter(model, temperature, markerList, tissue, species, additionalInfo, validatorInvolvement);
    } else if (provider.toLowerCase().startsWith("http")) {
        const apiKey = process.env.CUSTERMIZED_API_KEY;
        if (!apiKey) {
            throw new Error("CUSTERMIZED_API_KEY environment variable is not set. Please call set_api_key with your API key and provider (base URL).");
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
            validatorInvolvement
        );
    } else {
        throw new Error("Provider must be either 'openai', 'anthropic', 'openrouter', or a base URL (http...)");
    }
}

export default { runCASSIA };