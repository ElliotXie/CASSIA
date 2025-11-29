import { callLLM } from './llm_utils.js';

// ----------------- Helper Functions and Prompts -----------------

/**
 * Extracts a JSON object from a string, typically the LLM's reply.
 * @param {string} reply - The LLM response string
 * @returns {object|null} Parsed JSON object or null if not found
 */
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

/**
 * Constructs the initial prompt for the cell type annotation task.
 * @param {object} jsonData - Input data containing marker list and metadata
 * @returns {string} Constructed prompt
 */
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

/**
 * Manages the conversation with the final annotation agent.
 * @param {Function} agent - The agent function
 * @param {string} prompt - Initial prompt
 * @returns {Array} Conversation history
 */
async function finalAnnotation(agent, prompt) {
    const conversation = [];
    let currentPrompt = prompt;
    
    while (true) {
        const response = await agent(currentPrompt, "user");
        conversation.push(["Final Annotation Agent", response]);
        
        if (response.includes("FINAL ANNOTATION COMPLETED")) {
            break;
        }
        
        currentPrompt = response;
    }
    
    return conversation;
}

/**
 * Constructs and sends the validation prompt to the coupling validator agent.
 * @param {Function} agent - The validator agent function
 * @param {string} annotationResult - Result from annotation
 * @param {object} onboardingData - Original input data
 * @returns {Promise<string>} Validation response
 */
async function couplingValidation(agent, annotationResult, onboardingData) {
    // Handle different marker_list formats
    let markerList = onboardingData.marker_list || [];
    
    // If marker_list is an array with one element that contains commas
    if (Array.isArray(markerList) && 
        markerList.length === 1 && 
        typeof markerList[0] === 'string' && 
        markerList[0].includes(',')) {
        
        // Split the single string into individual genes
        const markers = markerList[0].split(/,\s*/).map(m => m.trim()).filter(m => m);
        markerList = markers.join(', ');
    } else if (Array.isArray(markerList)) {
        // Normal case: array of individual genes
        markerList = markerList.join(', ');
    } else if (typeof markerList === 'string') {
        // Single string case: use as-is
        markerList = markerList;
    } else {
        // Fallback
        markerList = String(markerList);
    }
    
    const validationMessage = `Please validate the following annotation result:

Annotation Result:
${annotationResult}

Context:

Marker List: ${markerList}
Additional Info: ${onboardingData.additional_info || 'None'}

Validate the annotation based on this context.`;

    return await agent(validationMessage, "final_annotation");
}

/**
 * Calls the formatting agent and returns its raw response.
 * @param {Function} agent - The formatting agent function
 * @param {Array} finalAnnotations - Final annotation conversation
 * @returns {Promise<string>} Formatted response
 */
async function formatResults(agent, finalAnnotations) {
    const finalText = finalAnnotations.map(msg => msg[1]).join('\n\n');
    return await agent(finalText, "user");
}

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
`;

const couplingValidatorSystemV1 = `
You are a rigorous validation agent for single-cell RNA-seq cell type annotations.
Your role is to critically evaluate cell type annotations against the provided marker gene evidence.

**Validation Criteria:**

1. **Marker Consistency**: Do the provided markers strongly support the annotated cell type?
2. **Biological Plausibility**: Is this cell type reasonable given the tissue/species context?
3. **Evidence Strength**: Are there sufficient and specific markers to support the annotation?
4. **Alternative Explanations**: Could the markers better support a different cell type?

**Validation Process:**
- Examine each major marker gene and its known associations
- Check for any conflicting evidence or missing expected markers
- Consider if the annotation is too broad or too specific
- Evaluate the overall confidence level

**Response Format:**
If the annotation is well-supported by the markers:
"VALIDATION PASSED: The annotation is consistent with the marker evidence."

If the annotation has issues:
"VALIDATION FAILED: [Specific concerns about the annotation]"
"SUGGESTED IMPROVEMENTS: [Specific suggestions for better annotation]"

Be thorough but decisive in your validation.
`;

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
`;

// ----------------- Agent Classes -----------------

/**
 * Agent class for managing LLM conversations with system prompts
 */
class Agent {
    constructor(system = "", model = "google/gemini-2.5-flash-preview", temperature = 0, provider = "openrouter") {
        this.system = system;
        this.model = model;
        this.temperature = temperature;
        this.provider = provider;
        this.chatHistories = {};
    }
    
    async call(message, otherAgentId) {
        // Initialize conversation history for this agent pair if not exists
        if (!this.chatHistories[otherAgentId]) {
            this.chatHistories[otherAgentId] = [];
        }
        
        // Add user message to history
        this.chatHistories[otherAgentId].push({ role: "user", content: message });
        
        try {
            const response = await callLLM({
                prompt: message,
                systemPrompt: this.system,
                model: this.model,
                temperature: this.temperature,
                provider: this.provider,
                additionalParams: {
                    messages: this.chatHistories[otherAgentId]
                }
            });
            
            // Add assistant response to history
            this.chatHistories[otherAgentId].push({ role: "assistant", content: response });
            
            return response;
        } catch (error) {
            console.error(`Error calling LLM: ${error.message}`);
            throw error;
        }
    }
}

// ----------------- Core Analysis Logic -----------------

/**
 * Core analysis logic for single cell type annotation
 * @param {object} params - Analysis parameters
 * @returns {Promise<object>} Analysis results
 */
async function runAnalysisLogic({
    model = "google/gemini-2.5-flash-preview",
    temperature = 0,
    markerList = null,
    tissue = "lung",
    species = "human",
    additionalInfo = null,
    provider = "openrouter",
    validatorInvolvement = "v1",
    maxIterations = 3
}) {
    // Prepare input data
    const onboardingData = {
        marker_list: markerList,
        tissue_type: tissue,
        species: species,
        additional_info: additionalInfo || "No additional information provided."
    };
    
    // Initialize agents
    const finalAnnotationAgent = new Agent(finalAnnotationSystemV1, model, temperature, provider);
    const couplingValidatorAgent = new Agent(couplingValidatorSystemV1, model, temperature, provider);
    const formattingAgent = new Agent(formattingAgentSystemV1, model, temperature, provider);
    
    // Construct initial prompt
    const initialPrompt = constructPrompt(onboardingData);
    
    let validationPassed = false;
    let iteration = 0;
    let finalAnnotationConversation = [];
    let allConversations = [];
    
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
        
        // Run validation if validator involvement is enabled
        if (validatorInvolvement !== "v0") {
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
        } else {
            validationPassed = true; // Skip validation
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
            supporting_markers: markerList ? markerList.slice(0, 5) : [],
            reasoning: annotationText.substring(0, 200) + "..."
        };
    }
    
    // Add metadata
    formattedResult.iterations = iteration;
    formattedResult.validation_passed = validationPassed;
    formattedResult.provider = provider;
    formattedResult.model = model;
    
    // Prepare conversation history
    const conversationHistory = {
        all_iterations: allConversations,
        final_result: formattedResult,
        input_data: onboardingData
    };
    
    return { result: formattedResult, conversationHistory };
}

// ----------------- Provider-Specific Functions -----------------

/**
 * Run cell type analysis using OpenAI
 */
async function runCellTypeAnalysis(params) {
    return await runAnalysisLogic({ ...params, provider: "openai" });
}

/**
 * Run cell type analysis using Anthropic Claude
 */
async function runCellTypeAnalysisClaude(params) {
    return await runAnalysisLogic({ ...params, provider: "anthropic" });
}

/**
 * Run cell type analysis using OpenRouter
 */
async function runCellTypeAnalysisOpenRouter(params) {
    return await runAnalysisLogic({ ...params, provider: "openrouter" });
}

/**
 * Run cell type analysis using custom endpoint
 */
async function runCellTypeAnalysisCustom(params) {
    return await runAnalysisLogic({ ...params, provider: params.customEndpoint });
}

export {
    runCellTypeAnalysis,
    runCellTypeAnalysisClaude, 
    runCellTypeAnalysisOpenRouter,
    runCellTypeAnalysisCustom,
    runAnalysisLogic,
    Agent,
    extractJsonFromReply,
    constructPrompt,
    finalAnnotation,
    couplingValidation,
    formatResults
};