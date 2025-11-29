/**
 * CASSIA Symphony Compare - Multi-Model Cell Type Comparison with AI Consensus Building
 * 100% Python-compatible JavaScript implementation for the browser
 */

import { callLLM } from './llm_utils.js';
import { MODEL_PRESETS } from '../config/model-presets';

/**
 * Extract scores and reasoning for each celltype from the LLM response.
 */
export function extractCelltypeScores(responseText, celltypes) {
    console.log(`🔍 Extracting scores for: ${celltypes.join(', ')}`);
    console.log(`📄 Response length: ${responseText.length} characters`);
    console.log(`📝 Raw response preview: ${responseText.substring(0, 500)}...`);
    const results = {};

    for (const celltype of celltypes) {
        results[celltype] = { score: "No score found", reasoning: "No reasoning found" };
    }

    try {
        // First, try to parse the entire response as JSON, as some models might return a structured object.
        const cleanedResponse = responseText.replace(/```json|```/g, '').trim();
        const parsedJson = JSON.parse(cleanedResponse);

        if (typeof parsedJson === 'object' && parsedJson !== null) {
            for (const celltype of celltypes) {
                if (parsedJson[celltype] && typeof parsedJson[celltype] === 'object') {
                    results[celltype] = {
                        score: parsedJson[celltype].score !== undefined ? String(parsedJson[celltype].score) : "No score found",
                        reasoning: parsedJson[celltype].reasoning || "No reasoning found",
                    };
                }
            }
            // If we successfully parsed and extracted from JSON, we can return.
            const hasValidScores = Object.values(results).some(result => result.score !== "No score found");
            if (hasValidScores) {
                console.log(`✅ JSON parsing successful`);
                return results;
            }
        }
    } catch (e) {
        // Not a valid JSON, proceed with regex matching.
        console.log(`❌ JSON parsing failed, using XML/regex patterns`);
    }

    // Parse each celltype block individually
    for (const celltype of celltypes) {
        const escapedCelltype = celltype.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
        
        // Pattern 1: Standard format <celltype>NAME</celltype> <reasoning>...</reasoning> <score>N</score>
        const blockPattern = new RegExp(
            `<celltype>${escapedCelltype}</celltype>\\s*<reasoning>([\\s\\S]*?)</reasoning>\\s*<score>([\\d.]+)</score>`, 
            'i'
        );
        
        const blockMatch = responseText.match(blockPattern);
        
        if (blockMatch) {
            results[celltype].reasoning = blockMatch[1].trim();
            results[celltype].score = blockMatch[2];
            console.log(`✅ Successfully parsed ${celltype}: score=${blockMatch[2]}, reasoning length=${blockMatch[1].trim().length}`);
            continue;
        }
        
        console.log(`❌ Primary pattern failed for ${celltype}, trying alternative XML patterns...`);
        
        // Pattern 2: More flexible spacing and ordering
        const altPattern = new RegExp(
            `<celltype>\\s*${escapedCelltype}\\s*</celltype>[\\s\\S]*?<reasoning>\\s*([\\s\\S]*?)\\s*</reasoning>[\\s\\S]*?<score>\\s*([\\d.]+)\\s*</score>`,
            'i'
        );
        
        const altMatch = responseText.match(altPattern);
        if (altMatch) {
            results[celltype].reasoning = altMatch[1].trim();
            results[celltype].score = altMatch[2];
            console.log(`✅ Successfully parsed ${celltype} with alternative pattern: score=${altMatch[2]}`);
            continue;
        }
        
        console.log(`❌ All XML parsing failed for ${celltype}, trying fallback patterns...`);
        
        // Fallback patterns for non-XML responses
        const patterns = [
            new RegExp(`${escapedCelltype}[:\\s]*([0-9]{1,3})(?:\\/100|\\s*(?:out of|\\/)\\s*100)?`, 'i'),
            new RegExp(`${escapedCelltype}[\\s\\-]*(?:score|rating)[:\\s]*([0-9]{1,3})`, 'i'),
            new RegExp(`${escapedCelltype}[\\s\\S]{0,50}?\\b([0-9]{1,3})\\b`, 'i'),
            new RegExp(`(?:score|rating)[:\\s]*([0-9]{1,3})[\\s\\S]{0,30}?${escapedCelltype}`, 'i')
        ];
        
        for (const pattern of patterns) {
            const match = responseText.match(pattern);
            if (match && match[1]) {
                const score = parseInt(match[1]);
                if (score >= 0 && score <= 100) {
                    results[celltype].score = match[1];
                    console.log(`✅ Found fallback score for ${celltype}: ${match[1]}`);
                    break;
                }
            }
        }
        
        // Try to extract reasoning even without XML tags
        const reasoningPatterns = [
            new RegExp(`${escapedCelltype}[:\\s\\-]([^<>]{50,200})`, 'i'),
            new RegExp(`(?:reasoning|analysis|explanation)[\\s\\S]{0,20}${escapedCelltype}[:\\s\\-]([^<>]{50,200})`, 'i')
        ];
        
        for (const pattern of reasoningPatterns) {
            const match = responseText.match(pattern);
            if (match && match[1] && match[1].trim().length > 20) {
                results[celltype].reasoning = match[1].trim();
                console.log(`✅ Found fallback reasoning for ${celltype}`);
                break;
            }
        }
    }
    
    console.log(`📊 Final extraction results:`, results);
    return results;
}


/**
 * Extracts the content of the <discussion> tag from the response.
 */
export function extractDiscussion(responseText) {
    const discussionPattern = /<discussion>([\s\S]*?)<\/discussion>/i;
    const match = responseText.match(discussionPattern);
    return match ? match[1].trim() : "No discussion found";
}

/**
 * Helper function to make a single API call to a model.
 */
async function callModel(model, prompt, apiKey, provider) {
    console.log(`🔄 Calling model: ${model}`);
    try {
        const modelResponse = await callLLM(prompt, provider, model, apiKey, 0.5, 4000);
        console.log(`✅ Model ${model} responded with ${modelResponse.length} characters`);
        console.log(`📝 Response preview: ${modelResponse.substring(0, 150)}...`);
        return { response: modelResponse, status: 'success' };
    } catch (error) {
        console.error(`❌ Error calling model ${model}:`, error);
        return { response: `Error: ${error.message}`, status: 'error' };
    }
}


/**
 * Generate an HTML report for cell type comparison results.
 */
export function generateComparisonHtmlReport(allResults) {
    if (!allResults || allResults.length === 0) {
        return "<html><body><h1>No results to display</h1></body></html>";
    }

    const celltypes = [...new Set(allResults.flatMap(r => Object.keys(r.extracted_scores || {})))].sort();
    const rounds = [...new Set(allResults.map(r => r.round))].sort((a, b) => {
        if (a === 'initial') return -1;
        if (b === 'initial') return 1;
        const aNum = parseInt(a.split('_')[1] || 0);
        const bNum = parseInt(b.split('_')[1] || 0);
        return aNum - bNum;
    });

    let html = `
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>Symphony Compare Report</title>
        <style>
            body { font-family: sans-serif; margin: 2em; }
            h1, h2, h3 { color: #333; }
            table { border-collapse: collapse; width: 100%; margin-bottom: 2em; }
            th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
            th { background-color: #f2f2f2; }
            .discussion { border-left: 3px solid #007bff; padding-left: 1em; margin: 1em 0; background-color: #f8f9fa; }
        </style>
    </head>
    <body>
        <h1>Symphony Compare Analysis Report</h1>
    `;

    for (const round of rounds) {
        html += `<h2>${round.replace('_', ' ').toUpperCase()}</h2>`;
        const roundResults = allResults.filter(r => r.round === round);

        html += '<table><tr><th>Model</th>';
        for (const ct of celltypes) {
            html += `<th>${ct} Score</th>`;
        }
        html += '</tr>';

        for (const result of roundResults) {
            html += `<tr><td>${result.model}</td>`;
            for (const ct of celltypes) {
                const score = result.extracted_scores[ct]?.score || 'N/A';
                html += `<td>${score}</td>`;
            }
            html += '</tr>';
        }
        html += '</table>';

        if (round.startsWith('discussion')) {
            const discussionText = roundResults[0]?.discussion || "No discussion recorded for this round.";
            html += `<div class="discussion"><h3>Discussion Summary</h3><p>${discussionText}</p></div>`;
        }

        for (const result of roundResults) {
            html += `<h3>${result.model} Reasoning (${result.round})</h3>`;
            for (const ct of celltypes) {
                 html += `<h4>${ct}</h4><p>${result.extracted_scores[ct]?.reasoning || 'No reasoning provided.'}</p>`;
            }
        }
    }

    html += '</body></html>';
    return html;
}

/**
 * Main function to run the Symphony Compare analysis.
 * Updated to match Python version parameters and logic.
 */
export async function symphonyCompare({
    tissue,
    celltypes,
    markerSet,
    species = "human",
    modelPreset = "performance",
    customModels = null,
    enableDiscussion = true,
    maxDiscussionRounds = 2,
    consensusThreshold = 0.8,
    apiKey,
    provider = "openrouter"
}) {

    // Input validation
    if (!celltypes || celltypes.length < 2 || celltypes.length > 4) {
        throw new Error("Please provide 2-4 cell types to compare");
    }

    // Use centralized model presets
    const modelPresets = {};
    for (const [key, preset] of Object.entries(MODEL_PRESETS)) {
        modelPresets[key] = preset.models;
    }
    
    console.log(`🤖 Using model preset: ${modelPreset}`);

    // Researcher persona names for working models
    const modelPersonas = {
        "google/gemini-2.5-flash": "Dr. Ada Lovelace",
        "anthropic/claude-3-haiku": "Dr. Alan Turing", 
        "openai/gpt-4o-mini": "Dr. Marie Curie",
        "anthropic/claude-3.5-sonnet": "Dr. Claude Shannon",
        "openai/gpt-4o": "Dr. Albert Einstein",
        "meta-llama/llama-3.1-405b-instruct": "Dr. Rosalind Franklin"
    };

    // Select models based on preset or custom list
    let modelList;
    if (modelPreset === "custom" && customModels) {
        modelList = customModels;
        console.log(`🔧 Using custom models: ${customModels.join(', ')}`);
    } else if (modelPresets[modelPreset]) {
        modelList = modelPresets[modelPreset];
        console.log(`🎵 Using ${modelPreset} preset: ${modelList.join(', ')}`);
    } else {
        console.warn(`Unknown preset '${modelPreset}'. Using 'performance' preset.`);
        modelList = modelPresets["performance"];
        console.log(`⚠️ Fallback to performance preset: ${modelList.join(', ')}`);
    }

    // Get persona names
    const modelToPersona = {};
    for (const m of modelList) {
        modelToPersona[m] = modelPersonas[m] || `Researcher_${m.split('/').pop()}`;
    }
    console.log(`👥 Model personas:`, modelToPersona);

    console.log(`🎼 CASSIA Symphony Compare - Orchestrating ${modelList.length} AI Models`);
    console.log(`📍 Tissue: ${species} ${tissue}`);
    console.log(`🔬 Comparing: ${celltypes.join(', ')}`);
    console.log(`🧬 Markers: ${markerSet}`);
    console.log(`🤖 Models: ${modelList.map(m => modelToPersona[m]?.split(' ').pop() || m.split('/').pop()).join(', ')}`);
    console.log(`🔑 API Key length: ${apiKey ? apiKey.length : 0} characters`);
    console.log(`🌐 Provider: ${provider}`);

    // Construct initial prompt matching Python version
    const celltypesListStr = celltypes.map(ct => `- ${ct}`).join('\n');
    const initialPrompt = `You are a professional biologist. Your task is to analyze how well a given marker set matches a list of cell types from ${species} ${tissue}.

For EACH of the following cell types, you must provide your analysis in a specific structured format.
The cell types to analyze are:
${celltypesListStr}

CRITICAL: You must use EXACTLY this format for EACH cell type (no variations):

<celltype>B cells</celltype>
<reasoning>
Detailed reasoning analyzing the marker genes specifically for B cells. Consider which markers are expressed, which are specific, and overall match quality.
</reasoning>
<score>85</score>

<celltype>T cells</celltype>
<reasoning>
Detailed reasoning analyzing the marker genes specifically for T cells. Consider which markers are expressed, which are specific, and overall match quality.
</reasoning>
<score>92</score>

Replace "B cells" and "T cells" with the actual cell type names. The score must be a number from 0-100.

Please provide a complete block of <celltype>, <reasoning>, and <score> for every cell type listed above.

Marker genes to analyze: ${markerSet}`;
    
    console.log(`\n📋 Generated initial prompt:\n${initialPrompt}\n`);

    // Initialize results storage
    const allResults = [];
    let currentResults = [];
    let roundsPerformed = 0;
    let consensusReached = false;
    let finalConsensus = null;

    // --- Initial Analysis Round ---
    console.log("🎵 Movement I: Initial Analysis (Parallel Processing)");
    console.log(`📋 Initial prompt length: ${initialPrompt.length} characters`);
    console.log(`🤖 Models to query: ${modelList.join(', ')}`);

    const initialPromises = modelList.map(async (model) => {
        const { response, status } = await callModel(model, initialPrompt, apiKey, provider);
        if (status === 'success') {
            const extracted_scores = extractCelltypeScores(response, celltypes);
            return { 
                model, 
                response, 
                extracted_scores, 
                round: 'initial',
                researcher: modelToPersona[model] || model,
                tissue: tissue,
                species: species,
                status: 'success'
            };
        }
        return { 
            model, 
            response, 
            extracted_scores: {}, 
            round: 'initial', 
            researcher: modelToPersona[model] || model,
            tissue: tissue,
            species: species,
            status: 'error'
        };
    });

    const initialResults = await Promise.all(initialPromises);
    console.log(`📈 Initial API calls completed: ${initialResults.length} total`);
    console.log(`✅ Successful calls: ${initialResults.filter(r => r.status === 'success').length}`);
    console.log(`❌ Failed calls: ${initialResults.filter(r => r.status === 'error').length}`);
    
    // DEBUG: Show all raw model responses
    console.log('\n🔍 === RAW MODEL RESPONSES DEBUG ===');
    for (const result of initialResults) {
        if (result.status === 'success') {
            console.log(`\n📋 ${result.model} (${result.researcher}) Response:`);
            console.log('=' * 80);
            console.log(result.response);
            console.log('=' * 80);
            console.log(`🎯 Extracted scores for ${result.model}:`, result.extracted_scores);
        } else {
            console.log(`\n❌ ${result.model} failed:`, result.response);
        }
    }
    console.log('\n🔍 === END RAW RESPONSES ===\n');
    
    currentResults = initialResults.filter(r => r.status === 'success');
    allResults.push(...currentResults);
    
    if (currentResults.length === 0) {
        console.error('💥 No successful API calls - all models failed!');
        throw new Error('All API calls failed. Check your API key and model availability.');
    }

    // Check for initial consensus
    let winners = [];
    const validResults = currentResults.filter(r => r.status === 'success' && r.extracted_scores);
    
    for (const result of validResults) {
        const scores = {};
        for (const [celltype, data] of Object.entries(result.extracted_scores)) {
            try {
                scores[celltype] = parseFloat(data.score);
            } catch (e) {
                scores[celltype] = -1;
            }
        }
        if (Object.keys(scores).length > 0) {
            const winner = Object.keys(scores).reduce((a, b) => scores[a] > scores[b] ? a : b);
            winners.push(winner);
        }
    }

    // Calculate consensus
    if (winners.length > 0) {
        const winnerCounts = {};
        for (const winner of winners) {
            winnerCounts[winner] = (winnerCounts[winner] || 0) + 1;
        }
        const mostCommon = Object.keys(winnerCounts).reduce((a, b) => winnerCounts[a] > winnerCounts[b] ? a : b);
        const consensusRatio = winnerCounts[mostCommon] / validResults.length;
        
        if (consensusRatio >= consensusThreshold) {
            consensusReached = true;
            finalConsensus = mostCommon;
            console.log(`✅ Consensus reached! ${Math.round(consensusRatio * 100)}% agree on: ${finalConsensus}`);
        }
    }

    // --- Discussion Rounds ---
    if (enableDiscussion && !consensusReached && modelList.length > 1 && maxDiscussionRounds > 0) {
        console.log("🎵 Movement II: Discussion & Debate");
        
        for (let roundNum = 0; roundNum < maxDiscussionRounds; roundNum++) {
            if (consensusReached) {
                break;
            }
                
            roundsPerformed = roundNum + 1;
            console.log(`📢 Discussion Round ${roundsPerformed}/${maxDiscussionRounds}`);
            
            // Prepare discussion prompt
            let allResponsesText = "";
            for (const res of currentResults) {
                const researcher = res.researcher || res.model;
                allResponsesText += `\n--- Analysis from ${researcher} ---\n`;
                allResponsesText += `${res.response}\n`;
            }
            
            const discussionPromptTemplate = `You are a professional biologist participating in a panel discussion.
You are {persona_name}. Your colleagues' analyses are provided below. Review their arguments critically.

First, provide a brief critique of each colleague's analysis in a <discussion> block.
Then, provide your refined analysis for each cell type.

Original request:
{original_prompt}

Colleague analyses:
{all_responses}

Start with <discussion>, then provide your analysis for each cell type using <celltype>, <reasoning>, <score> format.`;
            
            // Run discussion round
            const discussionPromises = [];
            const roundName = `discussion_${roundsPerformed}`;
            
            for (const model of modelList) {
                const personaName = modelToPersona[model];
                const thisPrompt = discussionPromptTemplate
                    .replace('{persona_name}', personaName)
                    .replace('{original_prompt}', initialPrompt)
                    .replace('{all_responses}', allResponsesText);
                
                discussionPromises.push(
                    callModel(model, thisPrompt, apiKey, provider).then(({ response, status }) => {
                        if (status === 'success') {
                            const extracted_scores = extractCelltypeScores(response, celltypes);
                            const discussion = extractDiscussion(response);
                            return {
                                model,
                                response,
                                extracted_scores,
                                discussion,
                                round: roundName,
                                researcher: personaName,
                                tissue: tissue,
                                species: species,
                                status: 'success'
                            };
                        }
                        return {
                            model,
                            response,
                            extracted_scores: {},
                            round: roundName,
                            researcher: personaName,
                            tissue: tissue,
                            species: species,
                            status: 'error'
                        };
                    })
                );
            }
            
            const discussionResults = await Promise.all(discussionPromises);
            currentResults = discussionResults.filter(r => r.status === 'success');
            allResults.push(...currentResults);
            
            // Check consensus again
            winners = [];
            const validDiscussionResults = currentResults.filter(r => r.status === 'success' && r.extracted_scores);
            
            for (const result of validDiscussionResults) {
                const scores = {};
                for (const [celltype, data] of Object.entries(result.extracted_scores)) {
                    try {
                        scores[celltype] = parseFloat(data.score);
                    } catch (e) {
                        scores[celltype] = -1;
                    }
                }
                if (Object.keys(scores).length > 0) {
                    const winner = Object.keys(scores).reduce((a, b) => scores[a] > scores[b] ? a : b);
                    winners.push(winner);
                }
            }
            
            if (winners.length > 0) {
                const winnerCounts = {};
                for (const winner of winners) {
                    winnerCounts[winner] = (winnerCounts[winner] || 0) + 1;
                }
                const mostCommon = Object.keys(winnerCounts).reduce((a, b) => winnerCounts[a] > winnerCounts[b] ? a : b);
                const consensusRatio = winnerCounts[mostCommon] / validDiscussionResults.length;
                
                if (consensusRatio >= consensusThreshold) {
                    consensusReached = true;
                    finalConsensus = mostCommon;
                    console.log(`✅ Consensus reached! ${Math.round(consensusRatio * 100)}% agree on: ${finalConsensus}`);
                } else {
                    console.log(`⚡ No consensus yet (${Math.round(consensusRatio * 100)}% for ${mostCommon})`);
                }
            }
        }
    }

    // --- Generate Summary Statistics ---
    const summary = {
        total_rounds: 1 + roundsPerformed,
        models_used: modelList.length,
        consensus_reached: consensusReached,
        consensus_celltype: finalConsensus,
        consensus_confidence: 0.0
    };

    // Calculate final scores and voting information for each cell type
    const celltypeFinalScores = {};
    const celltypeVotes = {};
    
    // First, count votes based on which cell type each model chose
    const finalRound = roundsPerformed > 0 ? `discussion_${roundsPerformed}` : 'initial';
    const finalRoundResults = allResults.filter(r => r.round === finalRound && r.status === 'success');
    
    for (const result of finalRoundResults) {
        if (result.extracted_scores) {
            // Find which cell type this model voted for (highest score)
            let bestCellType = null;
            let bestScore = -1;
            
            for (const [celltype, data] of Object.entries(result.extracted_scores)) {
                try {
                    const score = parseFloat(data.score);
                    if (score > bestScore) {
                        bestScore = score;
                        bestCellType = celltype;
                    }
                } catch (e) {
                    // Skip invalid scores
                }
            }
            
            if (bestCellType) {
                celltypeVotes[bestCellType] = (celltypeVotes[bestCellType] || 0) + 1;
            }
        }
    }
    
    // Calculate statistics for each cell type
    for (const celltype of celltypes) {
        const scores = [];
        
        for (const result of finalRoundResults) {
            if (result.extracted_scores && celltype in result.extracted_scores) {
                try {
                    const score = parseFloat(result.extracted_scores[celltype].score);
                    scores.push(score);
                } catch (e) {
                    // Skip invalid scores
                }
            }
        }
        
        if (scores.length > 0) {
            const mean = scores.reduce((a, b) => a + b, 0) / scores.length;
            const min = Math.min(...scores);
            const max = Math.max(...scores);
            const variance = scores.reduce((acc, score) => acc + Math.pow(score - mean, 2), 0) / scores.length;
            const std = Math.sqrt(variance);
            
            celltypeFinalScores[celltype] = { 
                mean, 
                min, 
                max, 
                std,
                votes: celltypeVotes[celltype] || 0
            };
        }
    }

    summary.celltype_scores = celltypeFinalScores;

    if (finalConsensus && finalConsensus in celltypeFinalScores) {
        summary.consensus_confidence = celltypeFinalScores[finalConsensus].mean / 100.0;
    }

    const reportHtml = generateComparisonHtmlReport(allResults);

    console.log("🎼 Symphony Complete!");
    console.log(`📈 Performance Summary:`);
    console.log(`  • Models: ${summary.models_used} experts consulted`);
    console.log(`  • Rounds: ${summary.total_rounds} total (1 initial + ${roundsPerformed} discussion)`);
    console.log(`  • Consensus: ${consensusReached ? '✅ Yes' : '❌ No'}`);
    if (finalConsensus) {
        console.log(`  • Winner: ${finalConsensus} (confidence: ${Math.round(summary.consensus_confidence * 100)}%)`);
    }

    return {
        results: allResults,
        consensus: finalConsensus,
        confidence: summary.consensus_confidence,
        reportHtml: reportHtml,
        summary: summary,
        raw_responses: allResults.map(r => ({ 
            model: r.model, 
            researcher: r.researcher,
            response: r.response, 
            extracted: r.extracted_scores,
            round: r.round 
        })) // DEBUG: Include raw responses for debugging
    };
} 