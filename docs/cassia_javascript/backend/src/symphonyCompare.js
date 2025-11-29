/**
 * CASSIA Symphony Compare - Multi-Model Cell Type Comparison with AI Consensus Building
 * 100% Python-compatible JavaScript implementation
 * 
 * This module provides an advanced cell type comparison function that orchestrates multiple AI models
 * to analyze and compare cell types based on marker expression patterns. It features:
 * 
 * - Parallel multi-model analysis for diverse perspectives
 * - Automatic consensus detection and discussion rounds when models disagree
 * - Beautiful interactive HTML reports with score progression tracking
 * - Structured CSV output for downstream analysis
 * - Support for custom model configurations
 */

import fs from 'fs';
import path from 'path';

// Import utility functions
try {
    const { callLLM } = await import('./llm_utils.js');
    global.callLLM = callLLM;
} catch (error) {
    console.warn('Warning: Could not import llm_utils.js:', error.message);
}

/**
 * Extract scores and reasoning for each celltype from the LLM response.
 * Replicates Python extract_celltype_scores function.
 * 
 * @param {string} responseText - The raw response text from the LLM
 * @param {Array<string>} celltypes - List of cell types being compared
 * @returns {Object} Dictionary with celltype as key and dict of score/reasoning as value
 */
export function extractCelltypeScores(responseText, celltypes) {
    const results = {};
    
    // Try to find celltype-specific responses
    for (const celltype of celltypes) {
        // Look for celltype tags
        const escapedCelltype = celltype.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
        const celltypePattern = new RegExp(`<celltype>${escapedCelltype}</celltype>(.*?)(?=<celltype>|$)`, 'gis');
        const celltypeMatch = responseText.match(celltypePattern);
        
        let celltypeContent;
        if (celltypeMatch && celltypeMatch.length > 0) {
            // Extract the content part from the full match
            const fullMatch = celltypeMatch[0];
            const contentStart = fullMatch.indexOf('</celltype>') + '</celltype>'.length;
            celltypeContent = fullMatch.substring(contentStart);
        } else {
            // Fallback: look for mentions of the celltype
            celltypeContent = responseText;
        }
        
        // Extract reasoning
        const reasoningPattern = /<reasoning>(.*?)<\/reasoning>/gis;
        const reasoningMatch = celltypeContent.match(reasoningPattern);
        const reasoning = reasoningMatch && reasoningMatch.length > 0 ? reasoningMatch[0].replace(/<\/?reasoning>/g, '').trim() : "No reasoning found";
        
        // Extract score
        const scorePattern = /<score>(\d+(?:\.\d+)?)<\/score>/gi;
        const scoreMatch = celltypeContent.match(scorePattern);
        const score = scoreMatch && scoreMatch.length > 0 ? scoreMatch[0].replace(/<\/?score>/g, '') : "No score found";
        
        results[celltype] = {
            score: score,
            reasoning: reasoning
        };
    }
    
    // If no structured responses found, try to extract any scores mentioned
    const hasValidScores = Object.values(results).some(result => result.score !== "No score found");
    if (!hasValidScores) {
        // Fallback: look for any numbers that might be scores
        for (let i = 0; i < celltypes.length; i++) {
            const celltype = celltypes[i];
            const escapedCelltype = celltype.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
            const pattern = new RegExp(`(${escapedCelltype}[^0-9]*(\\d+(?:\\.\\d+)?))`, 'gi');
            const matches = [...responseText.matchAll(pattern)];
            if (matches.length > 0) {
                const firstMatch = matches[0];
                if (firstMatch.length > 2) {
                    results[celltype].score = firstMatch[2];
                    results[celltype].reasoning = `Extracted from: ${firstMatch[0]}`;
                }
            }
        }
    }
    
    return results;
}

/**
 * Extracts the content of the <discussion> tag from the response.
 * Replicates Python extract_discussion function.
 * 
 * @param {string} responseText - The raw response text from the LLM
 * @returns {string} Extracted discussion content
 */
export function extractDiscussion(responseText) {
    const discussionPattern = /<discussion>(.*?)<\/discussion>/gis;
    const discussionMatch = responseText.match(discussionPattern);
    if (discussionMatch && discussionMatch.length > 0) {
        // Extract content between the tags
        return discussionMatch[0].replace(/<\/?discussion>/gi, '').trim();
    }
    return "No discussion found";
}

/**
 * Helper function to make a single API call to a model.
 * Replicates Python _call_model function.
 * 
 * @param {string} model - Model identifier
 * @param {string} prompt - The prompt to send
 * @param {string} tissue - Tissue type being analyzed
 * @param {string} species - Species being analyzed
 * @param {Array<string>} celltypes - List of cell types
 * @param {string} roundName - Name of the current round
 * @param {string} apiKey - OpenRouter API key
 * @param {boolean} isDiscussionRound - Whether this is a discussion round
 * @returns {Promise<Object>} Result object with model response and extracted data
 */
export async function callModel(model, prompt, tissue, species, celltypes, roundName, apiKey, isDiscussionRound = false) {
    try {
        const response = await fetch("https://openrouter.ai/api/v1/chat/completions", {
            method: "POST",
            headers: {
                "Authorization": `Bearer ${apiKey}`,
                "HTTP-Referer": "https://elliotxie.github.io/CASSIA/",
                "X-Title": "CASSIA",
                "Content-Type": "application/json"
            },
            body: JSON.stringify({
                model: model,
                messages: [{ role: "user", content: prompt }]
            })
        });

        if (response.ok) {
            const responseData = await response.json();
            const modelResponse = responseData.choices[0].message.content;
            const extractedScores = extractCelltypeScores(modelResponse, celltypes);
            
            let discussion = null;
            if (isDiscussionRound) {
                discussion = extractDiscussion(modelResponse);
            }

            console.log(`Model (${roundName}): ${model} [OK]`);

            const resultDict = {
                model: model,
                tissue: tissue,
                species: species,
                cell_types: celltypes.join(', '),
                response: modelResponse,
                extracted_scores: extractedScores,
                status: 'success',
                round: roundName
            };

            if (discussion) {
                resultDict.discussion = discussion;
            }
            return resultDict;
        } else {
            console.log(`Model (${roundName}): ${model} [ERROR ${response.status}]`);
            return {
                model: model,
                tissue: tissue,
                species: species,
                cell_types: celltypes.join(', '),
                response: `Error: ${response.status}`,
                extracted_scores: {},
                status: 'error',
                round: roundName
            };
        }
    } catch (error) {
        console.log(`Model (${roundName}): ${model} [EXCEPTION]`);
        return {
            model: model,
            tissue: tissue,
            species: species,
            cell_types: celltypes.join(', '),
            response: `Exception: ${error.message}`,
            extracted_scores: {},
            status: 'error',
            round: roundName
        };
    }
}

/**
 * Generate an HTML report for cell type comparison results, including discussion progression.
 * Replicates Python generate_comparison_html_report function.
 * 
 * @param {Array<Object>} allResults - All results from model calls
 * @param {string} outputFile - Output file path (optional)
 * @returns {string} Generated HTML content
 */
export function generateComparisonHtmlReport(allResults, outputFile = null) {
    if (!allResults || allResults.length === 0) {
        return "<html><body><h1>No results to display</h1></body></html>";
    }
    
    // Get unique cell types, models, and researchers
    const celltypes = new Set();
    const models = new Set();
    const researchers = new Set();
    const rounds = new Set();
    
    for (const result of allResults) {
        if (result.extracted_scores) {
            Object.keys(result.extracted_scores).forEach(ct => celltypes.add(ct));
        }
        models.add(result.model || 'Unknown');
        researchers.add(result.researcher || result.model || 'Unknown');
        rounds.add(result.round || 'initial');
    }
    
    const celltypesList = Array.from(celltypes).sort();
    const modelsList = Array.from(models).sort();
    const researchersList = Array.from(researchers).sort();
    const roundsList = Array.from(rounds).sort((a, b) => {
        if (a === 'initial') return -1;
        if (b === 'initial') return 1;
        if (a.startsWith('discussion_') && b.startsWith('discussion_')) {
            const aNum = parseInt(a.split('_')[1]) || 0;
            const bNum = parseInt(b.split('_')[1]) || 0;
            return aNum - bNum;
        }
        return a.localeCompare(b);
    });

    // Group results by round
    const roundToResults = {};
    for (const round of roundsList) {
        roundToResults[round] = [];
    }
    for (const result of allResults) {
        const round = result.round || 'initial';
        if (roundToResults[round]) {
            roundToResults[round].push(result);
        }
    }

    // Generate HTML
    let htmlContent = `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CASSIA Symphony Compare Report</title>
    <style>
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f8f9fa;
            color: #333;
        }
        .container {
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 10px;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
            overflow: hidden;
        }
        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            text-align: center;
        }
        .header h1 {
            margin: 0;
            font-size: 2.5em;
            font-weight: 300;
        }
        .header .subtitle {
            margin-top: 10px;
            font-size: 1.2em;
            opacity: 0.9;
        }
        .content {
            padding: 30px;
        }
        .round-selector {
            margin-bottom: 30px;
            text-align: center;
        }
        .round-button {
            display: inline-block;
            margin: 5px;
            padding: 10px 20px;
            background: #e9ecef;
            border: none;
            border-radius: 25px;
            cursor: pointer;
            font-size: 14px;
            transition: all 0.3s;
        }
        .round-button:hover, .round-button.active {
            background: #667eea;
            color: white;
        }
        .summary-section {
            margin-bottom: 40px;
        }
        .summary-title {
            font-size: 1.5em;
            font-weight: bold;
            margin-bottom: 20px;
            color: #495057;
        }
        .score-table {
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 30px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .score-table th {
            background: #667eea;
            color: white;
            padding: 15px 10px;
            text-align: center;
            font-weight: 600;
        }
        .score-table td {
            padding: 12px 10px;
            text-align: center;
            border-bottom: 1px solid #dee2e6;
        }
        .score-table tr:hover {
            background: #f8f9fa;
        }
        .cell-type-name {
            font-weight: bold;
            text-align: left !important;
            background: #f8f9fa;
        }
        .high-score { background: #d4edda; color: #155724; }
        .medium-score { background: #fff3cd; color: #856404; }
        .low-score { background: #f8d7da; color: #721c24; }
        .majority-winner { background: #e7f3ff; }
        .details-title {
            font-size: 1.3em;
            font-weight: bold;
            margin: 30px 0 20px 0;
            color: #495057;
        }
        .round-section {
            margin-bottom: 40px;
            border: 1px solid #dee2e6;
            border-radius: 8px;
            overflow: hidden;
        }
        .round-header {
            background: #f8f9fa;
            padding: 15px 20px;
            font-weight: bold;
            border-bottom: 1px solid #dee2e6;
        }
        .result-item {
            padding: 20px;
            border-bottom: 1px solid #f1f3f4;
        }
        .result-item:last-child {
            border-bottom: none;
        }
        .model-name {
            font-weight: bold;
            color: #667eea;
            margin-bottom: 10px;
        }
        .celltype-result {
            margin: 15px 0;
            padding: 15px;
            background: #f8f9fa;
            border-radius: 5px;
            border-left: 4px solid #667eea;
        }
        .celltype-name {
            font-weight: bold;
            color: #495057;
            margin-bottom: 5px;
        }
        .score {
            font-size: 1.2em;
            font-weight: bold;
            margin-bottom: 10px;
        }
        .reasoning {
            line-height: 1.6;
            color: #6c757d;
        }
        .discussion-section {
            background: #fff3cd;
            border: 1px solid #ffeaa7;
            border-radius: 5px;
            padding: 15px;
            margin-top: 15px;
        }
        .discussion-title {
            font-weight: bold;
            color: #856404;
            margin-bottom: 10px;
        }
        .discussion-content {
            line-height: 1.6;
            color: #856404;
        }
    </style>
    <script>
        function showRound(roundName) {
            // Hide all round sections
            const sections = document.querySelectorAll('.summary-section');
            sections.forEach(section => section.style.display = 'none');
            
            // Show selected round
            const targetSection = document.getElementById('summary_table_' + roundName);
            if (targetSection) {
                targetSection.style.display = 'block';
            }
            
            // Update button states
            const buttons = document.querySelectorAll('.round-button');
            buttons.forEach(btn => btn.classList.remove('active'));
            
            const activeButton = document.querySelector('[onclick="showRound(\\''+roundName+'\\')"]');
            if (activeButton) {
                activeButton.classList.add('active');
            }
        }
    </script>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🎼 CASSIA Symphony Compare</h1>
            <div class="subtitle">Multi-Model Cell Type Comparison Report</div>
            <div class="subtitle">Generated: ${new Date().toLocaleString()}</div>
        </div>
        <div class="content">`;

    // Round selector buttons
    if (roundsList.length > 1) {
        htmlContent += '<div class="round-selector">';
        htmlContent += '<div style="margin-bottom: 15px; font-weight: bold;">View Results by Round:</div>';
        
        for (const roundName of roundsList) {
            const prettyRound = roundName === 'initial' ? 'Initial Analysis' : 
                               roundName.startsWith('discussion_') ? `Discussion Round ${roundName.split('_')[1]}` : 
                               roundName;
            const isActive = roundName === roundsList[roundsList.length - 1];
            const activeClass = isActive ? ' active' : '';
            htmlContent += `<button class="round-button${activeClass}" onclick="showRound('${roundName}')">${prettyRound}</button>`;
        }
        htmlContent += '</div>';
    }

    // Summary tables for each round
    for (const roundName of roundsList) {
        const prettyRound = roundName === 'initial' ? 'Initial Analysis' : 
                           roundName.startsWith('discussion_') ? `Discussion Round ${roundName.split('_')[1]}` : 
                           roundName;
        const isFinal = roundName === roundsList[roundsList.length - 1];
        
        htmlContent += `<div class="summary-section" id="summary_table_${roundName}" style="display: ${isFinal ? 'block' : 'none'};">`;
        htmlContent += `<div class="summary-title">${prettyRound} Summary Table</div>`;
        htmlContent += `
                <table class="score-table">
                    <thead>
                        <tr>
                            <th style="text-align: left;">Cell Type</th>`;
        
        // Get models for this round
        const roundModels = [...new Set(roundToResults[roundName].map(r => r.model || 'Unknown'))].sort();
        for (const model of roundModels) {
            const modelShort = model.includes('/') ? model.split('/').pop() : model;
            htmlContent += `<th>${modelShort}</th>`;
        }
        htmlContent += `
                            <th>Average</th>
                            <th>Majority Votes</th>
                        </tr>
                    </thead>
                    <tbody>`;
        
        // Calculate averages and majority voting for this round
        const celltypeScores = {};
        const modelWinners = [];
        
        for (const celltype of celltypesList) {
            const scores = [];
            for (const model of roundModels) {
                for (const result of roundToResults[roundName]) {
                    if (result.model === model && result.extracted_scores) {
                        if (celltype in result.extracted_scores) {
                            const scoreStr = result.extracted_scores[celltype].score;
                            try {
                                const scoreNum = parseFloat(scoreStr);
                                scores.push(scoreNum);
                            } catch (e) {
                                // Skip invalid scores
                            }
                            break;
                        }
                    }
                }
            }
            celltypeScores[celltype] = scores;
        }
        
        for (const model of roundModels) {
            const modelScores = {};
            for (const result of roundToResults[roundName]) {
                if (result.model === model && result.extracted_scores) {
                    for (const celltype of celltypesList) {
                        if (celltype in result.extracted_scores) {
                            const scoreStr = result.extracted_scores[celltype].score;
                            try {
                                const scoreNum = parseFloat(scoreStr);
                                modelScores[celltype] = scoreNum;
                            } catch (e) {
                                modelScores[celltype] = 0;
                            }
                        }
                    }
                    break;
                }
            }
            if (Object.keys(modelScores).length > 0) {
                const winner = Object.keys(modelScores).reduce((a, b) => modelScores[a] > modelScores[b] ? a : b);
                modelWinners.push(winner);
            }
        }
        
        const majorityVotes = {};
        for (const celltype of celltypesList) {
            majorityVotes[celltype] = modelWinners.filter(w => w === celltype).length;
        }
        
        for (const celltype of celltypesList) {
            const votes = majorityVotes[celltype] || 0;
            const rowClass = votes >= Math.floor(roundModels.length / 2) + 1 ? "majority-winner" : "";
            htmlContent += `<tr class="${rowClass}"><td class="cell-type-name">${celltype}</td>`;
            
            for (const model of roundModels) {
                let score = "N/A";
                for (const result of roundToResults[roundName]) {
                    if (result.model === model && result.extracted_scores) {
                        if (celltype in result.extracted_scores) {
                            score = result.extracted_scores[celltype].score;
                            break;
                        }
                    }
                }
                htmlContent += `<td>${score}</td>`;
            }
            
            const scores = celltypeScores[celltype] || [];
            if (scores.length > 0) {
                const avgScore = scores.reduce((a, b) => a + b, 0) / scores.length;
                let avgClass = "";
                if (avgScore >= 70) {
                    avgClass = "high-score";
                } else if (avgScore >= 40) {
                    avgClass = "medium-score";
                } else {
                    avgClass = "low-score";
                }
                htmlContent += `<td class="${avgClass}"><strong>${avgScore.toFixed(1)}</strong></td>`;
            } else {
                htmlContent += '<td>N/A</td>';
            }
            
            let voteClass = "";
            if (votes >= Math.floor(roundModels.length / 2) + 1) {
                voteClass = "high-score";
            } else if (votes > 0) {
                voteClass = "medium-score";
            } else {
                voteClass = "low-score";
            }
            htmlContent += `<td class="${voteClass}"><strong>${votes}/${roundModels.length}</strong></td>`;
            htmlContent += "</tr>";
        }
        htmlContent += `
                    </tbody>
                </table>
            </div>`;
    }

    // Discussion Progression Section
    htmlContent += '<div class="details-title">🗣️ Discussion Progression by Round</div>';
    for (const roundName of roundsList) {
        const prettyRound = roundName === 'initial' ? 'Initial Analysis' : 
                           roundName.startsWith('discussion_') ? `Discussion Round ${roundName.split('_')[1]}` : 
                           roundName;
        
        htmlContent += `<div class="round-section">`;
        htmlContent += `<div class="round-header">${prettyRound}</div>`;
        
        for (const result of roundToResults[roundName]) {
            htmlContent += '<div class="result-item">';
            const researcher = result.researcher || result.model || 'Unknown';
            htmlContent += `<div class="model-name">🤖 ${researcher}</div>`;
            
            // Show discussion first if it exists
            if (result.discussion && result.discussion !== "No discussion found") {
                htmlContent += '<div class="discussion-section">';
                htmlContent += '<div class="discussion-title">💬 Discussion:</div>';
                htmlContent += `<div class="discussion-content">${result.discussion}</div>`;
                htmlContent += '</div>';
            }
            
            // Show cell type results
            if (result.extracted_scores) {
                for (const celltype of celltypesList) {
                    if (celltype in result.extracted_scores) {
                        const data = result.extracted_scores[celltype];
                        htmlContent += '<div class="celltype-result">';
                        htmlContent += `<div class="celltype-name">${celltype}</div>`;
                        htmlContent += `<div class="score">Score: ${data.score}</div>`;
                        htmlContent += `<div class="reasoning">${data.reasoning}</div>`;
                        htmlContent += '</div>';
                    }
                }
            }
            
            htmlContent += '</div>';
        }
        
        htmlContent += '</div>'; // round-section
    }

    htmlContent += `
        </div>
    </body>
    </html>`;
    
    // Save to file if specified
    if (outputFile) {
        fs.writeFileSync(outputFile, htmlContent, 'utf-8');
        console.log(`HTML report saved to ${outputFile}`);
    }
    
    return htmlContent;
}

/**
 * Symphony Compare - Orchestrate multiple AI models to compare cell types with consensus building.
 * 100% Python-compatible implementation of the main symphonyCompare function.
 * 
 * @param {string} tissue - The tissue type being analyzed
 * @param {Array<string>} celltypes - List of 2-4 cell types to compare
 * @param {string} markerSet - Comma-separated string of gene markers to analyze
 * @param {string} species - Species being analyzed (default: "human")
 * @param {string} modelPreset - Preset model configuration
 * @param {Array<string>} customModels - Custom list of models to use
 * @param {string} outputDir - Directory to save results
 * @param {string} outputBasename - Base name for output files
 * @param {boolean} enableDiscussion - Enable automatic discussion rounds
 * @param {number} maxDiscussionRounds - Maximum discussion rounds to perform
 * @param {number} consensusThreshold - Fraction of models that must agree for consensus
 * @param {boolean} generateReport - Generate interactive HTML report
 * @param {string} apiKey - OpenRouter API key
 * @param {boolean} verbose - Print progress messages
 * @returns {Promise<Object>} Results object with consensus, files, and summary
 */
export async function symphonyCompare(
    tissue,
    celltypes,
    markerSet,
    species = "human",
    modelPreset = "symphony",
    customModels = null,
    outputDir = null,
    outputBasename = null,
    enableDiscussion = true,
    maxDiscussionRounds = 2,
    consensusThreshold = 0.8,
    generateReport = true,
    apiKey = null,
    verbose = true
) {
    // Get API key
    if (apiKey === null) {
        apiKey = process.env.OPENROUTER_API_KEY;
    }
    if (!apiKey) {
        throw new Error("OPENROUTER_API_KEY not found. Set it as an environment variable or pass it as apiKey parameter.");
    }
    
    // Input validation
    if (!celltypes || celltypes.length < 2 || celltypes.length > 4) {
        throw new Error("Please provide 2-4 cell types to compare");
    }
    
    // Set up output directory and filenames
    if (outputDir === null) {
        outputDir = process.cwd();
    }
    fs.mkdirSync(outputDir, { recursive: true });
    
    if (outputBasename === null) {
        const timestamp = new Date().toISOString().replace(/[-:.]/g, '').slice(0, 15);
        const celltypeStr = celltypes.map(ct => ct.replace(/[\s\/]/g, '_')).join('_vs_');
        outputBasename = `symphony_compare_${species}_${tissue.replace(/\s/g, '_')}_${celltypeStr}_${timestamp}`;
    }
    
    const csvFile = path.join(outputDir, `${outputBasename}.csv`);
    const htmlFile = generateReport ? path.join(outputDir, `${outputBasename}_report.html`) : null;
    
    // Define model presets
    const modelPresets = {
        "symphony": [
            "anthropic/claude-3.7-sonnet",
            "openai/o4-mini-high",
            "google/gemini-2.5-pro-preview"
        ],
        "quartet": [
            "anthropic/claude-3.7-sonnet",
            "openai/o4-mini-high", 
            "google/gemini-2.5-pro-preview",
            "meta-llama/llama-3.3-405b"
        ],
        "budget": [
            "google/gemini-2.5-flash",
            "deepseek/deepseek-chat-v3-0324",
            "x-ai/grok-3-mini-beta"
        ]
    };
    
    // Researcher persona names for each model
    const modelPersonas = {
        "google/gemini-2.5-flash": "Dr. Ada Lovelace",
        "deepseek/deepseek-chat-v3-0324": "Dr. Alan Turing", 
        "x-ai/grok-3-mini-beta": "Dr. Marie Curie",
        "anthropic/claude-3.7-sonnet": "Dr. Claude Shannon",
        "openai/o4-mini-high": "Dr. Albert Einstein",
        "google/gemini-2.5-pro-preview": "Dr. Emmy Noether",
        "meta-llama/llama-3.3-405b": "Dr. Rosalind Franklin"
    };
    
    // Select models based on preset or custom list
    let modelList;
    if (modelPreset === "custom" && customModels) {
        modelList = customModels;
    } else if (modelPresets[modelPreset]) {
        modelList = modelPresets[modelPreset];
    } else {
        if (verbose) {
            console.log(`Warning: Unknown preset '${modelPreset}'. Using 'symphony' preset.`);
        }
        modelList = modelPresets["symphony"];
    }
    
    // Get persona names
    const modelToPersona = {};
    for (const m of modelList) {
        modelToPersona[m] = modelPersonas[m] || `Researcher_${m.split('/').pop()}`;
    }
    
    if (verbose) {
        console.log(`\n🎼 CASSIA Symphony Compare - Orchestrating ${modelList.length} AI Models`);
        console.log(`${'='.repeat(60)}`);
        console.log(`📍 Tissue: ${species} ${tissue}`);
        console.log(`🔬 Comparing: ${celltypes.join(', ')}`);
        console.log(`🧬 Markers: ${markerSet}`);
        console.log(`🤖 Models: ${modelList.map(m => modelToPersona[m].split(' ').pop()).join(', ')}`);
        if (enableDiscussion) {
            console.log(`💬 Discussion: Enabled (max ${maxDiscussionRounds} rounds)`);
        }
        console.log(`${'='.repeat(60)}\n`);
    }
    
    // Construct initial prompt
    const celltypesListStr = celltypes.map(ct => `- ${ct}`).join('\n');
    const initialPrompt = `You are a professional biologist. Your task is to analyze how well a given marker set matches a list of cell types from ${species} ${tissue}.

For EACH of the following cell types, you must provide your analysis in a specific structured format.
The cell types to analyze are:
${celltypesListStr}

The required output format for EACH cell type is:
<celltype>cell type name</celltype>
<reasoning>
Your detailed reasoning for the match, considering each marker's relevance.
</reasoning>
<score>A score from 0-100 indicating the match quality.</score>

Please provide a complete block of <celltype>, <reasoning>, and <score> for every cell type listed above.

Ranked marker set: ${markerSet}`;
    
    // Initialize results storage
    const allResults = [];
    let currentResults = [];
    let roundsPerformed = 0;
    let consensusReached = false;
    let finalConsensus = null;
    
    // --- Initial Analysis Round ---
    if (verbose) {
        console.log("🎵 Movement I: Initial Analysis (Parallel Processing)");
    }
    
    const initialPromises = modelList.map(model => 
        callModel(model, initialPrompt, tissue, species, celltypes, 'initial', apiKey, false)
    );
    
    const initialResults = await Promise.all(initialPromises);
    
    for (const result of initialResults) {
        result.researcher = modelToPersona[result.model] || result.model;
        currentResults.push(result);
        allResults.push(result);
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
            if (verbose) {
                console.log(`\n✅ Consensus reached! ${Math.round(consensusRatio * 100)}% agree on: ${finalConsensus}`);
            }
        }
    }
    
    // --- Discussion Rounds ---
    if (enableDiscussion && !consensusReached && modelList.length > 1 && maxDiscussionRounds > 0) {
        if (verbose) {
            console.log(`\n🎵 Movement II: Discussion & Debate`);
        }
        
        for (let roundNum = 0; roundNum < maxDiscussionRounds; roundNum++) {
            if (consensusReached) {
                break;
            }
                
            roundsPerformed = roundNum + 1;
            if (verbose) {
                console.log(`\n  📢 Discussion Round ${roundsPerformed}/${maxDiscussionRounds}`);
            }
            
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
                
                discussionPromises.push(callModel(model, thisPrompt, tissue, species, celltypes, roundName, apiKey, true));
            }
            
            const discussionResults = await Promise.all(discussionPromises);
            
            for (const result of discussionResults) {
                result.researcher = modelToPersona[result.model] || result.model;
                allResults.push(result);
            }
            
            currentResults = discussionResults;
            
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
                    if (verbose) {
                        console.log(`\n  ✅ Consensus reached! ${Math.round(consensusRatio * 100)}% agree on: ${finalConsensus}`);
                    }
                } else if (verbose) {
                    console.log(`  ⚡ No consensus yet (${Math.round(consensusRatio * 100)}% for ${mostCommon})`);
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
    
    // Calculate final scores for each cell type
    const celltypeFinalScores = {};
    for (const celltype of celltypes) {
        const scores = [];
        const finalRound = roundsPerformed > 0 ? `discussion_${roundsPerformed}` : 'initial';
        
        for (const result of allResults.filter(r => r.round === finalRound)) {
            if (result.status === 'success' && result.extracted_scores && celltype in result.extracted_scores) {
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
            
            celltypeFinalScores[celltype] = { mean, min, max, std };
        }
    }
    
    summary.celltype_scores = celltypeFinalScores;
    
    if (finalConsensus && finalConsensus in celltypeFinalScores) {
        summary.consensus_confidence = celltypeFinalScores[finalConsensus].mean / 100.0;
    }
    
    // --- Save Results ---
    if (verbose) {
        console.log(`\n🎵 Movement III: Synthesis & Documentation`);
    }
    
    // Create CSV
    const csvData = [];
    for (const result of allResults) {
        const baseRow = {
            model: result.model,
            researcher: result.researcher || result.model,
            tissue: result.tissue,
            species: result.species, 
            round: result.round || 'initial',
            status: result.status
        };
        
        for (const celltype of celltypes) {
            if (result.extracted_scores && celltype in result.extracted_scores) {
                baseRow[`${celltype}_score`] = result.extracted_scores[celltype].score;
                baseRow[`${celltype}_reasoning`] = result.extracted_scores[celltype].reasoning;
            } else {
                baseRow[`${celltype}_score`] = 'N/A';
                baseRow[`${celltype}_reasoning`] = 'N/A';
            }
        }
        
        baseRow.discussion = result.discussion || 'N/A';
        csvData.push(baseRow);
    }
    
    // Convert to CSV format
    const headers = Object.keys(csvData[0]);
    const csvContent = [
        headers.join(','),
        ...csvData.map(row => 
            headers.map(header => 
                `"${String(row[header]).replace(/"/g, '""')}"`
            ).join(',')
        )
    ].join('\n');
    
    fs.writeFileSync(csvFile, csvContent);
    
    if (verbose) {
        console.log(`  📊 Results saved to: ${csvFile}`);
    }
    
    // Generate HTML report
    if (generateReport) {
        generateComparisonHtmlReport(allResults, htmlFile);
        if (verbose) {
            console.log(`  🎨 Interactive report: ${htmlFile}`);
        }
    }
    
    // --- Final Summary ---
    if (verbose) {
        console.log(`\n🎼 Symphony Complete!`);
        console.log(`${'='.repeat(60)}`);
        console.log(`📈 Performance Summary:`);
        console.log(`  • Models: ${summary.models_used} experts consulted`);
        console.log(`  • Rounds: ${summary.total_rounds} total (1 initial + ${roundsPerformed} discussion)`);
        console.log(`  • Consensus: ${consensusReached ? '✅ Yes' : '❌ No'}`);
        if (finalConsensus) {
            console.log(`  • Winner: ${finalConsensus} (confidence: ${Math.round(summary.consensus_confidence * 100)}%)`);
        }
        console.log(`\n📊 Detailed scores are available in the generated reports.`);
    }
    
    return {
        results: allResults,
        consensus: finalConsensus,
        confidence: summary.consensus_confidence,
        csv_file: csvFile,
        html_file: htmlFile,
        summary: summary,
        dataframe: csvData
    };
}

// Export all functions
export default {
    symphonyCompare,
    extractCelltypeScores,
    extractDiscussion,
    callModel,
    generateComparisonHtmlReport
};