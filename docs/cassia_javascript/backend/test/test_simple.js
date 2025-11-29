import { callLLM, runCASSIA, runCASSIABatch, setApiKey } from '../index.js';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

// ⚠️ SET YOUR API KEY HERE FOR TESTING ⚠️
const YOUR_API_KEY = "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"; // Put your OpenRouter API key here

async function quickTest() {
    console.log("🧪 CASSIA JavaScript Quick Test");
    console.log("===============================");
    
    // Check API key
    if (!YOUR_API_KEY && !process.env.OPENROUTER_API_KEY) {
        console.log("❌ No API key found!");
        console.log("Please either:");
        console.log("1. Set YOUR_API_KEY in this file, or");
        console.log("2. Set environment variable: export OPENROUTER_API_KEY=your_key");
        return;
    }
    
    // Set API key if provided
    if (YOUR_API_KEY) {
        setApiKey(YOUR_API_KEY, "openrouter");
        console.log("✅ API key set from file");
    } else {
        console.log("✅ Using API key from environment");
    }
    
    try {
        // Test 1: Simple LLM call
        console.log("\n📋 Test 1: Testing callLLM...");
        const response = await callLLM(
            "What are T cells? Answer in one sentence.",
            "openrouter", 
            "google/gemini-2.5-flash-preview",
            null,
            0.1,
            50
        );
        console.log("Response:", response.substring(0, 100) + "...");
        console.log("✅ Test 1 passed");
        
        // Test 2: Single cluster analysis
        console.log("\n📋 Test 2: Testing runCASSIA...");
        const markers = ["CD3D", "CD3E", "CD8A"];
        console.log(`Analyzing markers: ${markers.join(', ')}`);
        
        const [result, history] = await runCASSIA(
            "google/gemini-2.5-flash-preview",
            0,
            markers,
            "blood",
            "human",
            null,
            "openrouter",
            "v1"
        );
        
        console.log("Result:", result.main_cell_type);
        console.log("Iterations:", result.iterations);
        console.log("✅ Test 2 passed");
        
        console.log("\n🎉 All tests passed! CASSIA is working correctly.");
        
    } catch (error) {
        console.error("❌ Test failed:", error.message);
        if (error.message.includes("401") || error.message.includes("API key")) {
            console.log("This looks like an API key issue. Please check your key.");
        }
    }
}

// Run the test
quickTest();