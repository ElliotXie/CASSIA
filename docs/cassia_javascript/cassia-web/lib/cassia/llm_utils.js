import OpenAI from 'openai';
import axios from 'axios';

/**
 * Call an LLM from various providers and return the generated text.
 * 
 * @param {string} prompt - The user prompt to send to the LLM
 * @param {string} provider - One of "openai", "anthropic", or "openrouter", or custom HTTP URL
 * @param {string} model - Specific model from the provider to use (e.g., "gpt-4" for OpenAI)
 * @param {string} apiKey - API key for the provider (if null, gets from environment)
 * @param {number} temperature - Sampling temperature (0-1)
 * @param {number} maxTokens - Maximum tokens to generate
 * @param {string} systemPrompt - Optional system prompt for providers that support it
 * @param {object} additionalParams - Additional parameters to pass to the provider's API
 * @returns {Promise<string>} The generated text response
 */
export async function callLLM(
    prompt,
    provider = "openai",
    model = null,
    apiKey = null,
    temperature = 0.7,
    maxTokens = 7000,
    systemPrompt = null,
    additionalParams = null
) {
    provider = provider.toLowerCase();
    additionalParams = additionalParams || {};
    
    // Default models for each provider if not specified
    const defaultModels = {
        "openai": "gpt-3.5-turbo",
        "anthropic": "claude-3-sonnet-20240229",
        "openrouter": "google/gemini-2.5-flash",
    };
    
    // Use default model if not specified
    if (!model) {
        model = defaultModels[provider];
        if (!model) {
            throw new Error(`No model specified and no default available for provider: ${provider}`);
        }
    }
    
    // API key is required in browser environment
    if (!apiKey) {
        throw new Error(`API key is required for provider: ${provider}`);
    }
    
    // Prepare messages format
    let messages = [];
    if (systemPrompt) {
        messages.push({ role: "system", content: systemPrompt });
    }
    messages.push({ role: "user", content: prompt });
    
    // OpenAI API call
    if (provider === "openai") {
        try {
            const client = new OpenAI({ 
                apiKey,
                dangerouslyAllowBrowser: true 
            });
            
            // Handle message history properly for OpenAI
            let apiMessages = [...messages];
            
            // If additional_params contains message history, merge it properly
            if ('messages' in additionalParams) {
                // Use the full conversation history from additional_params instead
                const historyMessages = { ...additionalParams }.messages;
                delete additionalParams.messages;
                apiMessages = historyMessages;
                
                // Only add system prompt if it's not already in the history
                if (systemPrompt && !apiMessages.some(msg => msg.role === 'system')) {
                    apiMessages.unshift({ role: "system", content: systemPrompt });
                }
            }
            
            const response = await client.chat.completions.create({
                model,
                messages: apiMessages,
                temperature,
                max_tokens: maxTokens,
                ...additionalParams
            });
            
            return response.choices[0].message.content;
        } catch (error) {
            throw new Error(`OpenAI API error: ${error.message}`);
        }
    }
    
    // Custom OpenAI-compatible API call (base_url as provider)
    else if (provider.startsWith("http")) {
        try {
            if (!apiKey) {
                throw new Error("API key is required for custom endpoint");
            }
            
            const client = new OpenAI({ 
                apiKey: apiKey, 
                baseURL: provider,
                dangerouslyAllowBrowser: true 
            });
            
            // Handle message history properly
            let apiMessages = [...messages];
            
            // If additional_params contains message history, merge it properly
            if ('messages' in additionalParams) {
                // Use the full conversation history from additional_params instead
                const historyMessages = { ...additionalParams }.messages;
                delete additionalParams.messages;
                apiMessages = historyMessages;
                
                // Only add system prompt if it's not already in the history
                if (systemPrompt && !apiMessages.some(msg => msg.role === 'system')) {
                    apiMessages.unshift({ role: "system", content: systemPrompt });
                }
            }
            
            // Call the API with the proper message history
            const response = await client.chat.completions.create({
                model,
                messages: apiMessages,
                temperature,
                max_tokens: maxTokens,
                ...additionalParams
            });
            
            return response.choices[0].message.content;
        } catch (error) {
            throw new Error(`Custom API error: ${error.message}`);
        }
    }
    
    // Anthropic API call
    else if (provider === "anthropic") {
        try {
            // Using direct HTTP call to match Python implementation exactly
            const headers = {
                'Content-Type': 'application/json',
                'x-api-key': apiKey,
                'anthropic-version': '2023-06-01'
            };
            
            // Handle message history properly for Anthropic
            let apiMessages = messages;
            
            // If additional_params contains message history, use it
            if ('messages' in additionalParams) {
                const historyMessages = { ...additionalParams }.messages;
                delete additionalParams.messages;
                apiMessages = historyMessages;
            }
            
            // Convert messages to Anthropic format (text content structure)
            const anthropicMessages = apiMessages.map(msg => ({
                role: msg.role === 'system' ? 'user' : msg.role, // Anthropic doesn't have system role in messages
                content: [{ type: "text", text: msg.content }]
            })).filter(msg => msg.role !== 'user' || msg.content[0].text.trim() !== ''); // Remove empty user messages
            
            // Create the message params
            const messageParams = {
                model,
                max_tokens: maxTokens,
                temperature,
                messages: anthropicMessages
            };
            
            // Add system prompt if provided (Anthropic uses separate system field)
            if (systemPrompt) {
                messageParams.system = systemPrompt;
            }
            
            // Add any additional parameters
            Object.assign(messageParams, additionalParams);
            
            // Call the API
            const response = await axios.post(
                'https://api.anthropic.com/v1/messages',
                messageParams,
                { headers }
            );
            
            // Extract the text content from the response
            if (response.data.content && response.data.content.length > 0) {
                const contentBlock = response.data.content[0];
                if (contentBlock.text) {
                    return contentBlock.text;
                } else if (typeof contentBlock === 'object' && contentBlock.text) {
                    return contentBlock.text;
                } else {
                    return String(response.data.content);
                }
            } else {
                return "No content returned from Anthropic API";
            }
        } catch (error) {
            throw new Error(`Anthropic API error: ${error.message}`);
        }
    }
    
    // OpenRouter API call
    else if (provider === "openrouter") {
        try {
            const url = "https://openrouter.ai/api/v1/chat/completions";
            
            const headers = {
                "Authorization": `Bearer ${apiKey}`,
                "Content-Type": "application/json"
            };
            
            // Handle message history properly for OpenRouter
            let apiMessages = [...messages];
            
            // If additional_params contains message history, merge it properly
            if ('messages' in additionalParams) {
                // Use the full conversation history from additional_params instead
                const historyMessages = { ...additionalParams }.messages;
                delete additionalParams.messages;
                apiMessages = historyMessages;
                
                // Only add system prompt if it's not already in the history
                if (systemPrompt && !apiMessages.some(msg => msg.role === 'system')) {
                    apiMessages.unshift({ role: "system", content: systemPrompt });
                }
            }
            
            const data = {
                model,
                messages: apiMessages,
                temperature,
                max_tokens: maxTokens,
                ...additionalParams
            };
            
            const response = await axios.post(url, data, { headers });
            
            return response.data.choices[0].message.content;
        } catch (error) {
            if (error.response) {
                throw new Error(`OpenRouter API error: ${error.response.status} - ${error.response.data?.error?.message || error.message}`);
            } else {
                throw new Error(`OpenRouter API error: ${error.message}`);
            }
        }
    }
    
    else {
        throw new Error(`Unsupported provider: ${provider}`);
    }
}

export default { callLLM };