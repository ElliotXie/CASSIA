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
 * @param {object} reasoningConfig - Optional reasoning config (OpenAI, OpenRouter, custom HTTP)
 * @param {string} [reasoningConfig.effort] - "high"|"medium"|"low"|"minimal"|"none"
 * @param {number} [reasoningConfig.max_tokens] - Max tokens for reasoning
 * @param {boolean} [reasoningConfig.exclude] - Exclude reasoning from response
 * @param {boolean} [reasoningConfig.enabled] - Enable with default settings
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
    additionalParams = null,
    reasoningConfig = null
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
    
    // OpenAI API call - uses Chat Completions by default, Responses API for reasoning
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

            // Valid OpenAI reasoning effort values (excludes 'none')
            const VALID_REASONING_EFFORTS = ['high', 'medium', 'low'];

            // Only use Responses API if effort is a valid value
            const useResponsesAPI = reasoningConfig
                && reasoningConfig.effort
                && VALID_REASONING_EFFORTS.includes(reasoningConfig.effort.toLowerCase());

            if (useResponsesAPI) {
                // Responses API for reasoning models (GPT-5, o1, etc.)
                const requestOptions = {
                    model,
                    input: apiMessages,
                    reasoning: { effort: reasoningConfig.effort.toLowerCase() },
                    ...additionalParams
                };
                const response = await client.responses.create(requestOptions);
                return response.output_text;
            } else {
                // Chat Completions API (default - more stable, widely supported)
                // Note: Newer OpenAI models (gpt-4o, gpt-4o-mini, etc.) require max_completion_tokens
                const requestOptions = {
                    model,
                    messages: apiMessages,
                    temperature,
                    max_completion_tokens: maxTokens,
                    ...additionalParams
                };
                const response = await client.chat.completions.create(requestOptions);
                return response.choices[0].message.content;
            }
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

            // Check if we're in a browser environment
            const isBrowser = typeof window !== 'undefined';

            if (isBrowser) {
                // Use proxy to bypass CORS in browser
                const requestBody = {
                    baseUrl: provider,
                    apiKey: apiKey,
                    model,
                    messages: apiMessages,
                    temperature,
                    max_tokens: maxTokens,
                    ...additionalParams
                };

                const response = await fetch('/api/proxy', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify(requestBody),
                });

                if (!response.ok) {
                    const errorData = await response.json().catch(() => ({}));
                    throw new Error(errorData.error || `HTTP ${response.status}`);
                }

                const data = await response.json();
                return data.choices[0].message.content;
            } else {
                // Direct call in Node.js environment (no CORS issues)
                const client = new OpenAI({
                    apiKey: apiKey,
                    baseURL: provider,
                });

                // Use Responses API when reasoning is requested, Chat Completions otherwise
                if (reasoningConfig && reasoningConfig.effort) {
                    // Responses API for reasoning models
                    const requestOptions = {
                        model,
                        input: apiMessages,
                        reasoning: { effort: reasoningConfig.effort },
                        ...additionalParams
                    };
                    const response = await client.responses.create(requestOptions);
                    return response.output_text;
                } else {
                    // Chat Completions API (default - more compatible)
                    const requestOptions = {
                        model,
                        messages: apiMessages,
                        temperature,
                        max_tokens: maxTokens,
                        ...additionalParams
                    };
                    const response = await client.chat.completions.create(requestOptions);
                    return response.choices[0].message.content;
                }
            }
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

            // Add beta header if effort is being used (required for effort parameter)
            if (reasoningConfig && reasoningConfig.effort) {
                headers['anthropic-beta'] = 'effort-2025-11-24';
            }

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

            // Add effort config if provided (Anthropic uses output_config.effort)
            if (reasoningConfig && reasoningConfig.effort) {
                messageParams.output_config = { effort: reasoningConfig.effort };
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

            // Add reasoning config if provided (OpenRouter reasoning tokens feature)
            if (reasoningConfig) {
                data.reasoning = {};
                if (reasoningConfig.effort !== undefined) {
                    data.reasoning.effort = reasoningConfig.effort;
                }
                if (reasoningConfig.max_tokens !== undefined) {
                    data.reasoning.max_tokens = reasoningConfig.max_tokens;
                }
                if (reasoningConfig.exclude !== undefined) {
                    data.reasoning.exclude = reasoningConfig.exclude;
                }
                if (reasoningConfig.enabled !== undefined) {
                    data.reasoning.enabled = reasoningConfig.enabled;
                }
            }

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

/**
 * Cheapest models for each provider to use for API key testing.
 * These models are selected to minimize cost while still validating the API key.
 */
const CHEAPEST_MODELS = {
    openai: 'gpt-4o-mini',
    anthropic: 'claude-haiku-4-20250514',
    openrouter: 'google/gemini-2.5-flash',
    custom: 'deepseek-chat'
};

/**
 * Test an API key by sending a minimal request to the provider.
 * Uses the cheapest model available for each provider to minimize cost.
 *
 * @param {string} provider - One of "openai", "anthropic", "openrouter", or "custom"
 * @param {string} apiKey - The API key to test
 * @param {string} customBaseUrl - Base URL for custom provider (required if provider is "custom")
 * @returns {Promise<{success: boolean, error?: string}>} Result of the test
 */
export async function testApiKey(provider, apiKey, customBaseUrl = null) {
    if (!apiKey || !apiKey.trim()) {
        return { success: false, error: 'API key is required' };
    }

    try {
        // Determine the model and provider URL to use
        const model = CHEAPEST_MODELS[provider] || CHEAPEST_MODELS.custom;
        const providerUrl = provider === 'custom' ? customBaseUrl : provider;

        if (provider === 'custom' && !customBaseUrl) {
            return { success: false, error: 'Base URL is required for custom provider' };
        }

        // Send a minimal test request
        await callLLM(
            "say 'a'",      // Simple prompt
            providerUrl,    // Provider or custom URL
            model,          // Cheapest model
            apiKey,         // API key to test
            0.0,            // Temperature (deterministic)
            10,             // Minimal max tokens
            null,           // No system prompt
            null,           // No additional params
            null            // No reasoning config
        );

        return { success: true };
    } catch (error) {
        return { success: false, error: error.message };
    }
}

export default { callLLM, testApiKey };