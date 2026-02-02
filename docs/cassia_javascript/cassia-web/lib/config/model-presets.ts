export interface ModelPreset {
    name: string;
    description: string;
    models: string[];
}

export const MODEL_PRESETS: Record<string, ModelPreset> = {
    premium: {
        name: "Premium",
        description: "Top-tier models for maximum accuracy and reasoning",
        models: [
            "google/gemini-3-pro-preview",
            "anthropic/claude-sonnet-4.5",
            "openai/gpt-5.2",
            "x-ai/grok-4"
        ]
    },
    budget: {
        name: "Budget",
        description: "Cost-effective models with strong performance",
        models: [
            "deepseek/deepseek-v3.2",
            "x-ai/grok-4-fast",
            "moonshotai/kimi-k2-thinking",
            "google/gemini-3-flash-preview"
        ]
    }
};

export const MODEL_PERSONAS: Record<string, string> = {
    "google/gemini-3-pro-preview": "Dr. Emmy Noether",
    "anthropic/claude-sonnet-4.5": "Dr. Claude Shannon",
    "openai/gpt-5.2": "Dr. Albert Einstein",
    "x-ai/grok-4": "Dr. Marie Curie",
    "deepseek/deepseek-v3.2": "Dr. Alan Turing",
    "x-ai/grok-4-fast": "Dr. Nikola Tesla",
    "moonshotai/kimi-k2-thinking": "Dr. Ada Lovelace",
    "google/gemini-3-flash-preview": "Dr. Rosalind Franklin"
};

export const DEFAULT_MODEL_PRESET = "premium";

export function getModelPreset(presetName: string): ModelPreset | null {
    return MODEL_PRESETS[presetName] || null;
}

export function getAvailablePresets(): ModelPreset[] {
    return Object.values(MODEL_PRESETS);
}

export function getModelPersona(modelId: string): string {
    return MODEL_PERSONAS[modelId] || "Research Assistant";
}

// Reasoning effort configuration
export type ReasoningEffort = 'high' | 'medium' | 'low' | 'none';

/**
 * Get the default reasoning effort for a model based on provider and model type.
 *
 * Direct OpenAI provider:
 * - Always returns null (requires identity verification for reasoning models)
 * - User must opt-in via identity verification checkbox in UI
 *
 * OpenRouter/Anthropic providers:
 * - Claude models: "high"
 * - GPT-5 series via OpenRouter: "medium"
 *
 * Models WITHOUT reasoning effort (auto or not supported):
 * - GPT-4o, GPT-4, and older OpenAI models: null
 * - Gemini models: null (auto-selects effort internally)
 * - Grok models: null
 * - DeepSeek, Llama, etc.: null
 */
export function getDefaultReasoningEffort(provider: string, model: string): ReasoningEffort | null {
    const modelLower = model.toLowerCase();

    // Direct OpenAI provider: always return null (use Chat Completions API)
    // OpenAI requires identity verification for reasoning models
    // User must opt-in via identity verification checkbox
    if (provider === 'openai') {
        return null;
    }

    // Claude Opus only (direct Anthropic or via OpenRouter) - high
    // Note: Only Opus supports extended thinking/reasoning effort parameter
    if ((provider === 'anthropic' || modelLower.includes('claude')) && modelLower.includes('opus')) {
        return 'high';
    }

    // GPT-5 series via OpenRouter supports reasoning effort
    if (modelLower.includes('gpt-5') || modelLower.includes('gpt5')) {
        return 'medium';
    }

    // All other models: no reasoning effort configuration
    // - GPT-4o, GPT-4: no reasoning config
    // - Gemini: auto-selects effort internally
    // - Grok: no reasoning config
    // - DeepSeek, Llama, etc.: no reasoning config
    return null;
}

/**
 * Check if a model supports reasoning/effort configuration.
 */
export function modelSupportsReasoning(provider: string, model: string): boolean {
    return getDefaultReasoningEffort(provider, model) !== null;
}

/**
 * Get all available reasoning effort options.
 */
export function getReasoningEffortOptions(): { value: ReasoningEffort; label: string; description: string }[] {
    return [
        { value: 'high', label: 'High', description: 'Most thorough reasoning, slower' },
        { value: 'medium', label: 'Medium', description: 'Balanced reasoning and speed' },
        { value: 'low', label: 'Low', description: 'Faster with less reasoning' },
        { value: 'none', label: 'None', description: 'No extended reasoning' }
    ];
}