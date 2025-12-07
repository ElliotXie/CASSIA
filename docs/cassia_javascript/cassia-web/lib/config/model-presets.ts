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
            "openai/gpt-5.1",
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
            "google/gemini-2.5-flash"
        ]
    }
};

export const MODEL_PERSONAS: Record<string, string> = {
    "google/gemini-3-pro-preview": "Dr. Emmy Noether",
    "anthropic/claude-sonnet-4.5": "Dr. Claude Shannon",
    "openai/gpt-5.1": "Dr. Albert Einstein",
    "x-ai/grok-4": "Dr. Marie Curie",
    "deepseek/deepseek-v3.2": "Dr. Alan Turing",
    "x-ai/grok-4-fast": "Dr. Nikola Tesla",
    "moonshotai/kimi-k2-thinking": "Dr. Ada Lovelace",
    "google/gemini-2.5-flash": "Dr. Rosalind Franklin"
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
 * - Anthropic models (direct or via OpenRouter): "high"
 * - OpenAI models (GPT-5, o1): "medium"
 * - Other models with reasoning support (Gemini, DeepSeek): "high"
 * - Models without reasoning support: null
 */
export function getDefaultReasoningEffort(provider: string, model: string): ReasoningEffort | null {
    const modelLower = model.toLowerCase();

    // Anthropic models (direct or via OpenRouter) - high
    if (provider === 'anthropic' || modelLower.includes('claude')) {
        return 'high';
    }

    // OpenAI models (direct or via OpenRouter) - medium
    if (provider === 'openai' || modelLower.includes('gpt-5') || modelLower.includes('o1-') || modelLower.includes('o3-')) {
        return 'medium';
    }

    // Other models with reasoning support - high
    if (modelLower.includes('gemini') ||
        modelLower.includes('deepseek') ||
        modelLower.includes('grok') ||
        modelLower.includes('thinking') ||
        modelLower.includes('kimi')) {
        return 'high';
    }

    // Default: no reasoning
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