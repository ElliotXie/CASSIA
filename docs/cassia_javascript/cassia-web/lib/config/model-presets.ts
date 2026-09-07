import { SYMPHONY_PRESETS as _SYMPHONY_PRESETS, MODEL_PERSONAS as _MODEL_PERSONAS } from './model-data';

export interface ModelPreset {
    name: string;
    description: string;
    models: string[];
}

export const MODEL_PRESETS: Record<string, ModelPreset> = {
    premium: {
        name: "Premium",
        description: "Top-tier models for maximum accuracy and reasoning",
        models: [..._SYMPHONY_PRESETS.premium]
    },
    budget: {
        name: "Budget",
        description: "Cost-effective models with strong performance",
        models: [..._SYMPHONY_PRESETS.budget]
    }
};

export const MODEL_PERSONAS: Record<string, string> = { ..._MODEL_PERSONAS };

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
 * Defaults:
 * - GPT-6 Astra: "low" (recommended migration starting point)
 * - GPT-5.6 series: "medium"
 * - Claude Opus: "high"
 *
 * Models WITHOUT reasoning effort (auto or not supported):
 * - GPT-4o, GPT-4, and older OpenAI models: null
 * - Gemini models: null (auto-selects effort internally)
 * - Grok models: null
 * - DeepSeek, Llama, etc.: null
 */
export function getDefaultReasoningEffort(provider: string, model: string): ReasoningEffort | null {
    const modelLower = model.toLowerCase();

    if (modelLower.includes('gpt-6') || modelLower.includes('gpt6')) {
        return 'low';
    }

    // Claude Opus only (direct Anthropic or via OpenRouter) - high
    // Note: Only Opus supports extended thinking/reasoning effort parameter
    if ((provider === 'anthropic' || modelLower.includes('claude')) && modelLower.includes('opus')) {
        return 'high';
    }

    // GPT-5.6 series supports reasoning effort through direct OpenAI and OpenRouter.
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
