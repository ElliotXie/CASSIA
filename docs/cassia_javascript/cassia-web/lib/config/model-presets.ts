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