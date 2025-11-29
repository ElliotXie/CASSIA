import modelSettings from '../../public/examples/model_settings.json';

export interface ModelPreset {
    name: string;
    description: string;
    models: string[];
}

function getHighPerformanceModels(): string[] {
    const highCostModels = modelSettings.cost_tiers.high.models;
    return highCostModels.filter(model => {
        const provider = Object.values(modelSettings.providers).find(p => 
            Object.keys(p.models || {}).some(m => 
                m === model || p.models[m].actual_name === model
            )
        );
        return provider;
    });
}

function getBalancedModels(): string[] {
    const lowCostModels = modelSettings.cost_tiers.low.models;
    const veryLowCostModels = modelSettings.cost_tiers.very_low.models;
    return [...lowCostModels, ...veryLowCostModels].filter(model => {
        const provider = Object.values(modelSettings.providers).find(p => 
            Object.keys(p.models || {}).some(m => 
                m === model || p.models[m].actual_name === model
            )
        );
        return provider;
    });
}

export const MODEL_PRESETS: Record<string, ModelPreset> = {
    performance: {
        name: "Performance",
        description: "High-performance models for maximum accuracy",
        models: getHighPerformanceModels()
    },
    balance: {
        name: "Balance", 
        description: "Balanced selection of capable models",
        models: getBalancedModels()
    }
};

export const DEFAULT_MODEL_PRESET = "performance";

export function getModelPreset(presetName: string): ModelPreset | null {
    return MODEL_PRESETS[presetName] || null;
}

export function getAvailablePresets(): ModelPreset[] {
    return Object.values(MODEL_PRESETS);
}