// This file would contain utility functions for data management
// In a real application, this would likely connect to an API or database

export type Method = "gptcelltype" | "cassia"
export type Tissue = "Human Kidney" | "Human Lung" | "Human Large Intestine" | "Human Fetal Skin" | "Mouse Atlas"
export type Model = "GPT-4o" | "GPT-4" | "GPT-3.5" | "Claude-3" | "Llama-3" | "Mistral"

export interface BenchmarkResult {
  method: Method
  tissue: Tissue
  model: Model
  score: number
}

// In a real application, this would be a function to fetch data from an API
export async function getBenchmarkData(): Promise<BenchmarkResult[]> {
  // This would be replaced with an actual API call
  return [
    // Sample data would go here
  ]
}

// Function to add new benchmark data
export async function addBenchmarkData(data: BenchmarkResult): Promise<void> {
  // This would be replaced with an actual API call to save data
  console.log("Adding benchmark data:", data)
}

// Function to get available methods
export function getAvailableMethods(): Method[] {
  return ["gptcelltype", "cassia"]
}

// Function to get available tissues
export function getAvailableTissues(): Tissue[] {
  return ["Human Kidney", "Human Lung", "Human Large Intestine", "Human Fetal Skin", "Mouse Atlas"]
}

// Function to get available models
export function getAvailableModels(): Model[] {
  return ["GPT-4o", "GPT-4", "GPT-3.5", "Claude-3", "Llama-3", "Mistral"]
}
