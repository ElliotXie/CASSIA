"use client"

import { useState } from "react"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { ModelComparisonChart } from "@/components/model-comparison-chart"
import { ModelRankingChart } from "@/components/model-ranking-chart"
import { TissueSelector } from "@/components/tissue-selector"
import { Database, FileText, BarChart2, Beaker, Medal } from "lucide-react"

export default function CassiaPage() {
  const [selectedTissue, setSelectedTissue] = useState("All Tissues")

  return (
    <div className="container mx-auto py-10 space-y-8">
      <div>
        <h1 className="text-4xl font-extrabold tracking-tight cassia-title mb-2">CASSIA Benchmark</h1>
        <p className="text-xl text-gray-600">
          Performance comparison of different LLMs using the CASSIA method for single-cell annotation
        </p>
      </div>

      <Card className="card-hover border-0 shadow-lg overflow-hidden">
        <CardHeader className="bg-gradient-to-r from-blue-50 to-blue-100 border-b">
          <div className="flex flex-col md:flex-row justify-between md:items-center gap-4">
            <div>
              <CardTitle className="text-2xl flex items-center gap-2">
                <BarChart2 className="h-5 w-5 text-blue-500" />
                Model Performance by Tissue
              </CardTitle>
              <CardDescription className="text-gray-600">
                Comparing different LLM models using the CASSIA method across various tissues
              </CardDescription>
            </div>
            <TissueSelector selectedTissue={selectedTissue} onTissueChange={setSelectedTissue} />
          </div>
        </CardHeader>
        <CardContent className="p-6">
          <div className="chart-container bg-white p-4">
            <ModelComparisonChart method="cassia" selectedTissue={selectedTissue} />
          </div>
        </CardContent>
      </Card>

      <Card className="card-hover border-0 shadow-lg overflow-hidden">
        <CardHeader className="bg-gradient-to-r from-blue-50 to-blue-100 border-b">
          <div className="flex flex-col md:flex-row justify-between md:items-center gap-4">
            <div>
              <CardTitle className="text-2xl flex items-center gap-2">
                <Medal className="h-5 w-5 text-amber-500" />
                Cost-Performance Analysis
              </CardTitle>
              <CardDescription className="text-gray-600">
                Performance vs. cost ranking of models - optimal models appear in the top-left corner
              </CardDescription>
            </div>
          </div>
        </CardHeader>
        <CardContent className="p-6">
          <div className="chart-container bg-white p-4">
            <ModelRankingChart />
          </div>
        </CardContent>
      </Card>

      <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
        <Card className="card-hover border-0 shadow">
          <CardHeader className="bg-gradient-to-r from-blue-50 to-blue-100 border-b">
            <CardTitle className="flex items-center gap-2">
              <Beaker className="h-5 w-5 text-blue-500" />
              Method Description
            </CardTitle>
          </CardHeader>
          <CardContent className="p-6">
            <p className="text-gray-700 leading-relaxed">
              CASSIA (Collective Agent System for Single-cell Interpretable Annotation) is the first multi-agent LLM-based method
              for single-cell annotation. It enhances annotation accuracy across diverse datasets and rare cell types by
              integrating step-by-step reasoning, validation, quality scoring, and optional refinement or
              retrieval-augmented generation.
            </p>
            <p className="text-gray-700 leading-relaxed mt-4">
              The method leverages the collaboration of five basic agents and five advanced agents to provide
              comprehensive and interpretable cell type annotations with robust performance across different tissues.
            </p>
            <div className="mt-6 flex gap-4">
              <a 
                href="https://github.com/ElliotXie/CASSIA" 
                target="_blank" 
                rel="noopener noreferrer" 
                className="inline-flex items-center gap-2 px-4 py-2 bg-gray-100 hover:bg-gray-200 text-gray-800 rounded-md transition-colors"
              >
                <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" className="lucide lucide-github"><path d="M15 22v-4a4.8 4.8 0 0 0-1-3.5c3 0 6-2 6-5.5.08-1.25-.27-2.48-1-3.5.28-1.15.28-2.35 0-3.5 0 0-1 0-3 1.5-2.64-.5-5.36-.5-8 0C6 2 5 2 5 2c-.3 1.15-.3 2.35 0 3.5A5.403 5.403 0 0 0 4 9c0 3.5 3 5.5 6 5.5-.39.49-.68 1.05-.85 1.65-.17.6-.22 1.23-.15 1.85v4"></path><path d="M9 18c-4.51 2-5-2-7-2"></path></svg>
                GitHub Repository
              </a>
              <a 
                href="https://www.biorxiv.org/content/10.1101/2024.12.04.626476v2" 
                target="_blank" 
                rel="noopener noreferrer" 
                className="inline-flex items-center gap-2 px-4 py-2 bg-blue-100 hover:bg-blue-200 text-blue-800 rounded-md transition-colors"
              >
                <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" className="lucide lucide-file-text"><path d="M14.5 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V7.5L14.5 2z"></path><polyline points="14 2 14 8 20 8"></polyline><line x1="16" y1="13" x2="8" y2="13"></line><line x1="16" y1="17" x2="8" y2="17"></line><line x1="10" y1="9" x2="8" y2="9"></line></svg>
                View Preprint
              </a>
              <a 
                href="https://documentationeng.vercel.app/" 
                target="_blank" 
                rel="noopener noreferrer" 
                className="inline-flex items-center gap-2 px-4 py-2 bg-green-100 hover:bg-green-200 text-green-800 rounded-md transition-colors"
              >
                <svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" className="lucide lucide-book-open"><path d="M2 3h6a4 4 0 0 1 4 4v14a3 3 0 0 0-3-3H2z"></path><path d="M22 3h-6a4 4 0 0 0-4 4v14a3 3 0 0 1 3-3h7z"></path></svg>
                Documentation
              </a>
            </div>
          </CardContent>
        </Card>

        <Card className="card-hover border-0 shadow">
          <CardHeader className="bg-gradient-to-r from-blue-50 to-blue-100 border-b">
            <CardTitle className="flex items-center gap-2">
              <FileText className="h-5 w-5 text-blue-500" />
              Benchmark Details
            </CardTitle>
          </CardHeader>
          <CardContent className="p-6 space-y-6">
            <div>
              <h3 className="font-medium text-lg flex items-center gap-2 mb-2">
                <Database className="h-4 w-4 text-blue-500" />
                Dataset
              </h3>
              <p className="text-gray-700">
                100 cell types across 5 tissues: human kidney, human lung, human large intestine, human fetal skin, and
                whole mouse atlas.
              </p>
            </div>
            <div>
              <h3 className="font-medium text-lg flex items-center gap-2 mb-2">
                <BarChart2 className="h-4 w-4 text-blue-500" />
                Evaluation Metric
              </h3>
              <p className="text-gray-700">We built an agent that scores annotations by averaging similarity between predicted and gold standard cell types per tissue. The agent tends to underestimate accuracy, and although some clear errors in the gold standard were corrected, the true accuracy is still considered to be higher.</p>
            </div>
            <div>
              <h3 className="font-medium text-lg flex items-center gap-2 mb-2">
                <FileText className="h-4 w-4 text-blue-500" />
                Models Tested
              </h3>
              <div className="flex flex-wrap gap-2">
                <span className="badge badge-indigo">Llama 4 Maverick</span>
                <span className="badge badge-sky">GPT-4.1</span>
                <span className="badge badge-violet">Claude 3.7</span>
                <span className="badge badge-amber">Gemini 2.5 pro</span>
                <span className="badge badge-yellow">Gemini 2.5 flash</span>
                <span className="badge badge-blue">GPT-O4 Mini High</span>
                <span className="badge badge-pink">Deepseek v3</span>
                <span className="badge badge-teal">QWEN3-235b</span>
              </div>
            </div>
          </CardContent>
        </Card>
      </div>
    </div>
  )
}
