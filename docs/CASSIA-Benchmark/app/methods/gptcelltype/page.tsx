"use client"

import { useState } from "react"
import Link from "next/link"
import { Card, CardContent, CardDescription, CardHeader, CardTitle } from "@/components/ui/card"
import { ModelComparisonChart } from "@/components/model-comparison-chart"
import { TissueSelector } from "@/components/tissue-selector"
import { ArrowLeft, Database, FileText, BarChart2, Beaker } from "lucide-react"

export default function GPTcelltypePage() {
  const [selectedTissue, setSelectedTissue] = useState("All Tissues")

  return (
    <div className="container mx-auto py-10 space-y-8">
      <div>
        <Link href="/" className="flex items-center text-gray-500 hover:text-emerald-600 transition-colors gap-1 mb-4">
          <ArrowLeft className="h-4 w-4" />
          Back to Overview
        </Link>
        <h1 className="text-4xl font-extrabold tracking-tight method-title mb-2">GPTcelltype Benchmark</h1>
        <p className="text-xl text-gray-600">
          Performance comparison of different LLMs using the GPTcelltype method for single-cell annotation
        </p>
      </div>

      <Card className="card-hover border-0 shadow-lg overflow-hidden">
        <CardHeader className="bg-gradient-to-r from-emerald-50 to-emerald-100 border-b">
          <div className="flex flex-col md:flex-row justify-between md:items-center gap-4">
            <div>
              <CardTitle className="text-2xl flex items-center gap-2">
                <BarChart2 className="h-5 w-5 text-emerald-500" />
                Model Performance by Tissue
              </CardTitle>
              <CardDescription className="text-gray-600">
                Comparing different LLM models using the GPTcelltype method across various tissues
              </CardDescription>
            </div>
            <TissueSelector selectedTissue={selectedTissue} onTissueChange={setSelectedTissue} />
          </div>
        </CardHeader>
        <CardContent className="p-6">
          <div className="chart-container bg-white p-4">
            <ModelComparisonChart method="gptcelltype" selectedTissue={selectedTissue} />
          </div>
        </CardContent>
      </Card>

      <div className="grid grid-cols-1 md:grid-cols-2 gap-8">
        <Card className="card-hover border-0 shadow">
          <CardHeader className="bg-gradient-to-r from-emerald-50 to-emerald-100 border-b">
            <CardTitle className="flex items-center gap-2">
              <Beaker className="h-5 w-5 text-emerald-500" />
              Method Description
            </CardTitle>
          </CardHeader>
          <CardContent className="p-6">
            <p className="text-gray-700 leading-relaxed">
              GPTcelltype is the first method that leverages large language models for single-cell annotation. It
              processes all cell types simultaneously, making it fast and outperforming methods like scType and SingleR
              in various benchmark settings.
            </p>
            <p className="text-gray-700 leading-relaxed mt-4">
              The benchmark evaluates the performance of GPT-4o and GPT-4 on 100 cell types across various tissues,
              measuring the accuracy of cell type annotations.
            </p>
          </CardContent>
        </Card>

        <Card className="card-hover border-0 shadow">
          <CardHeader className="bg-gradient-to-r from-emerald-50 to-emerald-100 border-b">
            <CardTitle className="flex items-center gap-2">
              <FileText className="h-5 w-5 text-emerald-500" />
              Benchmark Details
            </CardTitle>
          </CardHeader>
          <CardContent className="p-6 space-y-6">
            <div>
              <h3 className="font-medium text-lg flex items-center gap-2 mb-2">
                <Database className="h-4 w-4 text-emerald-500" />
                Dataset
              </h3>
              <p className="text-gray-700">
                100 cell types across 5 tissues: human kidney, human lung, human large intestine, human fetal skin, and
                whole mouse atlas.
              </p>
            </div>
            <div>
              <h3 className="font-medium text-lg flex items-center gap-2 mb-2">
                <BarChart2 className="h-4 w-4 text-emerald-500" />
                Evaluation Metric
              </h3>
              <p className="text-gray-700">Accuracy score measuring the percentage of correctly annotated cells.</p>
            </div>
            <div>
              <h3 className="font-medium text-lg flex items-center gap-2 mb-2">
                <FileText className="h-4 w-4 text-emerald-500" />
                Models Tested
              </h3>
              <div className="flex flex-wrap gap-2">
                <span className="badge badge-blue">GPT-4o</span>
                <span className="badge badge-green">GPT-4</span>
              </div>
            </div>
          </CardContent>
        </Card>
      </div>
    </div>
  )
}
