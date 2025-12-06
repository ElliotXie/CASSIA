'use client'

import React, { useMemo } from 'react'
import { Calendar, BarChart3, Table, FileText } from 'lucide-react'
import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  CartesianGrid,
  Tooltip,
  ResponsiveContainer,
  Cell
} from 'recharts'
import { ScoreBadge } from './ScoreBadge'
import { cn } from '@/lib/utils'

export interface EvaluationResult {
  goldStandard: string
  prediction: string
  score: number
  reasoning?: string
  tissue?: string
  species?: string
}

export interface EvaluationMetrics {
  meanScore: number
  medianScore: number
  minScore: number
  maxScore: number
  stdScore: number
  count: number
  scoreRatio?: number
  // For discrete 0-5 scale
  perfectRatio?: number
  veryGoodRatio?: number
  goodRatio?: number
  partialRatio?: number
  poorRatio?: number
  nonsensicalRatio?: number
}

export interface EvaluationReportProps {
  data: EvaluationResult[]
  modelName?: string
  isSimilarityScale?: boolean // true for 0-100, false for 0-5
  className?: string
}

/**
 * Calculate metrics from evaluation results
 */
function calculateMetrics(data: EvaluationResult[], isSimilarityScale: boolean): EvaluationMetrics {
  const scores = data.map(d => d.score)
  const count = scores.length

  if (count === 0) {
    return {
      meanScore: 0,
      medianScore: 0,
      minScore: 0,
      maxScore: 0,
      stdScore: 0,
      count: 0,
      scoreRatio: 0
    }
  }

  const sorted = [...scores].sort((a, b) => a - b)
  const mean = scores.reduce((a, b) => a + b, 0) / count
  const median = count % 2 === 0
    ? (sorted[count / 2 - 1] + sorted[count / 2]) / 2
    : sorted[Math.floor(count / 2)]

  const variance = scores.reduce((acc, s) => acc + Math.pow(s - mean, 2), 0) / count
  const std = Math.sqrt(variance)

  const maxPossible = isSimilarityScale ? 100 : 5
  const scoreRatio = (scores.reduce((a, b) => a + b, 0) / (count * maxPossible)) * 100

  const metrics: EvaluationMetrics = {
    meanScore: mean,
    medianScore: median,
    minScore: Math.min(...scores),
    maxScore: Math.max(...scores),
    stdScore: std,
    count,
    scoreRatio
  }

  // Add discrete score ratios for 0-5 scale
  if (!isSimilarityScale) {
    metrics.perfectRatio = scores.filter(s => s === 5).length / count
    metrics.veryGoodRatio = scores.filter(s => s === 4).length / count
    metrics.goodRatio = scores.filter(s => s === 3).length / count
    metrics.partialRatio = scores.filter(s => s === 2).length / count
    metrics.poorRatio = scores.filter(s => s === 1).length / count
    metrics.nonsensicalRatio = scores.filter(s => s === 0).length / count
  }

  return metrics
}

/**
 * Calculate histogram data
 */
function calculateHistogram(data: EvaluationResult[], isSimilarityScale: boolean) {
  const scores = data.map(d => d.score)

  if (isSimilarityScale) {
    // Bins for 0-100 scale
    const bins = [
      { name: '0-20', min: 0, max: 20, count: 0 },
      { name: '20-40', min: 20, max: 40, count: 0 },
      { name: '40-60', min: 40, max: 60, count: 0 },
      { name: '60-80', min: 60, max: 80, count: 0 },
      { name: '80-100', min: 80, max: 100, count: 0 },
    ]

    scores.forEach(score => {
      // Handle score of exactly 100 (include in last bin)
      const bin = bins.find(b => score >= b.min && (score < b.max || (b.max === 100 && score === 100)))
      if (bin) bin.count++
    })

    return bins
  } else {
    // Discrete 0-5 scale
    const counts = [0, 0, 0, 0, 0, 0]
    scores.forEach(score => {
      const idx = Math.round(score)
      if (idx >= 0 && idx <= 5) counts[idx]++
    })

    return [
      { name: '0', count: counts[0] },
      { name: '1', count: counts[1] },
      { name: '2', count: counts[2] },
      { name: '3', count: counts[3] },
      { name: '4', count: counts[4] },
      { name: '5', count: counts[5] },
    ]
  }
}

/**
 * Calculate tissue/species breakdown
 */
function calculateBreakdown(data: EvaluationResult[]) {
  const groups: Record<string, { tissue: string; species: string; scores: number[] }> = {}

  data.forEach(item => {
    const key = `${item.tissue || 'Unknown'}_${item.species || 'Unknown'}`
    if (!groups[key]) {
      groups[key] = { tissue: item.tissue || 'Unknown', species: item.species || 'Unknown', scores: [] }
    }
    groups[key].scores.push(item.score)
  })

  return Object.values(groups).map(g => ({
    tissue: g.tissue,
    species: g.species,
    avgScore: g.scores.reduce((a, b) => a + b, 0) / g.scores.length,
    count: g.scores.length
  }))
}

/**
 * Get bar colors based on score range
 */
function getBarColor(index: number, isSimilarityScale: boolean): string {
  if (isSimilarityScale) {
    const colors = ['#ef4444', '#f97316', '#eab308', '#84cc16', '#22c55e']
    return colors[index]
  }
  const colors = ['#ef4444', '#f97316', '#eab308', '#84cc16', '#22c55e', '#14b8a6']
  return colors[index]
}

/**
 * Evaluation Report component
 * Displays score distribution histogram and statistics
 */
export function EvaluationReport({
  data,
  modelName,
  isSimilarityScale = false,
  className
}: EvaluationReportProps) {
  // Calculate metrics
  const metrics = useMemo(() => calculateMetrics(data, isSimilarityScale), [data, isSimilarityScale])
  const histogram = useMemo(() => calculateHistogram(data, isSimilarityScale), [data, isSimilarityScale])
  const breakdown = useMemo(() => calculateBreakdown(data), [data])

  // Timestamp
  const timestamp = useMemo(() => new Date().toLocaleString(), [])

  // Get sample results for each score range
  const sampleResults = useMemo(() => {
    const samples: EvaluationResult[] = []

    if (isSimilarityScale) {
      const bins = [[0, 20], [20, 40], [40, 60], [60, 80], [80, 100]]
      bins.forEach(([min, max]) => {
        const binData = data.filter(d => d.score >= min && d.score < max)
        samples.push(...binData.slice(0, 3))
      })
    } else {
      for (let score = 0; score <= 5; score++) {
        const scoreData = data.filter(d => Math.round(d.score) === score)
        samples.push(...scoreData.slice(0, 3))
      }
    }

    return samples
  }, [data, isSimilarityScale])

  if (data.length === 0) {
    return (
      <div className={cn("text-center py-12 text-gray-500", className)}>
        <BarChart3 className="h-12 w-12 mx-auto mb-4 text-gray-300" />
        <p>No evaluation data available to display.</p>
      </div>
    )
  }

  return (
    <div className={cn("space-y-6", className)}>
      {/* Header */}
      <header className="bg-gradient-to-r from-indigo-600 to-purple-600 text-white rounded-2xl p-8 shadow-lg">
        <h1 className="text-3xl font-bold mb-2">
          LLM Cell Type Annotation Evaluation Report
          {modelName && <span className="text-indigo-200"> - {modelName}</span>}
        </h1>
        <div className="flex flex-wrap gap-4 mt-4 text-sm">
          <span className="flex items-center gap-2 bg-white/20 px-3 py-1 rounded-lg">
            <Calendar className="h-4 w-4" />
            Generated: {timestamp}
          </span>
          <span className="flex items-center gap-2 bg-white/20 px-3 py-1 rounded-lg">
            <Table className="h-4 w-4" />
            Total Samples: {metrics.count}
          </span>
          <span className="flex items-center gap-2 bg-white/20 px-3 py-1 rounded-lg">
            <BarChart3 className="h-4 w-4" />
            Score Type: {isSimilarityScale ? 'Similarity (0-100)' : 'Discrete (0-5)'}
          </span>
        </div>
      </header>

      {/* Summary Metrics */}
      <div className="bg-white rounded-xl border border-gray-200 p-6 shadow-sm">
        <h2 className="text-xl font-bold text-gray-800 mb-4 flex items-center gap-2">
          <FileText className="h-5 w-5 text-indigo-600" />
          Summary Metrics
        </h2>
        <div className="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-6 gap-4">
          <div className="bg-indigo-50 rounded-lg p-4 text-center">
            <div className="text-2xl font-bold text-indigo-700">{metrics.meanScore.toFixed(2)}</div>
            <div className="text-sm text-gray-600">Mean Score</div>
          </div>
          <div className="bg-indigo-50 rounded-lg p-4 text-center">
            <div className="text-2xl font-bold text-indigo-700">{metrics.medianScore.toFixed(2)}</div>
            <div className="text-sm text-gray-600">Median Score</div>
          </div>
          <div className="bg-indigo-50 rounded-lg p-4 text-center">
            <div className="text-2xl font-bold text-indigo-700">{metrics.minScore.toFixed(2)}</div>
            <div className="text-sm text-gray-600">Min Score</div>
          </div>
          <div className="bg-indigo-50 rounded-lg p-4 text-center">
            <div className="text-2xl font-bold text-indigo-700">{metrics.maxScore.toFixed(2)}</div>
            <div className="text-sm text-gray-600">Max Score</div>
          </div>
          <div className="bg-indigo-50 rounded-lg p-4 text-center">
            <div className="text-2xl font-bold text-indigo-700">{metrics.stdScore.toFixed(2)}</div>
            <div className="text-sm text-gray-600">Std Dev</div>
          </div>
          <div className="bg-green-50 rounded-lg p-4 text-center">
            <div className="text-2xl font-bold text-green-700">{(metrics.scoreRatio || 0).toFixed(1)}%</div>
            <div className="text-sm text-gray-600">Score Ratio</div>
          </div>
        </div>

        {/* Discrete score ratios */}
        {!isSimilarityScale && (
          <div className="mt-4 grid grid-cols-3 md:grid-cols-6 gap-2">
            <div className="bg-red-50 rounded p-2 text-center text-sm">
              <div className="font-semibold text-red-700">{((metrics.nonsensicalRatio || 0) * 100).toFixed(1)}%</div>
              <div className="text-xs text-gray-500">Score 0</div>
            </div>
            <div className="bg-orange-50 rounded p-2 text-center text-sm">
              <div className="font-semibold text-orange-700">{((metrics.poorRatio || 0) * 100).toFixed(1)}%</div>
              <div className="text-xs text-gray-500">Score 1</div>
            </div>
            <div className="bg-yellow-50 rounded p-2 text-center text-sm">
              <div className="font-semibold text-yellow-700">{((metrics.partialRatio || 0) * 100).toFixed(1)}%</div>
              <div className="text-xs text-gray-500">Score 2</div>
            </div>
            <div className="bg-lime-50 rounded p-2 text-center text-sm">
              <div className="font-semibold text-lime-700">{((metrics.goodRatio || 0) * 100).toFixed(1)}%</div>
              <div className="text-xs text-gray-500">Score 3</div>
            </div>
            <div className="bg-green-50 rounded p-2 text-center text-sm">
              <div className="font-semibold text-green-700">{((metrics.veryGoodRatio || 0) * 100).toFixed(1)}%</div>
              <div className="text-xs text-gray-500">Score 4</div>
            </div>
            <div className="bg-teal-50 rounded p-2 text-center text-sm">
              <div className="font-semibold text-teal-700">{((metrics.perfectRatio || 0) * 100).toFixed(1)}%</div>
              <div className="text-xs text-gray-500">Score 5</div>
            </div>
          </div>
        )}
      </div>

      {/* Histogram */}
      <div className="bg-white rounded-xl border border-gray-200 p-6 shadow-sm">
        <h2 className="text-xl font-bold text-gray-800 mb-4 flex items-center gap-2">
          <BarChart3 className="h-5 w-5 text-indigo-600" />
          Score Distribution
        </h2>
        <div className="h-80">
          <ResponsiveContainer width="100%" height="100%">
            <BarChart data={histogram} margin={{ top: 20, right: 30, left: 20, bottom: 5 }}>
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis dataKey="name" />
              <YAxis />
              <Tooltip
                formatter={(value: number) => [value, 'Count']}
                contentStyle={{ borderRadius: 8, border: '1px solid #e5e7eb' }}
              />
              <Bar dataKey="count" radius={[4, 4, 0, 0]}>
                {histogram.map((_, index) => (
                  <Cell key={`cell-${index}`} fill={getBarColor(index, isSimilarityScale)} />
                ))}
              </Bar>
            </BarChart>
          </ResponsiveContainer>
        </div>
      </div>

      {/* Tissue/Species Breakdown */}
      {breakdown.length > 1 && (
        <div className="bg-white rounded-xl border border-gray-200 p-6 shadow-sm">
          <h2 className="text-xl font-bold text-gray-800 mb-4 flex items-center gap-2">
            <Table className="h-5 w-5 text-indigo-600" />
            Breakdown by Tissue and Species
          </h2>
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead>
                <tr className="bg-indigo-600 text-white">
                  <th className="px-4 py-3 text-left font-semibold">Tissue</th>
                  <th className="px-4 py-3 text-left font-semibold">Species</th>
                  <th className="px-4 py-3 text-left font-semibold">Average Score</th>
                  <th className="px-4 py-3 text-left font-semibold">Count</th>
                </tr>
              </thead>
              <tbody>
                {breakdown.map((row, i) => (
                  <tr key={i} className={i % 2 === 0 ? 'bg-white' : 'bg-gray-50'}>
                    <td className="px-4 py-3">{row.tissue}</td>
                    <td className="px-4 py-3">{row.species}</td>
                    <td className="px-4 py-3">
                      <ScoreBadge score={row.avgScore} size="sm" />
                    </td>
                    <td className="px-4 py-3">{row.count}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      )}

      {/* Sample Results */}
      <div className="bg-white rounded-xl border border-gray-200 p-6 shadow-sm">
        <h2 className="text-xl font-bold text-gray-800 mb-4 flex items-center gap-2">
          <Table className="h-5 w-5 text-indigo-600" />
          Sample Results
        </h2>
        <div className="overflow-x-auto">
          <table className="w-full text-sm">
            <thead>
              <tr className="bg-indigo-600 text-white">
                <th className="px-4 py-3 text-left font-semibold">Gold Standard</th>
                <th className="px-4 py-3 text-left font-semibold">Prediction</th>
                <th className="px-4 py-3 text-left font-semibold">Score</th>
                <th className="px-4 py-3 text-left font-semibold">Explanation</th>
              </tr>
            </thead>
            <tbody>
              {sampleResults.map((row, i) => (
                <tr key={i} className={cn(
                  "border-b border-gray-100 hover:bg-indigo-50",
                  i % 2 === 0 ? 'bg-white' : 'bg-gray-50'
                )}>
                  <td className="px-4 py-3 font-medium">{row.goldStandard}</td>
                  <td className="px-4 py-3">{row.prediction}</td>
                  <td className="px-4 py-3">
                    <ScoreBadge score={row.score} size="sm" />
                  </td>
                  <td className="px-4 py-3 text-gray-600 max-w-md">
                    <p className="line-clamp-3">{row.reasoning || 'N/A'}</p>
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
        <p className="text-sm text-gray-500 mt-4 italic">
          Showing up to 3 examples for each {isSimilarityScale ? 'score bin (0-20, 20-40, etc.)' : 'score (0-5)'}.
        </p>
      </div>
    </div>
  )
}

export default EvaluationReport
