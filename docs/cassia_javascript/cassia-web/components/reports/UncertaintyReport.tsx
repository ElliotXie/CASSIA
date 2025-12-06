'use client'

import React, { useState, useMemo } from 'react'
import { Calendar, Gauge, PieChart as PieChartIcon, ChevronDown, ChevronUp, AlertTriangle, CheckCircle } from 'lucide-react'
import {
  PieChart,
  Pie,
  Cell,
  ResponsiveContainer,
  Tooltip,
  Legend
} from 'recharts'
import { Button } from '@/components/ui/button'
import { ReportFilters, useReportFilters } from './ReportFilters'
import { ScoreBadge } from './ScoreBadge'
import { cn } from '@/lib/utils'

export interface RoundResult {
  iteration: number
  mainType: string
  subType: string
  agreement?: 'full' | 'partial' | 'none'
}

export interface ClusterUncertainty {
  clusterId: string
  consensusScore: number // 0-100 (llm_generated_consensus_score_llm in Python)
  consensusMainType: string // general_celltype_llm in Python
  consensusSubType: string // sub_celltype_llm in Python
  possibleMixedTypes?: string[]
  llmReasoning?: string
  // Original per-iteration results (before unification)
  originalResults?: RoundResult[]
  // Unified results from LLM unification
  unifiedResultsLlm?: RoundResult[]
  // Legacy field for backward compatibility (maps to originalResults)
  roundResults?: RoundResult[]
}

export interface UncertaintyReportProps {
  data: ClusterUncertainty[]
  title?: string
  className?: string
}

// Colors for pie charts
const PIE_COLORS = [
  '#14b8a6', '#0d9488', '#06b6d4', '#0ea5e9', '#3b82f6',
  '#6366f1', '#8b5cf6', '#a855f7', '#d946ef', '#ec4899'
]

/**
 * Get color and label based on consensus score
 */
function getScoreInfo(score: number): { color: string; label: string } {
  if (score >= 80) return { color: '#22c55e', label: 'High Confidence' }
  if (score >= 50) return { color: '#eab308', label: 'Moderate Confidence' }
  return { color: '#ef4444', label: 'Low Confidence' }
}

/**
 * Calculate type distribution from round results
 */
function calculateTypeDistribution(results: RoundResult[], field: 'mainType' | 'subType') {
  const counts: Record<string, number> = {}

  results.forEach(r => {
    const value = r[field] || 'Unknown'
    counts[value] = (counts[value] || 0) + 1
  })

  return Object.entries(counts)
    .map(([name, value]) => ({ name, value }))
    .sort((a, b) => b.value - a.value)
}

/**
 * Single cluster uncertainty card
 */
function ClusterUncertaintyCard({
  cluster,
  expanded,
  onToggle
}: {
  cluster: ClusterUncertainty
  expanded: boolean
  onToggle: () => void
}) {
  const scoreInfo = getScoreInfo(cluster.consensusScore)
  // Use originalResults with fallback to roundResults for backward compatibility
  const results = cluster.originalResults || cluster.roundResults || []
  const mainTypeDistribution = useMemo(
    () => calculateTypeDistribution(results, 'mainType'),
    [results]
  )
  const subTypeDistribution = useMemo(
    () => calculateTypeDistribution(results, 'subType'),
    [results]
  )

  // Count agreements
  const agreements = useMemo(() => {
    let full = 0, partial = 0, none = 0
    results.forEach(r => {
      const mainMatch = r.mainType.toLowerCase() === cluster.consensusMainType.toLowerCase()
      const subMatch = r.subType.toLowerCase() === cluster.consensusSubType.toLowerCase()

      if (mainMatch && subMatch) full++
      else if (mainMatch) partial++
      else none++
    })
    return { full, partial, none, total: results.length }
  }, [results, cluster.consensusMainType, cluster.consensusSubType])

  return (
    <div className="bg-white border border-gray-200 rounded-xl shadow-sm overflow-hidden">
      {/* Header */}
      <div className="p-5 flex items-start justify-between gap-4">
        <div className="flex-1">
          <h3 className="text-lg font-bold text-gray-800">{cluster.clusterId}</h3>
          <div className="mt-2 flex items-center gap-4">
            {/* Consensus Score */}
            <div className="flex items-center gap-2">
              <Gauge className="h-5 w-5" style={{ color: scoreInfo.color }} />
              <span className="text-2xl font-bold" style={{ color: scoreInfo.color }}>
                {cluster.consensusScore.toFixed(0)}%
              </span>
              <span className="text-sm text-gray-500">{scoreInfo.label}</span>
            </div>
          </div>
        </div>

        {/* Consensus Types */}
        <div className="text-right">
          <div className="text-sm text-gray-500">Consensus</div>
          <div className="font-semibold text-teal-700">{cluster.consensusMainType}</div>
          <div className="text-sm text-gray-600">{cluster.consensusSubType}</div>
        </div>
      </div>

      {/* Mixed Types */}
      {cluster.possibleMixedTypes && cluster.possibleMixedTypes.length > 0 && (
        <div className="px-5 pb-3">
          <span className="text-xs font-semibold text-gray-500 uppercase">Possible Mixed Types:</span>
          <div className="flex flex-wrap gap-1 mt-1">
            {cluster.possibleMixedTypes.map((type, i) => (
              <span key={i} className="px-2 py-0.5 bg-amber-50 text-amber-700 rounded text-xs">
                {type}
              </span>
            ))}
          </div>
        </div>
      )}

      {/* Agreement Stats */}
      <div className="px-5 pb-4 flex gap-4 text-sm">
        <div className="flex items-center gap-1">
          <CheckCircle className="h-4 w-4 text-green-500" />
          <span className="text-gray-600">{agreements.full}/{agreements.total} full</span>
        </div>
        <div className="flex items-center gap-1">
          <span className="w-3 h-3 rounded-full bg-yellow-400"></span>
          <span className="text-gray-600">{agreements.partial} partial</span>
        </div>
        <div className="flex items-center gap-1">
          <AlertTriangle className="h-4 w-4 text-red-400" />
          <span className="text-gray-600">{agreements.none} disagree</span>
        </div>
      </div>

      {/* Expand/Collapse Button */}
      <button
        onClick={onToggle}
        className="w-full px-5 py-3 bg-gray-50 hover:bg-gray-100 transition-colors flex items-center justify-center gap-2 text-sm text-gray-600 font-medium"
      >
        {expanded ? (
          <>
            <ChevronUp className="h-4 w-4" />
            Hide Details
          </>
        ) : (
          <>
            <ChevronDown className="h-4 w-4" />
            Show Details
          </>
        )}
      </button>

      {/* Expanded Content */}
      {expanded && (
        <div className="p-5 border-t border-gray-100 space-y-6">
          {/* LLM Reasoning */}
          {cluster.llmReasoning && (
            <div>
              <h4 className="text-sm font-semibold text-gray-700 mb-2">LLM Reasoning</h4>
              <div className="bg-teal-50 rounded-lg p-4 text-sm text-gray-700">
                {cluster.llmReasoning}
              </div>
            </div>
          )}

          {/* Pie Charts */}
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            {/* Main Type Distribution */}
            <div>
              <h4 className="text-sm font-semibold text-gray-700 mb-2 flex items-center gap-2">
                <PieChartIcon className="h-4 w-4 text-teal-600" />
                Main Type Distribution
              </h4>
              <div className="h-48">
                <ResponsiveContainer width="100%" height="100%">
                  <PieChart>
                    <Pie
                      data={mainTypeDistribution}
                      cx="50%"
                      cy="50%"
                      innerRadius={40}
                      outerRadius={70}
                      paddingAngle={2}
                      dataKey="value"
                    >
                      {mainTypeDistribution.map((_, index) => (
                        <Cell key={`cell-${index}`} fill={PIE_COLORS[index % PIE_COLORS.length]} />
                      ))}
                    </Pie>
                    <Tooltip />
                    <Legend
                      layout="vertical"
                      align="right"
                      verticalAlign="middle"
                      formatter={(value) => <span className="text-xs">{value}</span>}
                    />
                  </PieChart>
                </ResponsiveContainer>
              </div>
            </div>

            {/* Sub Type Distribution */}
            <div>
              <h4 className="text-sm font-semibold text-gray-700 mb-2 flex items-center gap-2">
                <PieChartIcon className="h-4 w-4 text-teal-600" />
                Sub Type Distribution
              </h4>
              <div className="h-48">
                <ResponsiveContainer width="100%" height="100%">
                  <PieChart>
                    <Pie
                      data={subTypeDistribution}
                      cx="50%"
                      cy="50%"
                      innerRadius={40}
                      outerRadius={70}
                      paddingAngle={2}
                      dataKey="value"
                    >
                      {subTypeDistribution.map((_, index) => (
                        <Cell key={`cell-${index}`} fill={PIE_COLORS[index % PIE_COLORS.length]} />
                      ))}
                    </Pie>
                    <Tooltip />
                    <Legend
                      layout="vertical"
                      align="right"
                      verticalAlign="middle"
                      formatter={(value) => <span className="text-xs">{value}</span>}
                    />
                  </PieChart>
                </ResponsiveContainer>
              </div>
            </div>
          </div>

          {/* Per-Round Results Table */}
          <div>
            <h4 className="text-sm font-semibold text-gray-700 mb-2">Per-Round Results</h4>
            <div className="overflow-x-auto">
              <table className="w-full text-sm">
                <thead>
                  <tr className="bg-teal-600 text-white">
                    <th className="px-3 py-2 text-left font-semibold">Round</th>
                    <th className="px-3 py-2 text-left font-semibold">Main Type</th>
                    <th className="px-3 py-2 text-left font-semibold">Sub Type</th>
                    <th className="px-3 py-2 text-left font-semibold">Agreement</th>
                  </tr>
                </thead>
                <tbody>
                  {results.map((round, i) => {
                    const mainMatch = round.mainType.toLowerCase() === cluster.consensusMainType.toLowerCase()
                    const subMatch = round.subType.toLowerCase() === cluster.consensusSubType.toLowerCase()
                    const agreement = mainMatch && subMatch ? 'full' : mainMatch ? 'partial' : 'none'

                    return (
                      <tr
                        key={i}
                        className={cn(
                          "border-b border-gray-100",
                          i % 2 === 0 ? "bg-white" : "bg-gray-50"
                        )}
                      >
                        <td className="px-3 py-2 font-medium">{round.iteration}</td>
                        <td className={cn(
                          "px-3 py-2",
                          mainMatch ? "text-green-700 font-medium" : "text-gray-700"
                        )}>
                          {round.mainType}
                        </td>
                        <td className={cn(
                          "px-3 py-2",
                          subMatch ? "text-green-700 font-medium" : "text-gray-700"
                        )}>
                          {round.subType}
                        </td>
                        <td className="px-3 py-2">
                          <span className={cn(
                            "px-2 py-0.5 rounded text-xs font-medium",
                            agreement === 'full' && "bg-green-100 text-green-700",
                            agreement === 'partial' && "bg-yellow-100 text-yellow-700",
                            agreement === 'none' && "bg-red-100 text-red-700"
                          )}>
                            {agreement === 'full' ? 'Full Match' :
                             agreement === 'partial' ? 'Partial Match' : 'No Match'}
                          </span>
                        </td>
                      </tr>
                    )
                  })}
                </tbody>
              </table>
            </div>
          </div>
        </div>
      )}
    </div>
  )
}

/**
 * Uncertainty Report component
 * Displays consensus scores and per-round analysis
 */
export function UncertaintyReport({
  data,
  title = 'CASSIA Uncertainty Quantification Report',
  className
}: UncertaintyReportProps) {
  // Track expanded clusters
  const [expandedClusters, setExpandedClusters] = useState<Set<string>>(new Set())

  // Search filter
  const searchFields: (keyof ClusterUncertainty)[] = ['clusterId', 'consensusMainType', 'consensusSubType']
  const {
    searchValue,
    setSearchValue,
    filteredData,
    handleClearAll,
    totalCount,
    visibleCount
  } = useReportFilters(data, searchFields)

  // Timestamp
  const timestamp = useMemo(() => new Date().toLocaleString(), [])

  // Summary stats
  const stats = useMemo(() => {
    const scores = data.map(d => d.consensusScore)
    const avgScore = scores.length > 0 ? scores.reduce((a, b) => a + b, 0) / scores.length : 0
    const highConfidence = scores.filter(s => s >= 80).length
    const lowConfidence = scores.filter(s => s < 50).length

    return {
      avgScore,
      highConfidence,
      lowConfidence,
      totalClusters: data.length
    }
  }, [data])

  // Toggle cluster expansion
  const toggleCluster = (clusterId: string) => {
    setExpandedClusters(prev => {
      const next = new Set(prev)
      if (next.has(clusterId)) {
        next.delete(clusterId)
      } else {
        next.add(clusterId)
      }
      return next
    })
  }

  // Expand/collapse all
  const expandAll = () => {
    setExpandedClusters(new Set(data.map(d => d.clusterId)))
  }

  const collapseAll = () => {
    setExpandedClusters(new Set())
  }

  if (data.length === 0) {
    return (
      <div className={cn("text-center py-12 text-gray-500", className)}>
        <Gauge className="h-12 w-12 mx-auto mb-4 text-gray-300" />
        <p>No uncertainty data available to display.</p>
      </div>
    )
  }

  return (
    <div className={cn("space-y-6", className)}>
      {/* Header */}
      <header className="bg-gradient-to-r from-teal-600 to-cyan-600 text-white rounded-2xl p-8 shadow-lg">
        <div className="flex items-center gap-3 mb-2">
          <Gauge className="h-8 w-8" />
          <h1 className="text-3xl font-bold">{title}</h1>
        </div>
        <p className="text-teal-100 flex items-center gap-2">
          <Calendar className="h-4 w-4" />
          Generated: {timestamp}
        </p>
      </header>

      {/* Summary Stats */}
      <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
        <div className="bg-white rounded-xl border border-gray-200 p-5 text-center shadow-sm">
          <div className="text-3xl font-bold text-teal-600">{stats.totalClusters}</div>
          <div className="text-sm text-gray-500 uppercase tracking-wider mt-1">Total Clusters</div>
        </div>
        <div className="bg-white rounded-xl border border-gray-200 p-5 text-center shadow-sm">
          <div className="text-3xl font-bold text-teal-600">{stats.avgScore.toFixed(1)}%</div>
          <div className="text-sm text-gray-500 uppercase tracking-wider mt-1">Avg Score</div>
        </div>
        <div className="bg-white rounded-xl border border-gray-200 p-5 text-center shadow-sm">
          <div className="text-3xl font-bold text-green-600">{stats.highConfidence}</div>
          <div className="text-sm text-gray-500 uppercase tracking-wider mt-1">High Confidence</div>
        </div>
        <div className="bg-white rounded-xl border border-gray-200 p-5 text-center shadow-sm">
          <div className="text-3xl font-bold text-red-500">{stats.lowConfidence}</div>
          <div className="text-sm text-gray-500 uppercase tracking-wider mt-1">Low Confidence</div>
        </div>
      </div>

      {/* Filters */}
      <div className="flex flex-col md:flex-row gap-4 items-start md:items-center">
        <div className="flex-1 w-full">
          <ReportFilters
            searchPlaceholder="Search clusters, cell types..."
            searchValue={searchValue}
            onSearchChange={setSearchValue}
            onClearAll={handleClearAll}
            totalCount={totalCount}
            visibleCount={visibleCount}
          />
        </div>
        <div className="flex gap-2">
          <Button variant="outline" size="sm" onClick={expandAll}>
            Expand All
          </Button>
          <Button variant="outline" size="sm" onClick={collapseAll}>
            Collapse All
          </Button>
        </div>
      </div>

      {/* Cluster Cards */}
      <div className="space-y-4">
        {filteredData.map((cluster) => (
          <ClusterUncertaintyCard
            key={cluster.clusterId}
            cluster={cluster}
            expanded={expandedClusters.has(cluster.clusterId)}
            onToggle={() => toggleCluster(cluster.clusterId)}
          />
        ))}

        {filteredData.length === 0 && (
          <div className="text-center py-12 text-gray-500">
            <Gauge className="h-12 w-12 mx-auto mb-4 text-gray-300" />
            <p>No clusters match your search.</p>
          </div>
        )}
      </div>
    </div>
  )
}

export default UncertaintyReport
