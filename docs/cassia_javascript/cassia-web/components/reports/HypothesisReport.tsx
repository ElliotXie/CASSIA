'use client'

import React, { useState, useMemo } from 'react'
import { Calendar, ChevronDown, ChevronUp, FileText, Lightbulb } from 'lucide-react'
import { Button } from '@/components/ui/button'
import {
  Accordion,
  AccordionContent,
  AccordionItem,
  AccordionTrigger,
} from '@/components/ui/accordion'
import { ReportFilters, useReportFilters } from './ReportFilters'
import { cn } from '@/lib/utils'

export interface Hypothesis {
  rank: number | string
  cellType: string
  reasoning: string
}

export interface ClusterHypothesis {
  clusterName: string
  hypotheses: Hypothesis[]
  rawResponse?: string
  error?: string
}

export interface HypothesisReportProps {
  data: ClusterHypothesis[]
  title?: string
  sourceFile?: string
  className?: string
}

/**
 * Hypothesis Report component
 * Displays ranked cell type hypotheses for each cluster
 */
export function HypothesisReport({
  data,
  title = 'CASSIA Hypothesis Generation Report',
  sourceFile,
  className
}: HypothesisReportProps) {
  // Track which raw outputs are expanded
  const [expandedRaw, setExpandedRaw] = useState<Set<string>>(new Set())

  // Search filter
  const searchFields: (keyof ClusterHypothesis)[] = ['clusterName']
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

  // Toggle raw output visibility
  const toggleRaw = (clusterName: string) => {
    setExpandedRaw(prev => {
      const next = new Set(prev)
      if (next.has(clusterName)) {
        next.delete(clusterName)
      } else {
        next.add(clusterName)
      }
      return next
    })
  }

  if (data.length === 0) {
    return (
      <div className={cn("text-center py-12 text-gray-500", className)}>
        <Lightbulb className="h-12 w-12 mx-auto mb-4 text-gray-300" />
        <p>No hypothesis data available to display.</p>
      </div>
    )
  }

  return (
    <div className={cn("space-y-6", className)}>
      {/* Header */}
      <header className="bg-gradient-to-r from-blue-600 to-indigo-600 text-white rounded-2xl p-8 shadow-lg">
        <div className="text-center">
          <h1 className="text-3xl font-bold mb-2">{title}</h1>
          <p className="text-blue-100 flex items-center justify-center gap-2">
            <Calendar className="h-4 w-4" />
            Generated on: {timestamp}
          </p>
          {sourceFile && (
            <p className="text-blue-100 flex items-center justify-center gap-2 mt-1">
              <FileText className="h-4 w-4" />
              Source: {sourceFile}
            </p>
          )}
        </div>
      </header>

      {/* Filters */}
      <ReportFilters
        searchPlaceholder="Search clusters..."
        searchValue={searchValue}
        onSearchChange={setSearchValue}
        onClearAll={handleClearAll}
        totalCount={totalCount}
        visibleCount={visibleCount}
      />

      {/* Cluster Sections */}
      <div className="space-y-6">
        {filteredData.map((cluster, index) => (
          <div
            key={cluster.clusterName || index}
            className="bg-white border border-gray-200 rounded-xl shadow-sm overflow-hidden"
          >
            {/* Cluster Header */}
            <div className="bg-gradient-to-r from-blue-50 to-indigo-50 px-6 py-4 border-b border-gray-200">
              <h2 className="text-xl font-bold text-gray-800 flex items-center gap-2">
                <Lightbulb className="h-5 w-5 text-blue-500" />
                Cluster: {cluster.clusterName}
              </h2>
            </div>

            {/* Hypotheses Table */}
            <div className="p-6">
              <h3 className="font-semibold text-gray-700 mb-4">Top Hypotheses</h3>

              {cluster.error ? (
                <div className="bg-red-50 border border-red-200 rounded-lg p-4 text-red-700">
                  Error: {cluster.error}
                </div>
              ) : cluster.hypotheses.length > 0 ? (
                <div className="overflow-x-auto">
                  <table className="w-full text-sm">
                    <thead>
                      <tr className="bg-blue-600 text-white">
                        <th className="px-4 py-3 text-left font-semibold w-20">Rank</th>
                        <th className="px-4 py-3 text-left font-semibold w-1/4">Predicted Cell Type</th>
                        <th className="px-4 py-3 text-left font-semibold">Reasoning</th>
                      </tr>
                    </thead>
                    <tbody>
                      {cluster.hypotheses.map((hyp, i) => (
                        <tr
                          key={i}
                          className={cn(
                            "border-b border-gray-100 hover:bg-blue-50 transition-colors",
                            i % 2 === 0 ? "bg-white" : "bg-gray-50"
                          )}
                        >
                          <td className="px-4 py-3">
                            <span className={cn(
                              "inline-flex items-center justify-center w-8 h-8 rounded-full font-bold text-sm",
                              i === 0 && "bg-gradient-to-br from-yellow-400 to-amber-500 text-white",
                              i === 1 && "bg-gradient-to-br from-gray-300 to-gray-400 text-gray-800",
                              i === 2 && "bg-gradient-to-br from-amber-600 to-amber-700 text-white",
                              i > 2 && "bg-gray-200 text-gray-600"
                            )}>
                              {hyp.rank}
                            </span>
                          </td>
                          <td className="px-4 py-3 font-medium text-gray-800">
                            {hyp.cellType}
                          </td>
                          <td className="px-4 py-3 text-gray-600">
                            {hyp.reasoning}
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              ) : (
                <p className="text-gray-500 italic">No hypotheses generated for this cluster.</p>
              )}

              {/* Raw LLM Output (Collapsible) */}
              {cluster.rawResponse && (
                <div className="mt-6">
                  <button
                    onClick={() => toggleRaw(cluster.clusterName)}
                    className="flex items-center gap-2 text-blue-600 hover:text-blue-700 font-semibold text-sm"
                  >
                    {expandedRaw.has(cluster.clusterName) ? (
                      <ChevronUp className="h-4 w-4" />
                    ) : (
                      <ChevronDown className="h-4 w-4" />
                    )}
                    View Raw LLM Output
                  </button>

                  {expandedRaw.has(cluster.clusterName) && (
                    <pre className="mt-3 bg-gray-800 text-gray-100 p-4 rounded-lg text-sm overflow-x-auto whitespace-pre-wrap">
                      {cluster.rawResponse}
                    </pre>
                  )}
                </div>
              )}
            </div>
          </div>
        ))}

        {filteredData.length === 0 && (
          <div className="text-center py-12 text-gray-500">
            <Lightbulb className="h-12 w-12 mx-auto mb-4 text-gray-300" />
            <p>No clusters match your search.</p>
          </div>
        )}
      </div>

      {/* Footer */}
      <footer className="text-center text-gray-500 text-sm py-4">
        Report generated by CASSIA
      </footer>
    </div>
  )
}

export default HypothesisReport
