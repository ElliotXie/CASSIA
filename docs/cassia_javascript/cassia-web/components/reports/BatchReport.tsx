'use client'

import React, { useState, useMemo } from 'react'
import { Calendar, Beaker, Globe, Download, Printer } from 'lucide-react'
import { Button } from '@/components/ui/button'
import { ReportFilters, useReportFilters, FilterConfig } from './ReportFilters'
import { ClusterCard } from './ClusterCard'
import { ConversationModal, ClusterData } from './ConversationModal'
import { cn } from '@/lib/utils'

export interface BatchReportData {
  clusterId: string
  mainType: string
  subTypes?: string | string[]
  mixedTypes?: string | string[]
  markers?: string
  markerCount?: number
  score?: number | string
  tissue?: string
  species?: string
  model?: string
  provider?: string
  iterations?: number
  conversationHistory?: string
  mergedGrouping1?: string
  mergedGrouping2?: string
  mergedGrouping3?: string
}

export interface BatchReportProps {
  data: BatchReportData[]
  title?: string
  subtitle?: string
  className?: string
  onExportHTML?: () => void
  onExportPDF?: () => void
}

/**
 * Transform raw data to ClusterData format
 */
function transformToClusterData(item: BatchReportData): ClusterData {
  // Handle subTypes - can be string (comma-separated) or array
  let subTypes: string[] = []
  if (Array.isArray(item.subTypes)) {
    subTypes = item.subTypes
  } else if (typeof item.subTypes === 'string' && item.subTypes) {
    subTypes = item.subTypes.split(',').map(s => s.trim()).filter(Boolean)
  }

  // Handle mixedTypes - can be string (comma-separated) or array
  let mixedTypes: string[] = []
  if (Array.isArray(item.mixedTypes)) {
    mixedTypes = item.mixedTypes
  } else if (typeof item.mixedTypes === 'string' && item.mixedTypes) {
    mixedTypes = item.mixedTypes.split(',').map(s => s.trim()).filter(Boolean)
  }

  return {
    clusterId: item.clusterId,
    mainType: item.mainType,
    subTypes,
    mixedTypes,
    markers: item.markers,
    markerCount: item.markerCount,
    score: item.score,
    tissue: item.tissue,
    species: item.species,
    model: item.model,
    provider: item.provider,
    iterations: item.iterations,
    conversationHistory: item.conversationHistory,
    mergedGrouping1: item.mergedGrouping1,
    mergedGrouping2: item.mergedGrouping2,
    mergedGrouping3: item.mergedGrouping3
  }
}

/**
 * Main Batch Analysis Report component
 * Displays cluster cards with filtering and conversation modals
 */
export function BatchReport({
  data,
  title = 'CASSIA Batch Analysis Report',
  subtitle = 'Comprehensive Cell Type Annotation Analysis',
  className,
  onExportHTML,
  onExportPDF
}: BatchReportProps) {
  // Modal state
  const [selectedCluster, setSelectedCluster] = useState<ClusterData | null>(null)
  const [modalOpen, setModalOpen] = useState(false)

  // Transform data
  const transformedData = useMemo(() => data.map(transformToClusterData), [data])

  // Filter configuration
  const filterFields: (keyof BatchReportData)[] = ['tissue', 'species', 'model', 'provider']
  const searchFields: (keyof BatchReportData)[] = ['clusterId', 'mainType', 'markers']

  const {
    searchValue,
    setSearchValue,
    filterValues,
    filterOptions,
    filteredData,
    handleFilterChange,
    handleClearAll,
    totalCount,
    visibleCount
  } = useReportFilters(data, searchFields, filterFields)

  // Build filter configs
  const filters: FilterConfig[] = useMemo(() => {
    return filterFields.map(field => ({
      key: String(field),
      label: String(field).charAt(0).toUpperCase() + String(field).slice(1),
      options: filterOptions[String(field)] || [],
      placeholder: `All ${String(field).charAt(0).toUpperCase() + String(field).slice(1)}s`
    }))
  }, [filterOptions])

  // Get filtered transformed data
  const displayData = useMemo(() => {
    const filteredIds = new Set(filteredData.map(d => d.clusterId))
    return transformedData.filter(d => filteredIds.has(d.clusterId))
  }, [transformedData, filteredData])

  // Statistics
  const stats = useMemo(() => {
    const tissues = new Set(data.map(d => d.tissue).filter(Boolean))
    const models = new Set(data.map(d => d.model).filter(Boolean))
    return {
      totalClusters: data.length,
      uniqueTissues: tissues.size,
      uniqueModels: models.size
    }
  }, [data])

  // Generate timestamp
  const timestamp = useMemo(() => {
    return new Date().toLocaleString()
  }, [])

  // Handle card click
  const handleCardClick = (cluster: ClusterData) => {
    setSelectedCluster(cluster)
    setModalOpen(true)
  }

  // Handle print/export
  const handlePrint = () => {
    if (onExportPDF) {
      onExportPDF()
    } else {
      window.print()
    }
  }

  if (data.length === 0) {
    return (
      <div className={cn("text-center py-12 text-gray-500", className)}>
        <p>No batch data available to display.</p>
      </div>
    )
  }

  return (
    <div className={cn("space-y-6", className)}>
      {/* Header */}
      <header className="bg-gradient-to-r from-teal-600 to-emerald-600 text-white rounded-2xl p-8 shadow-lg relative overflow-hidden">
        <div className="absolute top-0 right-0 w-1/2 h-full bg-gradient-radial from-white/10 to-transparent pointer-events-none" />
        <div className="relative z-10">
          <h1 className="text-3xl font-bold mb-2">{title}</h1>
          <p className="text-teal-100 text-lg">{subtitle}</p>
          <div className="flex flex-wrap gap-4 mt-5">
            <span className="flex items-center gap-2 bg-white/20 backdrop-blur-sm px-4 py-2 rounded-lg text-sm">
              <Calendar className="h-4 w-4" />
              {timestamp}
            </span>
            {data[0]?.tissue && (
              <span className="flex items-center gap-2 bg-white/20 backdrop-blur-sm px-4 py-2 rounded-lg text-sm">
                <Beaker className="h-4 w-4" />
                {data[0].tissue}
              </span>
            )}
            {data[0]?.species && (
              <span className="flex items-center gap-2 bg-white/20 backdrop-blur-sm px-4 py-2 rounded-lg text-sm">
                <Globe className="h-4 w-4" />
                {data[0].species}
              </span>
            )}
          </div>
        </div>
      </header>

      {/* Stats Bar */}
      <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
        <div className="bg-white/90 backdrop-blur-sm border border-teal-100 rounded-xl p-5 text-center shadow-sm">
          <div className="text-3xl font-bold text-teal-600">{stats.totalClusters}</div>
          <div className="text-sm text-gray-500 uppercase tracking-wider mt-1">Total Clusters</div>
        </div>
        <div className="bg-white/90 backdrop-blur-sm border border-teal-100 rounded-xl p-5 text-center shadow-sm">
          <div className="text-3xl font-bold text-teal-600">{visibleCount}</div>
          <div className="text-sm text-gray-500 uppercase tracking-wider mt-1">Showing</div>
        </div>
        <div className="bg-white/90 backdrop-blur-sm border border-teal-100 rounded-xl p-5 text-center shadow-sm">
          <div className="text-3xl font-bold text-teal-600">{stats.uniqueTissues}</div>
          <div className="text-sm text-gray-500 uppercase tracking-wider mt-1">Tissues</div>
        </div>
        <div className="bg-white/90 backdrop-blur-sm border border-teal-100 rounded-xl p-5 text-center shadow-sm">
          <div className="text-3xl font-bold text-teal-600">{stats.uniqueModels}</div>
          <div className="text-sm text-gray-500 uppercase tracking-wider mt-1">Models</div>
        </div>
      </div>

      {/* Filters */}
      <div className="flex flex-col md:flex-row gap-4 items-start md:items-center">
        <div className="flex-1 w-full">
          <ReportFilters
            searchValue={searchValue}
            onSearchChange={setSearchValue}
            filters={filters}
            filterValues={filterValues}
            onFilterChange={handleFilterChange}
            onClearAll={handleClearAll}
            totalCount={totalCount}
            visibleCount={visibleCount}
          />
        </div>
        <div className="flex gap-2 print:hidden">
          {onExportHTML && (
            <Button variant="outline" onClick={onExportHTML} className="border-teal-200 hover:border-teal-400">
              <Download className="h-4 w-4 mr-2" />
              Export HTML
            </Button>
          )}
          <Button onClick={handlePrint} className="bg-gradient-to-r from-teal-600 to-emerald-600 text-white">
            <Printer className="h-4 w-4 mr-2" />
            Export PDF
          </Button>
        </div>
      </div>

      {/* Cluster Grid */}
      {displayData.length > 0 ? (
        <div className="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-3 gap-6">
          {displayData.map((cluster, index) => (
            <ClusterCard
              key={cluster.clusterId || index}
              data={cluster}
              onClick={() => handleCardClick(cluster)}
            />
          ))}
        </div>
      ) : (
        <div className="text-center py-16 text-gray-500">
          <div className="w-16 h-16 mx-auto mb-4 rounded-full bg-gray-100 flex items-center justify-center">
            <svg className="w-8 h-8 text-gray-400" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
              <circle cx="11" cy="11" r="8" />
              <path d="M21 21l-4.35-4.35" />
            </svg>
          </div>
          <h3 className="text-lg font-semibold text-gray-700 mb-2">No clusters found</h3>
          <p className="text-gray-400">Try adjusting your search or filters</p>
        </div>
      )}

      {/* Conversation Modal */}
      <ConversationModal
        open={modalOpen}
        onOpenChange={setModalOpen}
        data={selectedCluster}
      />
    </div>
  )
}

export default BatchReport
