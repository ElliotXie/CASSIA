'use client'

import React, { useState, useMemo } from 'react'
import { Calendar, Layers, ChevronDown, ChevronUp, FileText } from 'lucide-react'
import { Button } from '@/components/ui/button'
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
} from '@/components/ui/dialog'
import { ScrollArea } from '@/components/ui/scroll-area'
import { ReportFilters, useReportFilters } from './ReportFilters'
import { cn } from '@/lib/utils'

export interface SubclusterResult {
  clusterId: string
  mainCellType: string
  subCellType?: string
  keyMarkers?: string
  reason?: string
  parentCluster?: string
}

export interface SubclusteringReportProps {
  data: SubclusterResult[]
  modelName?: string
  className?: string
}

/**
 * Marker Modal for displaying full marker list
 */
function MarkerModal({
  open,
  onOpenChange,
  clusterId,
  markers
}: {
  open: boolean
  onOpenChange: (open: boolean) => void
  clusterId: string
  markers: string
}) {
  return (
    <Dialog open={open} onOpenChange={onOpenChange}>
      <DialogContent className="max-w-2xl">
        <DialogHeader>
          <DialogTitle className="text-lg font-bold text-blue-800">
            Markers for {clusterId}
          </DialogTitle>
        </DialogHeader>
        <ScrollArea className="max-h-96">
          <div className="p-4 bg-blue-50 rounded-lg">
            <p className="text-sm text-gray-700 leading-relaxed whitespace-pre-wrap">
              {markers || 'No markers available'}
            </p>
          </div>
        </ScrollArea>
      </DialogContent>
    </Dialog>
  )
}

/**
 * Subclustering Report component
 * Displays hierarchical annotation results with marker modals
 */
export function SubclusteringReport({
  data,
  modelName,
  className
}: SubclusteringReportProps) {
  // Modal state
  const [selectedCluster, setSelectedCluster] = useState<{ id: string; markers: string } | null>(null)
  const [markerModalOpen, setMarkerModalOpen] = useState(false)

  // Search filter
  const searchFields: (keyof SubclusterResult)[] = ['clusterId', 'mainCellType', 'subCellType']
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

  // Group by parent cluster if available
  const groupedData = useMemo(() => {
    const groups: Record<string, SubclusterResult[]> = {}

    filteredData.forEach(item => {
      const parent = item.parentCluster || 'Main Clusters'
      if (!groups[parent]) groups[parent] = []
      groups[parent].push(item)
    })

    return groups
  }, [filteredData])

  // Handle marker view
  const handleShowMarkers = (clusterId: string, markers: string) => {
    setSelectedCluster({ id: clusterId, markers })
    setMarkerModalOpen(true)
  }

  if (data.length === 0) {
    return (
      <div className={cn("text-center py-12 text-gray-500", className)}>
        <Layers className="h-12 w-12 mx-auto mb-4 text-gray-300" />
        <p>No subclustering data available to display.</p>
      </div>
    )
  }

  return (
    <div className={cn("space-y-6", className)}>
      {/* Header */}
      <header className="bg-gradient-to-r from-blue-600 to-cyan-600 text-white rounded-2xl p-8 shadow-lg">
        <div className="flex items-center gap-3 mb-2">
          <Layers className="h-8 w-8" />
          <h1 className="text-3xl font-bold">
            Subclustering Annotation Report
            {modelName && <span className="text-blue-200"> - {modelName}</span>}
          </h1>
        </div>
        <p className="text-blue-100 flex items-center gap-2">
          <Calendar className="h-4 w-4" />
          Generated: {timestamp}
        </p>
      </header>

      {/* Filters */}
      <ReportFilters
        searchPlaceholder="Search clusters, cell types..."
        searchValue={searchValue}
        onSearchChange={setSearchValue}
        onClearAll={handleClearAll}
        totalCount={totalCount}
        visibleCount={visibleCount}
      />

      {/* Results Table */}
      <div className="bg-white rounded-xl border border-gray-200 shadow-sm overflow-hidden">
        <div className="overflow-x-auto">
          <table className="w-full text-sm">
            <thead>
              <tr className="bg-blue-600 text-white">
                <th className="px-4 py-3 text-left font-semibold w-24">Cluster</th>
                <th className="px-4 py-3 text-left font-semibold w-1/4">Annotation</th>
                <th className="px-4 py-3 text-left font-semibold w-32">Markers</th>
                <th className="px-4 py-3 text-left font-semibold">Reasoning</th>
              </tr>
            </thead>
            <tbody>
              {Object.entries(groupedData).map(([parentCluster, clusters]) => (
                <React.Fragment key={parentCluster}>
                  {/* Parent cluster header if there are multiple groups */}
                  {Object.keys(groupedData).length > 1 && (
                    <tr className="bg-blue-50">
                      <td colSpan={4} className="px-4 py-2">
                        <span className="font-semibold text-blue-700 flex items-center gap-2">
                          <Layers className="h-4 w-4" />
                          {parentCluster}
                        </span>
                      </td>
                    </tr>
                  )}

                  {clusters.map((row, i) => (
                    <tr
                      key={row.clusterId || i}
                      className={cn(
                        "border-b border-gray-100 hover:bg-blue-50 transition-colors",
                        i % 2 === 0 ? "bg-white" : "bg-gray-50"
                      )}
                    >
                      <td className="px-4 py-3">
                        <span className="font-bold text-blue-600 text-lg">
                          {row.clusterId}
                        </span>
                      </td>
                      <td className="px-4 py-3">
                        <div className="font-semibold text-gray-800">
                          {row.mainCellType}
                        </div>
                        {row.subCellType && (
                          <div className="text-gray-500 text-sm mt-0.5">
                            {row.subCellType}
                          </div>
                        )}
                      </td>
                      <td className="px-4 py-3">
                        {row.keyMarkers ? (
                          <Button
                            variant="outline"
                            size="sm"
                            onClick={() => handleShowMarkers(row.clusterId, row.keyMarkers || '')}
                            className="text-blue-600 border-blue-200 hover:bg-blue-50"
                          >
                            Show Markers
                          </Button>
                        ) : (
                          <span className="text-gray-400 text-sm">N/A</span>
                        )}
                      </td>
                      <td className="px-4 py-3">
                        <div className="bg-gray-50 rounded-lg p-3 text-gray-600 text-sm max-w-lg">
                          {row.reason || 'No reasoning available'}
                        </div>
                      </td>
                    </tr>
                  ))}
                </React.Fragment>
              ))}

              {filteredData.length === 0 && (
                <tr>
                  <td colSpan={4} className="px-4 py-12 text-center text-gray-500">
                    No clusters match your search.
                  </td>
                </tr>
              )}
            </tbody>
          </table>
        </div>
      </div>

      {/* Summary Stats */}
      <div className="bg-white rounded-xl border border-gray-200 p-6 shadow-sm">
        <h2 className="text-lg font-bold text-gray-800 mb-4 flex items-center gap-2">
          <FileText className="h-5 w-5 text-blue-600" />
          Summary
        </h2>
        <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
          <div className="bg-blue-50 rounded-lg p-4 text-center">
            <div className="text-2xl font-bold text-blue-700">{data.length}</div>
            <div className="text-sm text-gray-600">Total Subclusters</div>
          </div>
          <div className="bg-blue-50 rounded-lg p-4 text-center">
            <div className="text-2xl font-bold text-blue-700">
              {new Set(data.map(d => d.mainCellType)).size}
            </div>
            <div className="text-sm text-gray-600">Unique Main Types</div>
          </div>
          <div className="bg-blue-50 rounded-lg p-4 text-center">
            <div className="text-2xl font-bold text-blue-700">
              {new Set(data.map(d => d.subCellType).filter(Boolean)).size}
            </div>
            <div className="text-sm text-gray-600">Unique Sub Types</div>
          </div>
          <div className="bg-blue-50 rounded-lg p-4 text-center">
            <div className="text-2xl font-bold text-blue-700">
              {new Set(data.map(d => d.parentCluster).filter(Boolean)).size || 1}
            </div>
            <div className="text-sm text-gray-600">Parent Clusters</div>
          </div>
        </div>
      </div>

      {/* Marker Modal */}
      <MarkerModal
        open={markerModalOpen}
        onOpenChange={setMarkerModalOpen}
        clusterId={selectedCluster?.id || ''}
        markers={selectedCluster?.markers || ''}
      />
    </div>
  )
}

export default SubclusteringReport
