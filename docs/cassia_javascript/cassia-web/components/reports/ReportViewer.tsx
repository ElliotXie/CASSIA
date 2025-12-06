'use client'

import React, { useState } from 'react'
import { FileText, X, Maximize2, Minimize2, Download } from 'lucide-react'
import { Button } from '@/components/ui/button'
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
} from '@/components/ui/dialog'
import { Tabs, TabsContent, TabsList, TabsTrigger } from '@/components/ui/tabs'
import { BatchReport, type BatchReportData } from './BatchReport'
import { HypothesisReport, type ClusterHypothesis } from './HypothesisReport'
import { EvaluationReport, type EvaluationResult } from './EvaluationReport'
import { SubclusteringReport, type SubclusterResult } from './SubclusteringReport'
import { UncertaintyReport, type ClusterUncertainty } from './UncertaintyReport'
import { cn } from '@/lib/utils'

export type ReportType = 'batch' | 'hypothesis' | 'evaluation' | 'subclustering' | 'uncertainty'

export interface ReportData {
  batch?: BatchReportData[]
  hypothesis?: ClusterHypothesis[]
  evaluation?: EvaluationResult[]
  subclustering?: SubclusterResult[]
  uncertainty?: ClusterUncertainty[]
}

export interface ReportViewerProps {
  data: ReportData
  availableReports?: ReportType[]
  defaultReport?: ReportType
  title?: string
  className?: string
  onExportHTML?: (reportType: ReportType) => void
}

/**
 * Unified Report Viewer component
 * Displays multiple report types with tab navigation
 */
export function ReportViewer({
  data,
  availableReports,
  defaultReport,
  title = 'Analysis Reports',
  className,
  onExportHTML
}: ReportViewerProps) {
  // Determine available reports based on data
  const reports = availableReports || Object.entries(data)
    .filter(([_, value]) => value && (Array.isArray(value) ? value.length > 0 : true))
    .map(([key]) => key as ReportType)

  const [activeReport, setActiveReport] = useState<ReportType>(
    defaultReport || reports[0] || 'batch'
  )

  const reportLabels: Record<ReportType, string> = {
    batch: 'Batch Analysis',
    hypothesis: 'Hypotheses',
    evaluation: 'Evaluation',
    subclustering: 'Subclustering',
    uncertainty: 'Uncertainty'
  }

  if (reports.length === 0) {
    return (
      <div className={cn("text-center py-12 text-gray-500", className)}>
        <FileText className="h-12 w-12 mx-auto mb-4 text-gray-300" />
        <p>No report data available.</p>
      </div>
    )
  }

  return (
    <div className={cn("space-y-4", className)}>
      {/* Header */}
      <div className="flex items-center justify-between">
        <h2 className="text-xl font-bold text-gray-800 flex items-center gap-2">
          <FileText className="h-5 w-5 text-teal-600" />
          {title}
        </h2>
        {onExportHTML && (
          <Button
            variant="outline"
            size="sm"
            onClick={() => onExportHTML(activeReport)}
            className="border-teal-200 hover:border-teal-400"
          >
            <Download className="h-4 w-4 mr-2" />
            Export {reportLabels[activeReport]}
          </Button>
        )}
      </div>

      {/* Report Tabs */}
      {reports.length > 1 ? (
        <Tabs value={activeReport} onValueChange={(v) => setActiveReport(v as ReportType)}>
          <TabsList className="bg-teal-50 border border-teal-100">
            {reports.map((report) => (
              <TabsTrigger
                key={report}
                value={report}
                className="data-[state=active]:bg-teal-600 data-[state=active]:text-white"
              >
                {reportLabels[report]}
              </TabsTrigger>
            ))}
          </TabsList>

          {data.batch && (
            <TabsContent value="batch" className="mt-4">
              <BatchReport data={data.batch} />
            </TabsContent>
          )}

          {data.hypothesis && (
            <TabsContent value="hypothesis" className="mt-4">
              <HypothesisReport data={data.hypothesis} />
            </TabsContent>
          )}

          {data.evaluation && (
            <TabsContent value="evaluation" className="mt-4">
              <EvaluationReport data={data.evaluation} />
            </TabsContent>
          )}

          {data.subclustering && (
            <TabsContent value="subclustering" className="mt-4">
              <SubclusteringReport data={data.subclustering} />
            </TabsContent>
          )}

          {data.uncertainty && (
            <TabsContent value="uncertainty" className="mt-4">
              <UncertaintyReport data={data.uncertainty} />
            </TabsContent>
          )}
        </Tabs>
      ) : (
        // Single report - no tabs needed
        <div className="mt-4">
          {activeReport === 'batch' && data.batch && (
            <BatchReport data={data.batch} />
          )}
          {activeReport === 'hypothesis' && data.hypothesis && (
            <HypothesisReport data={data.hypothesis} />
          )}
          {activeReport === 'evaluation' && data.evaluation && (
            <EvaluationReport data={data.evaluation} />
          )}
          {activeReport === 'subclustering' && data.subclustering && (
            <SubclusteringReport data={data.subclustering} />
          )}
          {activeReport === 'uncertainty' && data.uncertainty && (
            <UncertaintyReport data={data.uncertainty} />
          )}
        </div>
      )}
    </div>
  )
}

/**
 * Modal version of Report Viewer for full-screen display
 */
export function ReportViewerModal({
  open,
  onOpenChange,
  data,
  availableReports,
  defaultReport,
  title,
  onExportHTML
}: ReportViewerProps & {
  open: boolean
  onOpenChange: (open: boolean) => void
}) {
  const [isMaximized, setIsMaximized] = useState(true)

  return (
    <Dialog open={open} onOpenChange={onOpenChange}>
      <DialogContent
        className={cn(
          "p-0 overflow-hidden",
          isMaximized
            ? "max-w-[95vw] w-[95vw] max-h-[95vh] h-[95vh]"
            : "max-w-4xl max-h-[80vh]"
        )}
      >
        <DialogHeader className="px-6 py-4 border-b bg-gradient-to-r from-teal-50 to-emerald-50 flex flex-row items-center justify-between">
          <DialogTitle className="flex items-center gap-2 text-lg font-bold text-gray-800">
            <FileText className="h-5 w-5 text-teal-600" />
            {title || 'Analysis Reports'}
          </DialogTitle>
          <div className="flex items-center gap-2">
            <Button
              variant="ghost"
              size="sm"
              onClick={() => setIsMaximized(!isMaximized)}
              className="h-8 w-8 p-0"
            >
              {isMaximized ? (
                <Minimize2 className="h-4 w-4" />
              ) : (
                <Maximize2 className="h-4 w-4" />
              )}
            </Button>
          </div>
        </DialogHeader>

        <div className="overflow-auto p-6" style={{ maxHeight: isMaximized ? 'calc(95vh - 80px)' : 'calc(80vh - 80px)' }}>
          <ReportViewer
            data={data}
            availableReports={availableReports}
            defaultReport={defaultReport}
            title=""
            onExportHTML={onExportHTML}
          />
        </div>
      </DialogContent>
    </Dialog>
  )
}

export default ReportViewer
