// CASSIA Report Components
// Comprehensive report generation for all CASSIA analysis types

// Shared Components
export { ScoreBadge, type ScoreBadgeProps } from './ScoreBadge'
export { ReportFilters, useReportFilters, type ReportFiltersProps, type FilterConfig, type FilterOption } from './ReportFilters'
export { ConversationModal, type ConversationModalProps, type ClusterData, type ConversationSection } from './ConversationModal'
export { ClusterCard, type ClusterCardProps } from './ClusterCard'

// Report Components
export { BatchReport, type BatchReportProps, type BatchReportData } from './BatchReport'
export { HypothesisReport, type HypothesisReportProps, type ClusterHypothesis, type Hypothesis } from './HypothesisReport'
export { EvaluationReport, type EvaluationReportProps, type EvaluationResult, type EvaluationMetrics } from './EvaluationReport'
export { SubclusteringReport, type SubclusteringReportProps, type SubclusterResult } from './SubclusteringReport'
export { UncertaintyReport, type UncertaintyReportProps, type ClusterUncertainty, type RoundResult } from './UncertaintyReport'

// Unified Report Viewer
export { ReportViewer, ReportViewerModal, type ReportViewerProps, type ReportData, type ReportType } from './ReportViewer'

// HTML Export Utilities
export {
  exportBatchReportHTML,
  exportHypothesisReportHTML,
  exportEvaluationReportHTML,
  exportSubclusteringReportHTML,
  exportUncertaintyReportHTML,
  downloadHTML
} from './exportHTML'
