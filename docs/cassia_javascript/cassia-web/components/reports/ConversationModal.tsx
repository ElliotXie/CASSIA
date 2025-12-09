'use client'

import React from 'react'
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
} from '@/components/ui/dialog'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Search, CheckCircle, FileText, Star, Users } from 'lucide-react'
import { ScoreBadge } from './ScoreBadge'
import { cn } from '@/lib/utils'

export interface ConversationSection {
  annotations: string[]    // List of all annotation attempts
  validators: string[]     // List of all validation attempts
  formatting: string
  scoring: string
}

export interface ClusterData {
  clusterId: string
  mainType: string
  subTypes?: string[]
  mixedTypes?: string[]
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

export interface ConversationModalProps {
  open: boolean
  onOpenChange: (open: boolean) => void
  data: ClusterData | null
}

/**
 * Parse conversation history into structured sections
 * Collects all annotation and validation attempts (for multi-iteration scenarios)
 */
function parseConversationHistory(history: string): ConversationSection {
  const result: ConversationSection = {
    annotations: [],
    validators: [],
    formatting: '',
    scoring: ''
  }

  if (!history) return result

  const sections = history.split(' | ')

  for (const section of sections) {
    const trimmed = section.trim()
    if (trimmed.startsWith('Final Annotation Agent:')) {
      result.annotations.push(trimmed.replace('Final Annotation Agent:', '').trim())
    } else if (trimmed.startsWith('Coupling Validator:')) {
      result.validators.push(trimmed.replace('Coupling Validator:', '').trim())
    } else if (trimmed.startsWith('Formatting Agent:')) {
      result.formatting = trimmed.replace('Formatting Agent:', '').trim()
    } else if (trimmed.startsWith('Scoring Agent:')) {
      result.scoring = trimmed.replace('Scoring Agent:', '').trim()
    }
  }

  return result
}

/**
 * Format text with basic markdown-like rendering
 */
function formatAnalysisText(text: string): React.ReactNode {
  if (!text) return <p className="text-gray-500 italic">No content available.</p>

  // Split by lines and process
  const lines = text.split('\n')
  const elements: React.ReactNode[] = []

  lines.forEach((line, index) => {
    const trimmed = line.trim()
    if (!trimmed) {
      elements.push(<br key={`br-${index}`} />)
      return
    }

    // Headers
    if (trimmed.startsWith('### ')) {
      elements.push(
        <h4 key={index} className="text-base font-semibold text-gray-800 mt-4 mb-2">
          {trimmed.slice(4)}
        </h4>
      )
      return
    }
    if (trimmed.startsWith('## ')) {
      elements.push(
        <h3 key={index} className="text-lg font-semibold text-gray-800 mt-4 mb-2">
          {trimmed.slice(3)}
        </h3>
      )
      return
    }
    if (trimmed.startsWith('# ')) {
      elements.push(
        <h2 key={index} className="text-xl font-bold text-gray-800 mt-4 mb-2">
          {trimmed.slice(2)}
        </h2>
      )
      return
    }

    // Step headers
    if (/^Step \d+[\.\:]/.test(trimmed)) {
      elements.push(
        <h4 key={index} className="text-base font-semibold text-teal-700 mt-4 mb-2 pb-1 border-b-2 border-teal-100">
          {trimmed}
        </h4>
      )
      return
    }

    // Bullet points
    if (trimmed.startsWith('* ') || trimmed.startsWith('- ')) {
      elements.push(
        <li key={index} className="ml-4 mb-1 text-gray-700">
          {formatInlineMarkdown(trimmed.slice(2))}
        </li>
      )
      return
    }

    // Numbered lists
    if (/^\d+\. /.test(trimmed)) {
      const content = trimmed.replace(/^\d+\. /, '')
      elements.push(
        <li key={index} className="ml-4 mb-1 text-gray-700 list-decimal">
          {formatInlineMarkdown(content)}
        </li>
      )
      return
    }

    // Regular paragraph
    elements.push(
      <p key={index} className="text-gray-700 mb-2">
        {formatInlineMarkdown(trimmed)}
      </p>
    )
  })

  return <div className="space-y-1">{elements}</div>
}

/**
 * Format inline markdown (bold, italic)
 */
function formatInlineMarkdown(text: string): React.ReactNode {
  // Simple bold/italic parsing
  const parts = text.split(/(\*\*[^*]+\*\*|\*[^*]+\*)/g)

  return parts.map((part, i) => {
    if (part.startsWith('**') && part.endsWith('**')) {
      return <strong key={i} className="font-semibold">{part.slice(2, -2)}</strong>
    }
    if (part.startsWith('*') && part.endsWith('*')) {
      return <em key={i}>{part.slice(1, -1)}</em>
    }
    return part
  })
}

/**
 * Parse and display JSON summary from formatting agent
 */
function FormattingSummary({ text }: { text: string }) {
  try {
    const jsonMatch = text.match(/\{[\s\S]*\}/)
    if (jsonMatch) {
      const data = JSON.parse(jsonMatch[0])
      const mainType = data.main_cell_type || 'Not specified'
      const subTypes = data.sub_cell_types || []
      const mixedTypes = data.possible_mixed_cell_types || []

      return (
        <div className="space-y-4">
          <div>
            <span className="text-xs font-semibold uppercase tracking-wider text-amber-600">
              Main Cell Type
            </span>
            <p className="text-lg font-bold text-teal-800 mt-1">{mainType}</p>
          </div>

          <div>
            <span className="text-xs font-semibold uppercase tracking-wider text-amber-600">
              Sub Cell Types
            </span>
            {subTypes.length > 0 ? (
              <ul className="mt-2 space-y-1.5">
                {subTypes.map((type: string, i: number) => (
                  <li
                    key={i}
                    className={cn(
                      "px-3 py-2 rounded-md text-sm",
                      i === 0 && "bg-teal-100 text-teal-800 font-semibold border-l-3 border-teal-500",
                      i === 1 && "bg-teal-50 text-teal-700 border-l-3 border-teal-400",
                      i === 2 && "bg-emerald-50 text-emerald-700 border-l-3 border-emerald-300",
                      i > 2 && "bg-gray-50 text-gray-600"
                    )}
                  >
                    {type}
                  </li>
                ))}
              </ul>
            ) : (
              <p className="text-gray-500 italic mt-1">None identified</p>
            )}
          </div>

          <div>
            <span className="text-xs font-semibold uppercase tracking-wider text-amber-600">
              Possible Mixed Types
            </span>
            {mixedTypes.length > 0 ? (
              <ul className="mt-2 space-y-1">
                {mixedTypes.map((type: string, i: number) => (
                  <li key={i} className="px-3 py-1.5 bg-gray-50 rounded text-sm text-gray-700">
                    {type}
                  </li>
                ))}
              </ul>
            ) : (
              <p className="text-gray-500 italic mt-1">None identified</p>
            )}
          </div>
        </div>
      )
    }
  } catch {
    // Fall back to raw display
  }

  return (
    <pre className="bg-gray-800 text-gray-100 p-4 rounded-lg text-sm overflow-x-auto">
      {text}
    </pre>
  )
}

/**
 * Modal component for displaying full LLM conversation history
 */
export function ConversationModal({ open, onOpenChange, data }: ConversationModalProps) {
  if (!data) return null

  const sections = parseConversationHistory(data.conversationHistory || '')
  const finalValidator = sections.validators.length > 0 ? sections.validators[sections.validators.length - 1] : ''
  const validationPassed = finalValidator.toUpperCase().includes('VALIDATION PASSED')

  return (
    <Dialog open={open} onOpenChange={onOpenChange}>
      <DialogContent className="max-w-4xl max-h-[90vh] p-0 overflow-hidden">
        {/* Header */}
        <DialogHeader className="bg-gradient-to-r from-teal-600 to-emerald-600 text-white p-6">
          <DialogTitle className="text-xl font-bold">{data.clusterId}</DialogTitle>
          <p className="text-teal-100 mt-1">
            Predicted: <span className="font-semibold text-white">{data.mainType}</span>
          </p>
          {data.score !== undefined && (
            <div className="mt-2">
              <ScoreBadge score={data.score} size="lg" showLabel />
            </div>
          )}
        </DialogHeader>

        <ScrollArea className="h-[calc(90vh-140px)] p-6">
          {/* Full Marker List */}
          <div className="bg-teal-50 rounded-xl p-4 mb-6 border-l-4 border-teal-500">
            <h4 className="text-xs font-semibold uppercase tracking-wider text-teal-600 mb-2">
              Full Marker List ({data.markerCount || 0} markers)
            </h4>
            <p className="text-sm text-gray-700 break-words leading-relaxed">
              {data.markers || 'No markers available'}
            </p>
          </div>

          {/* Merged Groupings */}
          {(data.mergedGrouping1 || data.mergedGrouping2 || data.mergedGrouping3) && (
            <div className="bg-purple-50 rounded-xl p-4 mb-6 border border-purple-100">
              <div className="flex items-center gap-2 mb-3">
                <Users className="h-5 w-5 text-purple-600" />
                <h3 className="font-semibold text-purple-800">Merged Cell Type Groupings</h3>
              </div>
              <div className="space-y-2">
                {data.mergedGrouping1 && (
                  <div className="flex items-center gap-2">
                    <span className="text-xs font-semibold uppercase text-purple-600 w-16">Broad:</span>
                    <span className="px-2 py-1 bg-purple-100 text-purple-800 rounded text-sm">
                      {data.mergedGrouping1}
                    </span>
                  </div>
                )}
                {data.mergedGrouping2 && (
                  <div className="flex items-center gap-2">
                    <span className="text-xs font-semibold uppercase text-purple-600 w-16">Detailed:</span>
                    <span className="px-2 py-1 bg-purple-100 text-purple-800 rounded text-sm">
                      {data.mergedGrouping2}
                    </span>
                  </div>
                )}
                {data.mergedGrouping3 && (
                  <div className="flex items-center gap-2">
                    <span className="text-xs font-semibold uppercase text-purple-600 w-16">Specific:</span>
                    <span className="px-2 py-1 bg-purple-100 text-purple-800 rounded text-sm">
                      {data.mergedGrouping3}
                    </span>
                  </div>
                )}
              </div>
            </div>
          )}

          {/* Annotation Section */}
          <div className="mb-6 rounded-xl border border-teal-100 overflow-hidden">
            <div className="bg-gradient-to-r from-teal-50 to-teal-100/50 px-4 py-3 flex items-center gap-2">
              <Search className="h-5 w-5 text-teal-700" />
              <h3 className="font-semibold text-teal-800">Annotation Analysis</h3>
            </div>
            <div className="p-4 bg-white">
              {formatAnalysisText(sections.annotations.length > 0 ? sections.annotations[sections.annotations.length - 1] : '')}
            </div>
          </div>

          {/* Validator Section */}
          <div className="mb-6 rounded-xl border border-emerald-100 overflow-hidden">
            <div className="bg-gradient-to-r from-emerald-50 to-emerald-100/50 px-4 py-3 flex items-center gap-2">
              <CheckCircle className="h-5 w-5 text-emerald-700" />
              <h3 className="font-semibold text-emerald-800">Validation Check</h3>
            </div>
            <div className="p-4 bg-white">
              {sections.validators.length > 0 ? (
                <>
                  {/* Status badge */}
                  <div className={cn(
                    "inline-flex items-center gap-2 px-3 py-2 rounded-lg mb-3",
                    validationPassed
                      ? "bg-emerald-100 text-emerald-800"
                      : "bg-red-100 text-red-800"
                  )}>
                    <span>{validationPassed ? '✓' : '✗'}</span>
                    <span className="font-semibold">
                      Validation {validationPassed ? 'PASSED' : 'REVIEW NEEDED'}
                    </span>
                  </div>

                  {/* Failed attempts collapsed */}
                  {sections.validators.length > 1 && (
                    <details className="mb-3 p-3 bg-red-50 border-l-4 border-red-400 rounded">
                      <summary className="cursor-pointer text-red-700 font-semibold">
                        {validationPassed
                          ? `⚠️ ${sections.validators.length - 1} failed validation attempt(s) - click to expand`
                          : `⚠️ All ${sections.validators.length} validation attempts failed - click to expand previous attempts`}
                      </summary>
                      <div className="mt-3 space-y-3">
                        {sections.validators.slice(0, -1).map((v, i) => (
                          <div key={i} className="p-3 bg-white rounded border border-red-200">
                            <strong className="text-red-800">Attempt {i + 1}:</strong>
                            <div className="mt-2 text-sm text-gray-600">
                              {formatAnalysisText(v)}
                            </div>
                          </div>
                        ))}
                      </div>
                    </details>
                  )}

                  {/* Final validation with markdown formatting */}
                  <div className="text-gray-600 text-sm">
                    {formatAnalysisText(finalValidator)}
                  </div>
                </>
              ) : (
                <p className="text-gray-500 italic">No validation data available.</p>
              )}
            </div>
          </div>

          {/* Formatting Section */}
          <div className="mb-6 rounded-xl border border-amber-100 overflow-hidden">
            <div className="bg-gradient-to-r from-amber-50 to-amber-100/50 px-4 py-3 flex items-center gap-2">
              <FileText className="h-5 w-5 text-amber-700" />
              <h3 className="font-semibold text-amber-800">Structured Summary</h3>
            </div>
            <div className="p-4 bg-white">
              {sections.formatting ? (
                <FormattingSummary text={sections.formatting} />
              ) : (
                <p className="text-gray-500 italic">No formatting data available.</p>
              )}
            </div>
          </div>

          {/* Scoring Section */}
          {sections.scoring && (
            <div className="mb-6 rounded-xl border border-indigo-100 overflow-hidden">
              <div className="bg-gradient-to-r from-indigo-50 to-indigo-100/50 px-4 py-3 flex items-center gap-2">
                <Star className="h-5 w-5 text-indigo-700" />
                <h3 className="font-semibold text-indigo-800">Quality Score</h3>
              </div>
              <div className="p-4 bg-white">
                {formatAnalysisText(sections.scoring)}
              </div>
            </div>
          )}
        </ScrollArea>
      </DialogContent>
    </Dialog>
  )
}

export default ConversationModal
