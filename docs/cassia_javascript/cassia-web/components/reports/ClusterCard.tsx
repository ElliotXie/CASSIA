'use client'

import React from 'react'
import { Eye, Layers, Activity, Globe, Beaker } from 'lucide-react'
import { Button } from '@/components/ui/button'
import { ScoreBadge } from './ScoreBadge'
import { ClusterData } from './ConversationModal'
import { cn } from '@/lib/utils'

export interface ClusterCardProps {
  data: ClusterData
  onClick?: () => void
  className?: string
}

/**
 * Truncate text with ellipsis
 */
function truncate(text: string, maxLength: number = 100): string {
  if (!text) return ''
  if (text.length <= maxLength) return text
  return text.slice(0, maxLength) + '...'
}

/**
 * Card component for displaying individual cluster information
 * Part of the batch report view
 */
export function ClusterCard({ data, onClick, className }: ClusterCardProps) {
  const subTypes = data.subTypes || []
  const mixedTypes = data.mixedTypes || []

  return (
    <div
      className={cn(
        "bg-white/90 backdrop-blur-sm border border-teal-100 rounded-2xl p-5",
        "shadow-sm hover:shadow-lg hover:border-teal-300",
        "transition-all duration-300 hover:-translate-y-1",
        "flex flex-col gap-4",
        className
      )}
    >
      {/* Header */}
      <div className="flex justify-between items-start gap-3">
        <h3 className="font-bold text-gray-800 text-lg leading-tight">
          {data.clusterId}
        </h3>
        <div className="flex gap-1.5 flex-wrap justify-end">
          <span className="px-2 py-0.5 bg-gradient-to-r from-teal-500 to-emerald-500 text-white text-xs font-semibold rounded-md">
            {data.markerCount || 0} markers
          </span>
          {data.iterations && (
            <span className="px-2 py-0.5 bg-gradient-to-r from-amber-400 to-orange-400 text-white text-xs font-semibold rounded-md">
              {data.iterations} iter
            </span>
          )}
          {data.score !== undefined && data.score !== null && data.score !== 'N/A' && (
            <ScoreBadge score={data.score} size="sm" />
          )}
        </div>
      </div>

      {/* Main Type */}
      <div className="bg-gradient-to-r from-teal-50 to-emerald-50 rounded-xl p-4 border-l-4 border-teal-500">
        <span className="text-xs font-semibold uppercase tracking-wider text-teal-600 block mb-1">
          Predicted
        </span>
        <span className="text-lg font-bold text-teal-800">
          {data.mainType || 'Unknown'}
        </span>
      </div>

      {/* Sub Cell Types */}
      <div>
        <span className="text-xs font-semibold uppercase tracking-wider text-gray-500 block mb-2">
          Sub Cell Types
        </span>
        {subTypes.length > 0 ? (
          <ul className="space-y-1.5">
            {subTypes.slice(0, 3).map((type, i) => (
              <li
                key={i}
                className={cn(
                  "px-3 py-1.5 rounded-md text-sm",
                  i === 0 && "bg-teal-100/80 text-teal-800 font-semibold border-l-3 border-teal-500",
                  i === 1 && "bg-teal-50 text-teal-700 border-l-3 border-teal-400",
                  i === 2 && "bg-emerald-50/70 text-emerald-700 border-l-3 border-emerald-300",
                )}
              >
                {type}
              </li>
            ))}
            {subTypes.length > 3 && (
              <li className="text-gray-400 text-sm italic pl-3">
                +{subTypes.length - 3} more
              </li>
            )}
          </ul>
        ) : (
          <span className="text-gray-400 italic text-sm">None identified</span>
        )}
      </div>

      {/* Merged Groupings */}
      {(data.mergedGrouping1 || data.mergedGrouping2 || data.mergedGrouping3) && (
        <div className="mt-1">
          <span className="text-xs font-semibold uppercase tracking-wider text-gray-500 block mb-2">
            Merged Cell Type Groups
          </span>
          <div className="flex flex-col gap-1.5">
            {data.mergedGrouping1 && (
              <div className="flex items-center gap-2">
                <span className="text-[10px] font-bold uppercase text-gray-400 w-12">Broad:</span>
                <span className="px-2 py-0.5 bg-teal-100/70 text-teal-700 rounded text-xs border border-teal-200">
                  {truncate(data.mergedGrouping1, 40)}
                </span>
              </div>
            )}
            {data.mergedGrouping2 && (
              <div className="flex items-center gap-2">
                <span className="text-[10px] font-bold uppercase text-gray-400 w-12">Detail:</span>
                <span className="px-2 py-0.5 bg-emerald-100/70 text-emerald-700 rounded text-xs border border-emerald-200">
                  {truncate(data.mergedGrouping2, 40)}
                </span>
              </div>
            )}
            {data.mergedGrouping3 && (
              <div className="flex items-center gap-2">
                <span className="text-[10px] font-bold uppercase text-gray-400 w-12">Specific:</span>
                <span className="px-2 py-0.5 bg-amber-100/70 text-amber-700 rounded text-xs border border-amber-200">
                  {truncate(data.mergedGrouping3, 40)}
                </span>
              </div>
            )}
          </div>
        </div>
      )}

      {/* Mixed Types */}
      {mixedTypes.length > 0 && (
        <div className="flex items-start gap-2 text-sm text-gray-600">
          <span className="text-gray-400 mt-0.5">+</span>
          <div>
            <span className="font-medium">Mixed Types:</span>{' '}
            <span className="text-gray-500">{truncate(mixedTypes.join(', '), 60)}</span>
          </div>
        </div>
      )}

      {/* Metadata */}
      <div className="bg-gray-50/80 rounded-lg p-3 space-y-2">
        <div className="flex flex-wrap gap-x-4 gap-y-1">
          <span className="flex items-center gap-1.5 text-xs text-gray-500">
            <Layers className="h-3 w-3 text-teal-500" />
            {data.model || 'N/A'}
          </span>
          <span className="flex items-center gap-1.5 text-xs text-gray-500">
            <Activity className="h-3 w-3 text-teal-500" />
            {data.provider || 'N/A'}
          </span>
        </div>
        <div className="flex flex-wrap gap-x-4 gap-y-1">
          <span className="flex items-center gap-1.5 text-xs text-gray-500">
            <Beaker className="h-3 w-3 text-teal-500" />
            {data.tissue || 'N/A'}
          </span>
          <span className="flex items-center gap-1.5 text-xs text-gray-500">
            <Globe className="h-3 w-3 text-teal-500" />
            {data.species || 'N/A'}
          </span>
        </div>
      </div>

      {/* Marker Preview */}
      <div className="text-xs text-gray-500">
        <span className="font-medium">Markers:</span>{' '}
        <span className="text-gray-400">{truncate(data.markers || '', 80)}</span>
      </div>

      {/* View Analysis Button */}
      {onClick && (
        <Button
          onClick={onClick}
          className="w-full mt-auto bg-gradient-to-r from-teal-600 to-emerald-600 hover:from-teal-700 hover:to-emerald-700 text-white font-semibold"
        >
          <Eye className="h-4 w-4 mr-2" />
          View Full Analysis
        </Button>
      )}
    </div>
  )
}

export default ClusterCard
