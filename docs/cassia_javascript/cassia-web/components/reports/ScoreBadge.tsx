'use client'

import React from 'react'
import { cn } from '@/lib/utils'

export interface ScoreBadgeProps {
  score: number | string | null | undefined
  type?: 'score' | 'confidence'
  size?: 'sm' | 'md' | 'lg'
  showLabel?: boolean
  className?: string
}

/**
 * Color-coded score/confidence badge component
 * Used across all report types for visual quality indicators
 */
export function ScoreBadge({
  score,
  type = 'score',
  size = 'md',
  showLabel = false,
  className
}: ScoreBadgeProps) {
  // Handle null/undefined/empty
  if (score === null || score === undefined || score === '' || score === 'N/A') {
    return (
      <span className={cn(
        "inline-flex items-center rounded-md font-medium",
        "bg-gray-100 text-gray-600",
        size === 'sm' && "px-2 py-0.5 text-xs",
        size === 'md' && "px-2.5 py-1 text-sm",
        size === 'lg' && "px-3 py-1.5 text-base",
        className
      )}>
        {showLabel && <span className="mr-1">Score:</span>}
        N/A
      </span>
    )
  }

  // Parse numeric score
  const numericScore = typeof score === 'string' ? parseFloat(score) : score

  // Handle confidence level strings
  if (type === 'confidence' && typeof score === 'string') {
    const confidenceClasses: Record<string, string> = {
      'high': 'bg-green-100 text-green-800 border-green-200',
      'medium': 'bg-yellow-100 text-yellow-800 border-yellow-200',
      'low': 'bg-red-100 text-red-800 border-red-200',
    }
    const normalizedScore = score.toLowerCase()
    const colorClass = confidenceClasses[normalizedScore] || 'bg-gray-100 text-gray-600'

    return (
      <span className={cn(
        "inline-flex items-center rounded-md font-medium border",
        colorClass,
        size === 'sm' && "px-2 py-0.5 text-xs",
        size === 'md' && "px-2.5 py-1 text-sm",
        size === 'lg' && "px-3 py-1.5 text-base",
        className
      )}>
        {showLabel && <span className="mr-1">Confidence:</span>}
        {score}
      </span>
    )
  }

  // Numeric score styling (0-100 scale)
  if (isNaN(numericScore)) {
    return (
      <span className={cn(
        "inline-flex items-center rounded-md font-medium",
        "bg-gray-100 text-gray-600",
        size === 'sm' && "px-2 py-0.5 text-xs",
        size === 'md' && "px-2.5 py-1 text-sm",
        size === 'lg' && "px-3 py-1.5 text-base",
        className
      )}>
        {showLabel && <span className="mr-1">Score:</span>}
        {score}
      </span>
    )
  }

  // Color based on score thresholds
  let colorClass: string
  if (numericScore >= 90) {
    colorClass = 'bg-gradient-to-r from-green-500 to-emerald-500 text-white'
  } else if (numericScore >= 75) {
    colorClass = 'bg-gradient-to-r from-green-100 to-emerald-100 text-green-800 border border-green-200'
  } else if (numericScore >= 60) {
    colorClass = 'bg-gradient-to-r from-yellow-100 to-amber-100 text-yellow-800 border border-yellow-200'
  } else {
    colorClass = 'bg-gradient-to-r from-red-100 to-rose-100 text-red-800 border border-red-200'
  }

  return (
    <span className={cn(
      "inline-flex items-center rounded-md font-semibold",
      colorClass,
      size === 'sm' && "px-2 py-0.5 text-xs",
      size === 'md' && "px-2.5 py-1 text-sm",
      size === 'lg' && "px-3 py-1.5 text-base",
      className
    )}>
      {showLabel && <span className="mr-1 font-normal">Score:</span>}
      {Math.round(numericScore)}
    </span>
  )
}

export default ScoreBadge
