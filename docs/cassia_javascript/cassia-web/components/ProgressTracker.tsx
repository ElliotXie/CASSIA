'use client'

import React, { useEffect, useRef } from 'react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Progress } from '@/components/ui/progress'
import { Button } from '@/components/ui/button'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Terminal, Download, X, CheckCircle, Clock } from 'lucide-react'
import { useAnalysisStore } from '@/lib/stores/analysis-store'

export function ProgressTracker() {
  const { 
    isRunning, 
    progress, 
    currentStep, 
    logs, 
    downloadLinks,
    reset,
    fileMetadata 
  } = useAnalysisStore()
  
  const logsEndRef = useRef<HTMLDivElement>(null)

  // Calculate estimated time based on cluster count
  const getEstimatedTime = (clusterCount: number, currentProgress: number) => {
    if (!clusterCount || currentProgress <= 0) return 'Calculating...'
    
    // Base time is 30s for 5 clusters, +5s per additional cluster, max 90s at 30 clusters
    const baseTime = 30 // seconds for 5 clusters
    const additionalTime = Math.max(0, clusterCount - 5) * 5 // 5s per additional cluster
    const maxTime = 90 // maximum 1min 30s
    const totalEstimatedTime = Math.min(baseTime + additionalTime, maxTime)
    
    // Calculate remaining time based on progress
    const remainingTime = Math.max(0, totalEstimatedTime * (100 - currentProgress) / 100)
    
    if (remainingTime < 60) {
      return `${Math.round(remainingTime)}s`
    } else {
      const minutes = Math.floor(remainingTime / 60)
      const seconds = Math.round(remainingTime % 60)
      return seconds > 0 ? `${minutes}m ${seconds}s` : `${minutes}m`
    }
  }

  // Auto-scroll logs to bottom
  useEffect(() => {
    logsEndRef.current?.scrollIntoView({ behavior: 'smooth' })
  }, [logs])

  if (!isRunning && logs.length === 0) {
    return null
  }

  return (
    <Card className="w-full">
      <CardHeader>
        <div className="flex items-center justify-between">
          <CardTitle className="flex items-center space-x-2">
            <Terminal className="h-5 w-5" />
            <span>Analysis Progress</span>
          </CardTitle>
          {!isRunning && (
            <Button variant="ghost" size="sm" onClick={reset}>
              <X className="h-4 w-4" />
            </Button>
          )}
        </div>
      </CardHeader>
      <CardContent className="space-y-4">
        {/* Progress Bar */}
        {isRunning && (
          <div className="space-y-2">
            <div className="flex items-center justify-between text-sm">
              <span className="font-medium">{currentStep}</span>
              <span className="text-muted-foreground">{Math.round(progress)}%</span>
            </div>
            <Progress value={progress} className="h-2" />
          </div>
        )}

        {/* Status */}
        {!isRunning && progress === 100 && (
          <div className="flex items-center space-x-2 text-green-600">
            <CheckCircle className="h-5 w-5" />
            <span className="font-medium">Analysis Complete!</span>
          </div>
        )}

        {/* Download Links */}
        {downloadLinks.length > 0 && (
          <div className="space-y-2">
            <h4 className="font-medium text-sm flex items-center space-x-2">
              <Download className="h-4 w-4" />
              <span>Download Results</span>
            </h4>
            <div className="grid grid-cols-2 gap-2">
              {downloadLinks.map((link, index) => (
                <Button key={index} variant="outline" size="sm" asChild>
                  <a href={link} download>
                    Download {index + 1}
                  </a>
                </Button>
              ))}
            </div>
          </div>
        )}

        {/* Console Logs */}
        <div className="space-y-2">
          <h4 className="font-medium text-sm flex items-center space-x-2">
            <Terminal className="h-4 w-4" />
            <span>Console Output</span>
          </h4>
          <Card className="bg-gray-50 dark:bg-gray-900 border-2">
            <CardContent className="p-4">
              <ScrollArea className="h-64 w-full">
                <div className="font-mono text-xs space-y-1">
                  {logs.length === 0 ? (
                    <div className="text-gray-500 dark:text-gray-400">No logs yet...</div>
                  ) : (
                    logs.map((log, index) => (
                      <div key={index} className="flex items-start space-x-2">
                        <span className="text-gray-500 dark:text-gray-400 shrink-0">
                          {new Date().toLocaleTimeString()}
                        </span>
                        <span className="flex-1 break-all text-gray-900 dark:text-gray-100">{log}</span>
                      </div>
                    ))
                  )}
                  <div ref={logsEndRef} />
                </div>
              </ScrollArea>
            </CardContent>
          </Card>
        </div>

        {/* Time Estimate */}
        {isRunning && (
          <div className="flex items-center space-x-2 text-sm text-muted-foreground">
            <Clock className="h-4 w-4" />
            <span>
              Estimated time remaining: {getEstimatedTime(fileMetadata?.clusterCount || 5, progress)}
            </span>
          </div>
        )}
      </CardContent>
    </Card>
  )
}