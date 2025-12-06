'use client'

import React, { useState, useEffect } from 'react'
import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Badge } from '@/components/ui/badge'
import { 
  Download, 
  FileText, 
  Database, 
  Archive, 
  CheckCircle2, 
  FileDown,
  Folder,
  BarChart3,
  Zap
} from 'lucide-react'
import { createDownloadUrls, cleanupDownloadUrls } from '@/lib/cassia/runCASSIA_pipeline'

interface ResultsDownloaderProps {
  finalResults: any
  isVisible: boolean
}

interface FileGroup {
  title: string
  description: string
  icon: React.ReactNode
  files: Array<{
    key: string
    name: string
    size: number
    type: string
    downloadUrl?: string
  }>
}

export function ResultsDownloader({ finalResults, isVisible }: ResultsDownloaderProps) {
  const [downloadUrls, setDownloadUrls] = useState<any>(null)
  const [downloadingFiles, setDownloadingFiles] = useState<Set<string>>(new Set())

  // Create download URLs when finalResults is available
  useEffect(() => {
    if (finalResults && isVisible) {
      const urls = createDownloadUrls(finalResults)
      setDownloadUrls(urls)
      
      // Cleanup function
      return () => {
        if (urls) {
          cleanupDownloadUrls(urls)
        }
      }
    }
  }, [finalResults, isVisible])

  const handleDownload = async (fileKey: string, filename: string, url: string) => {
    setDownloadingFiles(prev => new Set([...prev, fileKey]))
    
    try {
      // Create a temporary link and trigger download
      const link = document.createElement('a')
      link.href = url
      link.download = filename
      document.body.appendChild(link)
      link.click()
      document.body.removeChild(link)
      
      // Small delay to show the downloading state
      await new Promise(resolve => setTimeout(resolve, 1000))
    } catch (error) {
      console.error('Download failed:', error)
    } finally {
      setDownloadingFiles(prev => {
        const newSet = new Set(prev)
        newSet.delete(fileKey)
        return newSet
      })
    }
  }

  const handleDownloadAll = async () => {
    if (!downloadUrls || !finalResults) return

    // Download all files individually (ZIP creation would require additional library)
    const allFiles = Object.entries(downloadUrls.individual)
    
    for (const [key, fileInfo] of allFiles) {
      await handleDownload(key, fileInfo.filename, fileInfo.url)
      // Add small delay between downloads
      await new Promise(resolve => setTimeout(resolve, 500))
    }
  }

  const formatFileSize = (bytes: number): string => {
    if (bytes === 0) return '0 Bytes'
    const k = 1024
    const sizes = ['Bytes', 'KB', 'MB', 'GB']
    const i = Math.floor(Math.log(bytes) / Math.log(k))
    return parseFloat((bytes / Math.pow(k, i)).toFixed(2)) + ' ' + sizes[i]
  }

  const getFileIcon = (type: string) => {
    if (type.includes('csv')) return <Database className="h-4 w-4 text-green-600" />
    if (type.includes('html')) return <FileText className="h-4 w-4 text-blue-600" />
    if (type.includes('text')) return <FileDown className="h-4 w-4 text-gray-600" />
    return <FileDown className="h-4 w-4 text-gray-600" />
  }

  const getBoostFileLabel = (fileName: string, type: string) => {
    const safeName = (fileName || '').toLowerCase()
    const safeType = (type || '').toLowerCase()

    if (safeName.includes('conversation') || safeName.includes('history') || safeType.includes('jsonl')) {
      return 'Conversation history'
    }
    if (safeName.includes('report') || safeName.includes('summary') || safeType.includes('html')) {
      return 'HTML report'
    }
    return 'Boost output'
  }

  if (!finalResults || !isVisible || !downloadUrls) {
    return null
  }

  // Organize files into logical groups
  const fileGroups: FileGroup[] = [
    {
      title: 'Core Results',
      description: 'Main annotation and scoring results',
      icon: <BarChart3 className="h-5 w-5 text-blue-600" />,
      files: [
        ...(finalResults.files.annotation ? [{
          key: 'annotation',
          name: finalResults.files.annotation.filename,
          size: finalResults.files.annotation.size,
          type: finalResults.files.annotation.type,
          downloadUrl: downloadUrls.individual.annotation?.url
        }] : []),
        ...(finalResults.files.scoring ? [{
          key: 'scoring',
          name: finalResults.files.scoring.filename,
          size: finalResults.files.scoring.size,
          type: finalResults.files.scoring.type,
          downloadUrl: downloadUrls.individual.scoring?.url
        }] : [])
      ]
    },
    {
      title: 'Analysis Reports',
      description: 'HTML reports and summaries',
      icon: <FileText className="h-5 w-5 text-purple-600" />,
      files: Object.entries(finalResults.files.reports || {}).map(([key, file]: [string, any]) => ({
        key: `report_${key}`,
        name: file.filename,
        size: file.size,
        type: file.type,
        downloadUrl: downloadUrls.individual[key]?.url
      }))
    },
    {
      title: 'Annotation Boost',
      description: 'Enhanced analysis for low-scoring clusters',
      icon: <Zap className="h-5 w-5 text-orange-600" />,
      files: Object.entries(finalResults.files.boost || {}).map(([key, file]: [string, any]) => ({
        key: `boost_${key}`,
        name: file.filename,
        size: file.size,
        type: file.type,
        downloadUrl: downloadUrls.individual[key]?.url
      }))
    }
  ].filter(group => group.files.length > 0)

  const totalFiles = fileGroups.reduce((sum, group) => sum + group.files.length, 0)
  const totalSize = fileGroups.reduce((sum, group) => 
    sum + group.files.reduce((groupSum, file) => groupSum + file.size, 0), 0
  )

  return (
    <Card className="w-full">
      <CardHeader>
        <div className="flex items-center justify-between">
          <CardTitle className="flex items-center space-x-2">
            <Download className="h-5 w-5" />
            <span>Download Results</span>
          </CardTitle>
          <Badge variant="secondary" className="flex items-center space-x-1">
            <Archive className="h-3 w-3" />
            <span>{totalFiles} files ({formatFileSize(totalSize)})</span>
          </Badge>
        </div>
      </CardHeader>
      <CardContent className="space-y-6">
        {/* Summary */}
        <div className="bg-muted/50 rounded-lg p-4">
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
            <div>
              <div className="font-medium text-muted-foreground">Total Clusters</div>
              <div className="text-lg font-bold">{finalResults.summary.totalClusters}</div>
            </div>
            <div>
              <div className="font-medium text-muted-foreground">Low Score Clusters</div>
              <div className="text-lg font-bold text-orange-600">{finalResults.summary.lowScoreClusters}</div>
            </div>
            <div>
              <div className="font-medium text-muted-foreground">Total Files</div>
              <div className="text-lg font-bold">{totalFiles}</div>
            </div>
            <div>
              <div className="font-medium text-muted-foreground">Total Size</div>
              <div className="text-lg font-bold">{formatFileSize(totalSize)}</div>
            </div>
          </div>
        </div>

        {/* Download All Button */}
        <div className="flex justify-center">
          <Button 
            onClick={handleDownloadAll}
            className="flex items-center space-x-2"
            size="lg"
          >
            <Download className="h-4 w-4" />
            <span>Download All Files</span>
          </Button>
        </div>

        {/* File Groups */}
        <div className="space-y-4">
          {fileGroups.map((group, groupIndex) => (
            <div key={groupIndex} className="border rounded-lg p-4">
              <div className="flex items-center space-x-3 mb-3">
                {group.icon}
                <div>
                  <h3 className="font-semibold">{group.title}</h3>
                  <p className="text-sm text-muted-foreground">{group.description}</p>
                </div>
                <Badge variant="outline" className="ml-auto">
                  {group.files.length} files
                </Badge>
              </div>

              {group.title === 'Annotation Boost' && (
                <div className="text-xs text-muted-foreground bg-muted/40 rounded p-2 mb-2">
                  <span className="font-semibold">What’s inside: </span>
                  HTML report = detailed results per cluster · Conversation history = full chat transcript used for boosting
                </div>
              )}
              
              <div className="grid grid-cols-1 md:grid-cols-2 gap-2">
                {group.files.map((file) => (
                  <div
                    key={file.key}
                    className="flex items-center justify-between p-3 bg-muted/30 rounded border"
                  >
                    <div className="flex items-center space-x-3 flex-1 min-w-0">
                      {getFileIcon(file.type)}
                      <div className="min-w-0 flex-1">
                        <div className="font-medium text-sm truncate">{file.name}</div>
                        <div className="flex items-center flex-wrap gap-2 text-xs text-muted-foreground">
                          <span>{formatFileSize(file.size)}</span>
                          {group.title === 'Annotation Boost' && (
                            <Badge variant="secondary" className="text-[11px] px-2 py-0">
                              {getBoostFileLabel(file.name, file.type)}
                            </Badge>
                          )}
                        </div>
                      </div>
                    </div>
                    
                    <Button
                      variant="outline"
                      size="sm"
                      onClick={() => file.downloadUrl && handleDownload(file.key, file.name, file.downloadUrl)}
                      disabled={!file.downloadUrl || downloadingFiles.has(file.key)}
                      className="ml-2 shrink-0"
                    >
                      {downloadingFiles.has(file.key) ? (
                        <CheckCircle2 className="h-4 w-4 text-green-600" />
                      ) : (
                        <Download className="h-4 w-4" />
                      )}
                    </Button>
                  </div>
                ))}
              </div>
            </div>
          ))}
        </div>

        {/* Analysis Configuration */}
        <div className="border-t pt-4">
          <h4 className="font-medium mb-3 flex items-center space-x-2">
            <Folder className="h-4 w-4" />
            <span>Analysis Configuration</span>
          </h4>
          <div className="bg-muted/30 rounded p-3 text-sm space-y-2">
            <div className="grid grid-cols-2 gap-4">
              <div>
                <span className="font-medium">Tissue:</span> {finalResults.summary.tissue}
              </div>
              <div>
                <span className="font-medium">Species:</span> {finalResults.summary.species}
              </div>
              <div>
                <span className="font-medium">Score Threshold:</span> {finalResults.metadata.scoreThreshold}
              </div>
              <div>
                <span className="font-medium">Max Workers:</span> {finalResults.metadata.configuration.maxWorkers}
              </div>
            </div>
            <div className="pt-2 border-t">
              <div className="font-medium mb-1">Models Used:</div>
              <div className="grid grid-cols-1 gap-1 text-xs">
                <div>Annotation: {finalResults.metadata.models.annotation.model}</div>
                <div>Scoring: {finalResults.metadata.models.scoring.model}</div>
                <div>Boost: {finalResults.metadata.models.annotationBoost.model}</div>
              </div>
            </div>
          </div>
        </div>
      </CardContent>
    </Card>
  )
}