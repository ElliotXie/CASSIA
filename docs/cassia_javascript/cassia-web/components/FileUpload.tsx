'use client'

import React, { useCallback, useState } from 'react'
import { useDropzone } from 'react-dropzone'
import { Upload, File, CheckCircle, AlertCircle, X } from 'lucide-react'
import { Card, CardContent } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { processFile, validateMarkerData, type FileData } from '@/lib/utils/file-processing'
import { useAnalysisStore } from '@/lib/stores/analysis-store'

interface FileUploadProps {
  onFileProcessed?: (file: File, data: FileData) => void
}

export function FileUpload({ onFileProcessed }: FileUploadProps) {
  const [isProcessing, setIsProcessing] = useState(false)
  const [fileData, setFileData] = useState<FileData | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [validationErrors, setValidationErrors] = useState<string[]>([])
  
  const { setFile, uploadedFile } = useAnalysisStore()

  const loadAndDownloadExample = async (filename: string) => {
    setIsProcessing(true)
    setError(null)
    setValidationErrors([])

    try {
      const response = await fetch(`/examples/${filename}`)
      if (!response.ok) {
        throw new Error(`Failed to fetch example file: ${response.status}`)
      }
      
      const text = await response.text()
      
      // Download the file for the user
      const blob = new Blob([text], { type: 'text/csv' })
      const url = URL.createObjectURL(blob)
      const link = document.createElement('a')
      link.href = url
      link.download = filename
      document.body.appendChild(link)
      link.click()
      document.body.removeChild(link)
      URL.revokeObjectURL(url)
      
      // Create a mock file object using Blob constructor instead of File
      const mockFile = new Blob([text], { type: 'text/csv' }) as any
      mockFile.name = filename
      mockFile.lastModified = Date.now()
      
      // Process the example file
      const data = await processFile(mockFile)
      
      if (!data || typeof data !== 'object') {
        throw new Error('Invalid file data structure returned')
      }
      
      if (!Array.isArray(data.data)) {
        throw new Error('File does not contain valid tabular data')
      }
      
      const validation = validateMarkerData(data.data)
      setValidationErrors(validation.errors || [])
      
      setFileData(data)
      setFile(mockFile, data.data, data)
      
      if (onFileProcessed) {
        onFileProcessed(mockFile, data)
      }
      
      console.log('Example file loaded and downloaded successfully:', filename)
    } catch (err) {
      console.error('Error loading example file:', err)
      setError(err instanceof Error ? err.message : 'Failed to load example file')
    } finally {
      setIsProcessing(false)
    }
  }

  const onDrop = useCallback(async (acceptedFiles: File[]) => {
    const file = acceptedFiles[0]
    if (!file) return

    setIsProcessing(true)
    setError(null)
    setValidationErrors([])

    try {
      console.log('Processing file:', file.name, 'Size:', file.size)
      
      // Process the file
      const data = await processFile(file)
      console.log('File processed successfully:', data)
      
      // Validate that we have the required data structure
      if (!data || typeof data !== 'object') {
        throw new Error('Invalid file data structure returned')
      }
      
      if (!Array.isArray(data.data)) {
        throw new Error('File does not contain valid tabular data')
      }
      
      // Validate the data
      const validation = validateMarkerData(data.data)
      setValidationErrors(validation.errors || [])
      
      if (!validation.isValid) {
        setError('File validation failed. Please check the errors below.')
        return
      }

      // Store in state
      setFileData(data)
      setFile(file, data.data, data)
      
      // Notify parent
      if (onFileProcessed) {
        onFileProcessed(file, data)
      }
      
      console.log('File upload completed successfully')
      
    } catch (err) {
      console.error('File upload error:', err)
      const errorMessage = err instanceof Error ? err.message : 'Failed to process file'
      setError(`Upload failed: ${errorMessage}`)
      
      // Clear any partial data on error
      setFileData(null)
    } finally {
      setIsProcessing(false)
    }
  }, [setFile, onFileProcessed])

  const { getRootProps, getInputProps, isDragActive } = useDropzone({
    onDrop,
    accept: {
      'text/csv': ['.csv'],
      'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet': ['.xlsx'],
      'application/vnd.ms-excel': ['.xls']
    },
    multiple: false,
    disabled: isProcessing
  })

  const clearFile = () => {
    try {
      setFileData(null)
      setError(null)
      setValidationErrors([])
      useAnalysisStore.getState().reset()
      console.log('File cleared successfully')
    } catch (err) {
      console.error('Error clearing file:', err)
      // Still try to reset UI state even if store reset fails
      setFileData(null)
      setError(null)
      setValidationErrors([])
    }
  }

  if (uploadedFile && fileData) {
    return (
      <Card className="w-full">
        <CardContent className="p-6">
          <div className="flex items-start justify-between">
            <div className="flex items-start space-x-3">
              <CheckCircle className="h-5 w-5 text-green-500 mt-0.5" />
              <div className="flex-1">
                <h3 className="font-medium text-sm">{uploadedFile.name}</h3>
                <div className="text-sm text-muted-foreground mt-1">
                  <div className="grid grid-cols-3 gap-4">
                    <div>{fileData.rowCount?.toLocaleString() || 0} rows</div>
                    <div>{fileData.clusterCount || 0} clusters</div>
                    <div>{fileData.geneCount || 0} genes</div>
                  </div>
                </div>
                
                {/* Show preview of headers */}
                {fileData.headers && fileData.headers.length > 0 && (
                  <div className="mt-2 text-xs text-muted-foreground">
                    <strong>Columns:</strong> {fileData.headers.slice(0, 5).join(', ')}
                    {fileData.headers.length > 5 && ` (+${fileData.headers.length - 5} more)`}
                  </div>
                )}
                
                {/* Show validation warnings */}
                {validationErrors.length > 0 && (
                  <div className="mt-2 space-y-1">
                    {validationErrors.map((error, index) => (
                      <div key={index} className={`text-xs flex items-center space-x-1 ${
                        error.startsWith('Warning:') ? 'text-amber-600' : 'text-red-600'
                      }`}>
                        <AlertCircle className="h-3 w-3" />
                        <span>{error}</span>
                      </div>
                    ))}
                  </div>
                )}
              </div>
            </div>
            <Button
              variant="ghost"
              size="sm"
              onClick={clearFile}
              className="ml-2"
            >
              <X className="h-4 w-4" />
            </Button>
          </div>
        </CardContent>
      </Card>
    )
  }

  return (
    <Card className="w-full">
      <CardContent className="p-6">
        <div
          {...getRootProps()}
          className={`border-2 border-dashed rounded-lg p-8 text-center cursor-pointer transition-colors ${
            isDragActive
              ? 'border-primary bg-primary/5'
              : 'border-muted-foreground/25 hover:border-primary/50'
          } ${isProcessing ? 'opacity-50 cursor-not-allowed' : ''}`}
        >
          <input {...getInputProps()} />
          <div className="space-y-4">
            <div className="mx-auto w-12 h-12 bg-primary/10 rounded-lg flex items-center justify-center">
              {isProcessing ? (
                <div className="w-6 h-6 border-2 border-primary border-t-transparent rounded-full animate-spin" />
              ) : (
                <Upload className="h-6 w-6 text-primary" />
              )}
            </div>
            
            <div className="space-y-2">
              <h3 className="font-medium">
                {isProcessing
                  ? 'Processing file...'
                  : isDragActive
                  ? 'Drop your file here'
                  : 'Upload marker data'
                }
              </h3>
              <p className="text-sm text-muted-foreground">
                {isProcessing
                  ? 'Please wait while we process your file'
                  : 'Drag & drop a CSV or XLSX file, or click to browse'
                }
              </p>
            </div>
            
            {!isProcessing && (
              <>
                <div className="text-xs text-muted-foreground space-y-1">
                  <div>📄 Supported formats: CSV, XLSX</div>
                  <div>📊 Expected: FindAllMarkers output with cluster and gene columns</div>
                </div>
                
                <div className="flex gap-2 justify-center pt-2">
                  <Button 
                    variant="outline" 
                    size="sm" 
                    onClick={(e) => {
                      e.stopPropagation();
                      loadAndDownloadExample('batch_raw_seurat_example.csv');
                    }}
                    disabled={isProcessing}
                  >
                    {isProcessing ? 'Loading...' : 'Load & Download Example'}
                  </Button>
                </div>
              </>
            )}
          </div>
        </div>
        
        {error && (
          <div className="mt-4 p-3 bg-red-50 border border-red-200 rounded-md">
            <div className="flex items-center space-x-2 text-red-800">
              <AlertCircle className="h-4 w-4" />
              <span className="text-sm font-medium">Error</span>
            </div>
            <p className="text-sm text-red-700 mt-1">{error}</p>
            
            {validationErrors.length > 0 && (
              <ul className="mt-2 text-sm text-red-700 space-y-1">
                {validationErrors.map((error, index) => (
                  <li key={index} className="flex items-center space-x-1">
                    <span>•</span>
                    <span>{error}</span>
                  </li>
                ))}
              </ul>
            )}
          </div>
        )}
      </CardContent>
    </Card>
  )
}