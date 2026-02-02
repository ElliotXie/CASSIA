'use client'

import React, { useCallback, useState } from 'react'
import { useDropzone } from 'react-dropzone'
import { Upload, File, CheckCircle, AlertCircle, X, RefreshCw, Eye } from 'lucide-react'
import { Card, CardContent } from '@/components/ui/card'
import { Button } from '@/components/ui/button'
import { Badge } from '@/components/ui/badge'
import { processFile, validateMarkerData, type FileData } from '@/lib/utils/file-processing'
import { useAnalysisStore } from '@/lib/stores/analysis-store'
import { GeneConversionPreview, type ConversionOptions } from './GeneConversionPreview'
import {
  generateConversionPreview,
  convertGeneIds,
  applyConversionToData,
  getIdTypeFromDetection,
  getSpeciesFromDetection,
  getConversionSummary,
  type ConversionResult
} from '@/lib/utils/gene-id-converter'
import { getFormatLabel } from '@/lib/utils/gene-id-detection'

interface FileUploadProps {
  onFileProcessed?: (file: File, data: FileData) => void
  showSimpleFormat?: boolean
}

export function FileUpload({ onFileProcessed, showSimpleFormat = false }: FileUploadProps) {
  const [isProcessing, setIsProcessing] = useState(false)
  const [fileData, setFileData] = useState<FileData | null>(null)
  const [error, setError] = useState<string | null>(null)
  const [validationErrors, setValidationErrors] = useState<string[]>([])

  // Gene conversion state
  const [showConversionPreview, setShowConversionPreview] = useState(false)
  const [pendingFile, setPendingFile] = useState<File | null>(null)
  const [pendingFileData, setPendingFileData] = useState<FileData | null>(null)
  const [conversionPreview, setConversionPreview] = useState<ConversionResult[]>([])
  const [isConverting, setIsConverting] = useState(false)
  const [conversionProgress, setConversionProgress] = useState(0)
  const [conversionStage, setConversionStage] = useState<'local' | 'api'>('local')
  const [conversionApplied, setConversionApplied] = useState(false)
  const [showFormatPreview, setShowFormatPreview] = useState(false)

  const { setFile, uploadedFile } = useAnalysisStore()

  // Helper function to finalize file upload
  const finalizeFileUpload = (file: File, data: FileData) => {
    setFileData(data)
    setFile(file, data.data, data)

    if (onFileProcessed) {
      onFileProcessed(file, data)
    }

    console.log('File upload completed successfully')
  }

  // Handle conversion confirmation from preview dialog
  const handleConversionConfirm = async (options: ConversionOptions) => {
    if (!pendingFile || !pendingFileData || !pendingFileData.geneIdDetection) {
      setShowConversionPreview(false)
      return
    }

    if (!options.proceed) {
      // User chose to skip conversion - proceed with original data
      finalizeFileUpload(pendingFile, pendingFileData)
      setShowConversionPreview(false)
      setPendingFile(null)
      setPendingFileData(null)
      return
    }

    // User confirmed conversion
    setIsConverting(true)
    setConversionProgress(0)
    setConversionStage('local')

    try {
      const idType = getIdTypeFromDetection(pendingFileData.geneIdDetection)
      const uniqueGenes = pendingFileData.uniqueGenes || []

      // Perform conversion
      const conversionResult = await convertGeneIds(uniqueGenes, {
        species: options.species,
        idType: idType,
        useApiFallback: options.useApiFallback,
        keepOriginalOnFailure: options.keepOriginalOnFailure,
        onProgress: (current, total, stage) => {
          setConversionProgress(Math.round((current / total) * 100))
          setConversionStage(stage)
        }
      })

      console.log('Conversion result:', getConversionSummary(conversionResult.stats))

      // Apply conversion to data
      if (pendingFileData.geneColumn) {
        const convertedData = applyConversionToData(
          pendingFileData.data,
          pendingFileData.geneColumn,
          conversionResult.results
        )

        // Create new FileData with converted data
        const convertedFileData: FileData = {
          ...pendingFileData,
          data: convertedData
        }

        // Add conversion info to validation warnings
        const conversionMessage = `Gene IDs converted: ${getConversionSummary(conversionResult.stats)}`
        setValidationErrors(prev => [...prev, `Info: ${conversionMessage}`])
        setConversionApplied(true)

        finalizeFileUpload(pendingFile, convertedFileData)
      } else {
        // No gene column found, proceed with original data
        finalizeFileUpload(pendingFile, pendingFileData)
      }

    } catch (err) {
      console.error('Gene conversion error:', err)
      setError(`Gene conversion failed: ${err instanceof Error ? err.message : 'Unknown error'}`)
      // Still proceed with original data on error
      finalizeFileUpload(pendingFile, pendingFileData)
    } finally {
      setIsConverting(false)
      setShowConversionPreview(false)
      setPendingFile(null)
      setPendingFileData(null)
    }
  }

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

      // Check if gene ID conversion is needed
      if (data.geneIdDetection?.needsConversion && data.uniqueGenes && data.uniqueGenes.length > 0) {
        console.log('Gene ID conversion detected, showing preview dialog')

        // Generate preview of first 20 genes
        try {
          const preview = await generateConversionPreview(
            data.uniqueGenes,
            data.geneIdDetection,
            20
          )
          setConversionPreview(preview)
        } catch (previewErr) {
          console.warn('Failed to generate conversion preview:', previewErr)
          setConversionPreview([])
        }

        // Store pending data and show preview dialog
        setPendingFile(file)
        setPendingFileData(data)
        setShowConversionPreview(true)
        setIsProcessing(false)
        return  // Wait for user confirmation
      }

      // No conversion needed - proceed normally
      finalizeFileUpload(file, data)

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
      setConversionApplied(false)
      setPendingFile(null)
      setPendingFileData(null)
      setShowConversionPreview(false)
      useAnalysisStore.getState().reset()
      console.log('File cleared successfully')
    } catch (err) {
      console.error('Error clearing file:', err)
      // Still try to reset UI state even if store reset fails
      setFileData(null)
      setError(null)
      setValidationErrors([])
      setConversionApplied(false)
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
                
                {/* Show gene conversion info if applied */}
                {conversionApplied && fileData.geneIdDetection && (
                  <div className="mt-2 flex items-center space-x-2">
                    <Badge variant="secondary" className="text-xs">
                      <RefreshCw className="h-3 w-3 mr-1" />
                      Gene IDs Converted
                    </Badge>
                    <span className="text-xs text-muted-foreground">
                      from {getFormatLabel(fileData.geneIdDetection.detectedFormat)}
                    </span>
                  </div>
                )}

                {/* Show validation warnings */}
                {validationErrors.length > 0 && (
                  <div className="mt-2 space-y-1">
                    {validationErrors.map((error, index) => (
                      <div key={index} className={`text-xs flex items-center space-x-1 ${
                        error.startsWith('Warning:') ? 'text-amber-600' :
                        error.startsWith('Info:') ? 'text-blue-600' : 'text-red-600'
                      }`}>
                        {error.startsWith('Info:') ? (
                          <CheckCircle className="h-3 w-3" />
                        ) : (
                          <AlertCircle className="h-3 w-3" />
                        )}
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
                      setShowFormatPreview(!showFormatPreview);
                    }}
                  >
                    <Eye className="h-4 w-4 mr-1" />
                    {showFormatPreview ? 'Hide Format' : 'View Example Format'}
                  </Button>
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

        {/* Example Format Preview */}
        {showFormatPreview && (
          <div className="mt-4 p-4 bg-muted/50 border border-muted-foreground/20 rounded-lg space-y-4" onClick={(e) => e.stopPropagation()}>
            {/* Format 1: Seurat/FindAllMarkers */}
            <div>
              <h4 className="text-sm font-semibold mb-1">Format 1: Seurat FindAllMarkers output (recommended)</h4>
              <p className="text-xs text-muted-foreground mb-2">
                Full differential expression output. CASSIA will automatically filter and rank genes.
              </p>
              <div className="overflow-x-auto">
                <table className="w-full text-xs border-collapse">
                  <thead>
                    <tr className="bg-muted">
                      <th className="border border-muted-foreground/20 px-2 py-1 text-left font-semibold">p_val</th>
                      <th className="border border-muted-foreground/20 px-2 py-1 text-left font-semibold bg-amber-100 dark:bg-amber-900/30">avg_log2FC</th>
                      <th className="border border-muted-foreground/20 px-2 py-1 text-left font-semibold">pct.1</th>
                      <th className="border border-muted-foreground/20 px-2 py-1 text-left font-semibold">pct.2</th>
                      <th className="border border-muted-foreground/20 px-2 py-1 text-left font-semibold">p_val_adj</th>
                      <th className="border border-muted-foreground/20 px-2 py-1 text-left font-semibold bg-green-100 dark:bg-green-900/30">cluster</th>
                      <th className="border border-muted-foreground/20 px-2 py-1 text-left font-semibold bg-blue-100 dark:bg-blue-900/30">gene</th>
                    </tr>
                  </thead>
                  <tbody>
                    <tr>
                      <td className="border border-muted-foreground/20 px-2 py-1">1.2e-50</td>
                      <td className="border border-muted-foreground/20 px-2 py-1 bg-amber-50 dark:bg-amber-900/10">2.45</td>
                      <td className="border border-muted-foreground/20 px-2 py-1">0.95</td>
                      <td className="border border-muted-foreground/20 px-2 py-1">0.12</td>
                      <td className="border border-muted-foreground/20 px-2 py-1">3.6e-48</td>
                      <td className="border border-muted-foreground/20 px-2 py-1 bg-green-50 dark:bg-green-900/10">0</td>
                      <td className="border border-muted-foreground/20 px-2 py-1 bg-blue-50 dark:bg-blue-900/10">CD3D</td>
                    </tr>
                    <tr>
                      <td className="border border-muted-foreground/20 px-2 py-1">4.5e-40</td>
                      <td className="border border-muted-foreground/20 px-2 py-1 bg-amber-50 dark:bg-amber-900/10">1.89</td>
                      <td className="border border-muted-foreground/20 px-2 py-1">0.82</td>
                      <td className="border border-muted-foreground/20 px-2 py-1">0.08</td>
                      <td className="border border-muted-foreground/20 px-2 py-1">1.1e-37</td>
                      <td className="border border-muted-foreground/20 px-2 py-1 bg-green-50 dark:bg-green-900/10">0</td>
                      <td className="border border-muted-foreground/20 px-2 py-1 bg-blue-50 dark:bg-blue-900/10">CD3E</td>
                    </tr>
                  </tbody>
                </table>
              </div>
              <p className="text-xs text-muted-foreground mt-1">
                <span className="inline-block w-3 h-3 bg-green-100 dark:bg-green-900/30 border border-muted-foreground/20 mr-1 align-middle rounded-sm"></span> Required: cluster &nbsp;
                <span className="inline-block w-3 h-3 bg-blue-100 dark:bg-blue-900/30 border border-muted-foreground/20 mr-1 align-middle rounded-sm"></span> Required: gene &nbsp;
                <span className="inline-block w-3 h-3 bg-amber-100 dark:bg-amber-900/30 border border-muted-foreground/20 mr-1 align-middle rounded-sm"></span> Required: avg_log2FC (default ranking)
              </p>
            </div>

            {showSimpleFormat && (
              <>
                <hr className="border-muted-foreground/20" />

                {/* Format 2: Simple pre-processed */}
                <div>
                  <h4 className="text-sm font-semibold mb-1">Format 2: Pre-processed marker list</h4>
                  <p className="text-xs text-muted-foreground mb-2">
                    Simple 2-column format with cluster name and comma-separated ordered genes. Used as-is without filtering.
                  </p>
                  <div className="overflow-x-auto">
                    <table className="w-full text-xs border-collapse">
                      <thead>
                        <tr className="bg-muted">
                          <th className="border border-muted-foreground/20 px-2 py-1 text-left font-semibold bg-green-100 dark:bg-green-900/30">cluster</th>
                          <th className="border border-muted-foreground/20 px-2 py-1 text-left font-semibold bg-blue-100 dark:bg-blue-900/30">markers</th>
                        </tr>
                      </thead>
                      <tbody>
                        <tr>
                          <td className="border border-muted-foreground/20 px-2 py-1 bg-green-50 dark:bg-green-900/10">0</td>
                          <td className="border border-muted-foreground/20 px-2 py-1 bg-blue-50 dark:bg-blue-900/10">CD3D, CD3E, IL7R, LEF1, TCF7</td>
                        </tr>
                        <tr>
                          <td className="border border-muted-foreground/20 px-2 py-1 bg-green-50 dark:bg-green-900/10">1</td>
                          <td className="border border-muted-foreground/20 px-2 py-1 bg-blue-50 dark:bg-blue-900/10">MS4A1, CD79A, CD79B, CD19, PAX5</td>
                        </tr>
                      </tbody>
                    </table>
                  </div>
                </div>
              </>
            )}
          </div>
        )}

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

        {/* Gene Conversion Preview Dialog */}
        {showConversionPreview && pendingFileData?.geneIdDetection && (
          <GeneConversionPreview
            isOpen={showConversionPreview}
            onClose={() => {
              setShowConversionPreview(false)
              setPendingFile(null)
              setPendingFileData(null)
            }}
            onConfirm={handleConversionConfirm}
            detectionResult={pendingFileData.geneIdDetection}
            previewData={conversionPreview}
            totalGenes={pendingFileData.geneCount || 0}
            isLoading={isConverting}
            conversionProgress={conversionProgress}
            conversionStage={conversionStage}
          />
        )}
      </CardContent>
    </Card>
  )
}