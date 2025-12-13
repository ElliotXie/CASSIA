'use client'

import React, { useState } from 'react'
import {
  Dialog,
  DialogContent,
  DialogHeader,
  DialogTitle,
  DialogDescription,
  DialogFooter
} from '@/components/ui/dialog'
import { Button } from '@/components/ui/button'
import { Badge } from '@/components/ui/badge'
import { Progress } from '@/components/ui/progress'
import { ScrollArea } from '@/components/ui/scroll-area'
import { Card, CardContent } from '@/components/ui/card'
import {
  CheckCircle,
  XCircle,
  AlertTriangle,
  ArrowRight,
  Loader2,
  RefreshCw,
  Globe
} from 'lucide-react'
import type { GeneIdDetectionResult } from '@/lib/utils/gene-id-detection'
import type { ConversionResult } from '@/lib/utils/gene-id-converter'
import { getFormatLabel } from '@/lib/utils/gene-id-detection'

export interface ConversionOptions {
  proceed: boolean
  species: 'human' | 'mouse'
  useApiFallback: boolean
  keepOriginalOnFailure: boolean
}

interface GeneConversionPreviewProps {
  isOpen: boolean
  onClose: () => void
  onConfirm: (options: ConversionOptions) => void
  detectionResult: GeneIdDetectionResult
  previewData: ConversionResult[]
  totalGenes: number
  isLoading: boolean
  conversionProgress?: number
  conversionStage?: 'local' | 'api'
}

export function GeneConversionPreview({
  isOpen,
  onClose,
  onConfirm,
  detectionResult,
  previewData,
  totalGenes,
  isLoading,
  conversionProgress,
  conversionStage
}: GeneConversionPreviewProps) {
  const [species, setSpecies] = useState<'human' | 'mouse'>(
    detectionResult.suggestedSpecies || 'human'
  )
  const [useApiFallback, setUseApiFallback] = useState(true)
  const [keepOriginalOnFailure, setKeepOriginalOnFailure] = useState(true)

  // Calculate preview statistics
  const previewStats = {
    converted: previewData.filter(r => r.status === 'converted' || r.status === 'api_fallback').length,
    notFound: previewData.filter(r => r.status === 'not_found').length,
    total: previewData.length
  }

  const estimatedSuccessRate = previewStats.total > 0
    ? Math.round((previewStats.converted / previewStats.total) * 100)
    : 0

  const handleConfirm = () => {
    onConfirm({
      proceed: true,
      species,
      useApiFallback,
      keepOriginalOnFailure
    })
  }

  const handleSkip = () => {
    onConfirm({
      proceed: false,
      species,
      useApiFallback,
      keepOriginalOnFailure
    })
  }

  // Get confidence badge variant
  const getConfidenceBadge = () => {
    if (detectionResult.confidence >= 0.95) {
      return { variant: 'default' as const, text: 'High confidence' }
    } else if (detectionResult.confidence >= 0.8) {
      return { variant: 'secondary' as const, text: 'Good confidence' }
    } else {
      return { variant: 'outline' as const, text: 'Low confidence' }
    }
  }

  const confidenceBadge = getConfidenceBadge()

  return (
    <Dialog open={isOpen} onOpenChange={(open) => !open && onClose()}>
      <DialogContent className="max-w-2xl max-h-[90vh] overflow-hidden flex flex-col">
        <DialogHeader>
          <DialogTitle className="flex items-center space-x-2">
            <RefreshCw className="h-5 w-5 text-blue-500" />
            <span>Gene ID Conversion Available</span>
          </DialogTitle>
          <DialogDescription>
            Your file contains gene identifiers that can be converted to gene symbols for better annotation results.
          </DialogDescription>
        </DialogHeader>

        {/* Detection Summary */}
        <Card className="bg-blue-50 dark:bg-blue-950 border-blue-200 dark:border-blue-800">
          <CardContent className="p-4">
            <div className="flex items-center justify-between">
              <div>
                <div className="font-medium text-blue-900 dark:text-blue-100">
                  Detected: {getFormatLabel(detectionResult.detectedFormat)}
                </div>
                <div className="text-sm text-blue-700 dark:text-blue-300">
                  {totalGenes.toLocaleString()} unique genes found
                </div>
              </div>
              <Badge variant={confidenceBadge.variant}>
                {Math.round(detectionResult.confidence * 100)}% {confidenceBadge.text}
              </Badge>
            </div>
          </CardContent>
        </Card>

        {/* Species Selection (for Entrez IDs or when needed) */}
        {(detectionResult.detectedFormat === 'entrez' || detectionResult.detectedFormat === 'mixed') && (
          <div className="space-y-2">
            <label className="text-sm font-medium">Select Species for Conversion</label>
            <div className="flex space-x-4">
              <label className="flex items-center space-x-2 cursor-pointer">
                <input
                  type="radio"
                  name="species"
                  value="human"
                  checked={species === 'human'}
                  onChange={() => setSpecies('human')}
                  className="w-4 h-4"
                  disabled={isLoading}
                />
                <span>Human (HGNC symbols)</span>
              </label>
              <label className="flex items-center space-x-2 cursor-pointer">
                <input
                  type="radio"
                  name="species"
                  value="mouse"
                  checked={species === 'mouse'}
                  onChange={() => setSpecies('mouse')}
                  className="w-4 h-4"
                  disabled={isLoading}
                />
                <span>Mouse (MGI symbols)</span>
              </label>
            </div>
          </div>
        )}

        {/* Preview Table */}
        <div className="flex-1 min-h-0">
          <div className="text-sm font-medium mb-2 flex items-center justify-between">
            <span>Conversion Preview (first {previewData.length} genes)</span>
            {previewStats.converted > 0 && (
              <span className="text-green-600 dark:text-green-400 text-xs">
                {previewStats.converted}/{previewStats.total} found in local mapping
              </span>
            )}
          </div>
          <ScrollArea className="h-48 border rounded-lg">
            <table className="w-full text-sm">
              <thead className="bg-muted sticky top-0">
                <tr>
                  <th className="px-3 py-2 text-left font-medium">Original ID</th>
                  <th className="px-3 py-2 text-center w-8"></th>
                  <th className="px-3 py-2 text-left font-medium">Gene Symbol</th>
                  <th className="px-3 py-2 text-center font-medium w-20">Status</th>
                </tr>
              </thead>
              <tbody>
                {previewData.map((result, index) => (
                  <tr key={index} className="border-t hover:bg-muted/50">
                    <td className="px-3 py-2 font-mono text-xs truncate max-w-[200px]" title={result.original}>
                      {result.original}
                    </td>
                    <td className="px-3 py-2 text-center">
                      <ArrowRight className="h-4 w-4 text-muted-foreground" />
                    </td>
                    <td className="px-3 py-2 font-mono text-xs">
                      {result.converted || (
                        <span className="text-muted-foreground italic">—</span>
                      )}
                    </td>
                    <td className="px-3 py-2 text-center">
                      {result.status === 'converted' || result.status === 'api_fallback' ? (
                        <CheckCircle className="h-4 w-4 text-green-500 mx-auto" />
                      ) : (
                        <XCircle className="h-4 w-4 text-red-400 mx-auto" />
                      )}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </ScrollArea>
        </div>

        {/* Statistics */}
        <div className="grid grid-cols-3 gap-4 text-center">
          <div className="p-3 bg-green-50 dark:bg-green-950 rounded-lg border border-green-200 dark:border-green-800">
            <div className="text-2xl font-bold text-green-600 dark:text-green-400">{previewStats.converted}</div>
            <div className="text-xs text-green-700 dark:text-green-300">Found</div>
          </div>
          <div className="p-3 bg-amber-50 dark:bg-amber-950 rounded-lg border border-amber-200 dark:border-amber-800">
            <div className="text-2xl font-bold text-amber-600 dark:text-amber-400">{previewStats.notFound}</div>
            <div className="text-xs text-amber-700 dark:text-amber-300">Not Found</div>
          </div>
          <div className="p-3 bg-blue-50 dark:bg-blue-950 rounded-lg border border-blue-200 dark:border-blue-800">
            <div className="text-2xl font-bold text-blue-600 dark:text-blue-400">{estimatedSuccessRate}%</div>
            <div className="text-xs text-blue-700 dark:text-blue-300">Est. Success</div>
          </div>
        </div>

        {/* Options */}
        <div className="space-y-3 p-4 bg-muted/50 rounded-lg">
          <div className="flex items-center justify-between">
            <label className="flex items-center space-x-2 cursor-pointer">
              <input
                type="checkbox"
                checked={useApiFallback}
                onChange={(e) => setUseApiFallback(e.target.checked)}
                className="w-4 h-4 rounded"
                disabled={isLoading}
              />
              <span className="text-sm flex items-center gap-1">
                <Globe className="h-3 w-3" />
                Try online lookup for unmapped genes
              </span>
            </label>
            <span className="text-xs text-muted-foreground">mygene.info</span>
          </div>
          <div className="flex items-center space-x-2">
            <input
              type="checkbox"
              checked={keepOriginalOnFailure}
              onChange={(e) => setKeepOriginalOnFailure(e.target.checked)}
              className="w-4 h-4 rounded"
              disabled={isLoading}
            />
            <span className="text-sm">Keep original ID if conversion fails</span>
          </div>
        </div>

        {/* Progress indicator during conversion */}
        {isLoading && conversionProgress !== undefined && (
          <div className="space-y-2">
            <div className="flex justify-between text-sm">
              <span>
                {conversionStage === 'api'
                  ? 'Querying online database...'
                  : 'Converting genes...'}
              </span>
              <span>{conversionProgress}%</span>
            </div>
            <Progress value={conversionProgress} className="h-2" />
          </div>
        )}

        {/* Warning for low success rate */}
        {estimatedSuccessRate < 50 && previewStats.total > 0 && !isLoading && (
          <div className="flex items-start space-x-2 p-3 bg-amber-50 dark:bg-amber-950 border border-amber-200 dark:border-amber-800 rounded-lg">
            <AlertTriangle className="h-4 w-4 text-amber-600 dark:text-amber-400 mt-0.5 flex-shrink-0" />
            <div className="text-sm text-amber-800 dark:text-amber-200">
              <strong>Low match rate detected.</strong> The gene IDs may be from a different species or format.
              {useApiFallback && ' Enabling online lookup may improve results.'}
            </div>
          </div>
        )}

        <DialogFooter className="flex-col sm:flex-row gap-2">
          <Button
            variant="outline"
            onClick={handleSkip}
            disabled={isLoading}
            className="sm:order-1"
          >
            Skip Conversion
          </Button>
          <Button
            variant="ghost"
            onClick={onClose}
            disabled={isLoading}
            className="sm:order-2"
          >
            Cancel
          </Button>
          <Button
            onClick={handleConfirm}
            disabled={isLoading}
            className="sm:order-3"
          >
            {isLoading ? (
              <>
                <Loader2 className="h-4 w-4 mr-2 animate-spin" />
                Converting...
              </>
            ) : (
              <>
                <CheckCircle className="h-4 w-4 mr-2" />
                Convert {totalGenes.toLocaleString()} Genes
              </>
            )}
          </Button>
        </DialogFooter>
      </DialogContent>
    </Dialog>
  )
}
