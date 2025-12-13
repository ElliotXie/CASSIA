/**
 * Gene ID Converter
 *
 * Converts gene IDs (Ensembl, Entrez) to gene symbols.
 * Uses bundled mapping files with optional API fallback.
 */

import { loadMappings, loadMappingsSafe, type Species, type IdType } from './gene-mapping-loader'
import { queryMyGene } from './mygene-api'
import type { GeneIdDetectionResult } from './gene-id-detection'

export interface ConversionResult {
  original: string
  converted: string | null
  status: 'converted' | 'unchanged' | 'not_found' | 'api_fallback'
}

export interface BatchConversionResult {
  results: ConversionResult[]
  stats: {
    total: number
    converted: number
    unchanged: number
    notFound: number
    apiFallback: number
  }
  unmappedGenes: string[]
}

export interface ConversionOptions {
  species: Species
  idType: IdType
  useApiFallback: boolean
  keepOriginalOnFailure: boolean
  onProgress?: (current: number, total: number, stage: 'local' | 'api') => void
}

/**
 * Normalize Ensembl ID by stripping version suffix
 * ENSG00000000003.15 -> ENSG00000000003
 */
function normalizeEnsemblId(id: string): string {
  return id.split('.')[0]
}

/**
 * Determine the ID type from detection result
 */
export function getIdTypeFromDetection(detection: GeneIdDetectionResult): IdType {
  if (detection.detectedFormat === 'ensembl_human' || detection.detectedFormat === 'ensembl_mouse') {
    return 'ensembl'
  }
  return 'entrez'
}

/**
 * Determine the species from detection result
 */
export function getSpeciesFromDetection(detection: GeneIdDetectionResult): Species {
  if (detection.detectedFormat === 'ensembl_human') {
    return 'human'
  }
  if (detection.detectedFormat === 'ensembl_mouse') {
    return 'mouse'
  }
  // For Entrez or mixed, use suggested species or default to human
  return detection.suggestedSpecies || 'human'
}

/**
 * Convert a batch of gene IDs to gene symbols
 *
 * @param genes - Array of gene IDs to convert
 * @param options - Conversion options
 * @returns Batch conversion result with individual results and statistics
 */
export async function convertGeneIds(
  genes: string[],
  options: ConversionOptions
): Promise<BatchConversionResult> {
  const { species, idType, useApiFallback, keepOriginalOnFailure, onProgress } = options

  const results: ConversionResult[] = []
  const unmapped: string[] = []
  const stats = {
    total: genes.length,
    converted: 0,
    unchanged: 0,
    notFound: 0,
    apiFallback: 0
  }

  if (genes.length === 0) {
    return { results, stats, unmappedGenes: [] }
  }

  // Load mappings (with fallback to empty if loading fails)
  let mappings: Record<string, string>
  try {
    mappings = await loadMappings(species, idType)
  } catch (error) {
    console.warn('Failed to load local mappings, will rely on API:', error)
    mappings = {}
  }

  // First pass: local mappings
  for (let i = 0; i < genes.length; i++) {
    const gene = genes[i].trim()
    const lookupKey = idType === 'ensembl' ? normalizeEnsemblId(gene) : gene
    const symbol = mappings[lookupKey]

    if (symbol) {
      results.push({
        original: gene,
        converted: symbol,
        status: 'converted'
      })
      stats.converted++
    } else {
      unmapped.push(gene)
      results.push({
        original: gene,
        converted: null,
        status: 'not_found'
      })
    }

    // Progress callback every 100 genes
    if (onProgress && (i + 1) % 100 === 0) {
      onProgress(i + 1, genes.length, 'local')
    }
  }

  if (onProgress) {
    onProgress(genes.length, genes.length, 'local')
  }

  // Second pass: API fallback for unmapped genes
  if (useApiFallback && unmapped.length > 0) {
    console.log(`Querying MyGene API for ${unmapped.length} unmapped genes...`)

    try {
      const apiResults = await queryMyGene(
        unmapped,
        species,
        idType,
        (completed, total) => {
          if (onProgress) {
            onProgress(completed, total, 'api')
          }
        }
      )

      // Update results with API findings
      for (const result of results) {
        if (result.status === 'not_found' && apiResults.has(result.original)) {
          const apiSymbol = apiResults.get(result.original)
          if (apiSymbol) {
            result.converted = apiSymbol
            result.status = 'api_fallback'
            stats.notFound--
            stats.apiFallback++
          }
        }
      }
    } catch (error) {
      console.warn('API fallback failed:', error)
    }
  }

  // Apply keepOriginalOnFailure
  if (keepOriginalOnFailure) {
    for (const result of results) {
      if (result.status === 'not_found' && result.converted === null) {
        result.converted = result.original
        result.status = 'unchanged'
        stats.notFound--
        stats.unchanged++
      }
    }
  }

  // Final count of not found
  stats.notFound = results.filter(r => r.status === 'not_found').length

  return {
    results,
    stats,
    unmappedGenes: results.filter(r => r.status === 'not_found' || r.status === 'unchanged').map(r => r.original)
  }
}

/**
 * Generate a quick preview of conversions (first N genes only)
 * Uses local mappings only for speed
 */
export async function generateConversionPreview(
  genes: string[],
  detection: GeneIdDetectionResult,
  previewCount: number = 20
): Promise<ConversionResult[]> {
  const species = getSpeciesFromDetection(detection)
  const idType = getIdTypeFromDetection(detection)
  const previewGenes = genes.slice(0, previewCount)

  // Try to load mappings, use empty if fails
  const mappings = await loadMappingsSafe(species, idType)

  return previewGenes.map(gene => {
    const lookupKey = idType === 'ensembl' ? normalizeEnsemblId(gene) : gene
    const symbol = mappings[lookupKey]

    return {
      original: gene,
      converted: symbol || null,
      status: symbol ? 'converted' : 'not_found'
    } as ConversionResult
  })
}

/**
 * Apply conversion results to data array
 *
 * @param data - Original data array
 * @param geneColumn - Name of the gene column
 * @param conversionResults - Results from convertGeneIds
 * @returns New data array with converted gene IDs
 */
export function applyConversionToData(
  data: any[],
  geneColumn: string,
  conversionResults: ConversionResult[]
): any[] {
  // Create lookup map for fast conversion
  const conversionMap = new Map<string, string>()
  for (const result of conversionResults) {
    if (result.converted) {
      conversionMap.set(result.original, result.converted)
    }
  }

  // Apply conversions to data
  return data.map(row => {
    const originalGene = row[geneColumn]
    const convertedGene = conversionMap.get(originalGene)

    if (convertedGene && convertedGene !== originalGene) {
      return {
        ...row,
        [geneColumn]: convertedGene,
        _originalGeneId: originalGene  // Preserve original for reference
      }
    }
    return row
  })
}

/**
 * Create a summary message for conversion results
 */
export function getConversionSummary(stats: BatchConversionResult['stats']): string {
  const parts: string[] = []

  if (stats.converted > 0) {
    parts.push(`${stats.converted} converted via local mapping`)
  }
  if (stats.apiFallback > 0) {
    parts.push(`${stats.apiFallback} converted via API`)
  }
  if (stats.unchanged > 0) {
    parts.push(`${stats.unchanged} kept original (not found)`)
  }
  if (stats.notFound > 0) {
    parts.push(`${stats.notFound} could not be converted`)
  }

  return parts.join(', ') || 'No genes processed'
}
