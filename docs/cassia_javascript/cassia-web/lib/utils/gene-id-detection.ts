/**
 * Gene ID Detection Module
 *
 * Analyzes gene identifiers to determine their format:
 * - Ensembl IDs (human: ENSG..., mouse: ENSMUSG...)
 * - Entrez IDs (numeric)
 * - Gene symbols (already converted)
 */

export type GeneIdFormat =
  | 'ensembl_human'
  | 'ensembl_mouse'
  | 'entrez'
  | 'symbol'
  | 'mixed'
  | 'unknown'

export interface GeneIdDetectionResult {
  detectedFormat: GeneIdFormat
  confidence: number  // 0-1
  sampleSize: number
  matchCounts: {
    ensembl_human: number
    ensembl_mouse: number
    entrez: number
    symbol: number
    unknown: number
  }
  needsConversion: boolean
  suggestedSpecies: 'human' | 'mouse' | null
  sampleGenes: string[]  // First few genes for preview
}

// Gene ID patterns
const PATTERNS = {
  // Human Ensembl: ENSG followed by 11 digits (with optional version suffix)
  ENSEMBL_HUMAN: /^ENSG\d{11}(\.\d+)?$/,

  // Mouse Ensembl: ENSMUSG followed by 11 digits (with optional version suffix)
  ENSEMBL_MOUSE: /^ENSMUSG\d{11}(\.\d+)?$/,

  // Entrez IDs: pure numeric strings (typically 1-9 digits)
  ENTREZ: /^\d{1,9}$/,

  // Gene symbols: alphanumeric with possible hyphens, periods, underscores
  // Examples: TP53, CD3D, HLA-A, IL2, Cd3d (mouse), orf1ab
  SYMBOL: /^[A-Za-z][A-Za-z0-9\-_\.]*$/
}

/**
 * Classify a single gene identifier
 */
function classifyGeneId(gene: string): keyof typeof PATTERNS | 'unknown' {
  const trimmed = gene.trim()

  if (PATTERNS.ENSEMBL_HUMAN.test(trimmed)) {
    return 'ENSEMBL_HUMAN'
  }
  if (PATTERNS.ENSEMBL_MOUSE.test(trimmed)) {
    return 'ENSEMBL_MOUSE'
  }
  if (PATTERNS.ENTREZ.test(trimmed)) {
    return 'ENTREZ'
  }
  if (PATTERNS.SYMBOL.test(trimmed)) {
    return 'SYMBOL'
  }
  return 'unknown'
}

/**
 * Detect the gene ID format from a list of genes
 *
 * @param genes - Array of gene identifiers to analyze
 * @param sampleLimit - Maximum number of genes to sample (default: 200)
 * @returns Detection result with format, confidence, and conversion recommendation
 */
export function detectGeneIdFormat(
  genes: string[],
  sampleLimit: number = 200
): GeneIdDetectionResult {
  // Handle empty input
  if (!genes || genes.length === 0) {
    return {
      detectedFormat: 'unknown',
      confidence: 0,
      sampleSize: 0,
      matchCounts: {
        ensembl_human: 0,
        ensembl_mouse: 0,
        entrez: 0,
        symbol: 0,
        unknown: 0
      },
      needsConversion: false,
      suggestedSpecies: null,
      sampleGenes: []
    }
  }

  // Sample genes for analysis
  const sampleSize = Math.min(genes.length, sampleLimit)
  const sample = genes.slice(0, sampleSize)

  // Count matches for each pattern
  const counts = {
    ensembl_human: 0,
    ensembl_mouse: 0,
    entrez: 0,
    symbol: 0,
    unknown: 0
  }

  for (const gene of sample) {
    const classification = classifyGeneId(gene)
    switch (classification) {
      case 'ENSEMBL_HUMAN':
        counts.ensembl_human++
        break
      case 'ENSEMBL_MOUSE':
        counts.ensembl_mouse++
        break
      case 'ENTREZ':
        counts.entrez++
        break
      case 'SYMBOL':
        counts.symbol++
        break
      default:
        counts.unknown++
    }
  }

  // Determine dominant format
  const validTotal = sampleSize - counts.unknown
  const threshold = 0.8  // 80% threshold for confident detection

  let detectedFormat: GeneIdFormat = 'unknown'
  let confidence = 0
  let suggestedSpecies: 'human' | 'mouse' | null = null
  let needsConversion = false

  if (validTotal === 0) {
    detectedFormat = 'unknown'
    confidence = 0
  } else if (counts.ensembl_human / validTotal >= threshold) {
    detectedFormat = 'ensembl_human'
    confidence = counts.ensembl_human / validTotal
    suggestedSpecies = 'human'
    needsConversion = true
  } else if (counts.ensembl_mouse / validTotal >= threshold) {
    detectedFormat = 'ensembl_mouse'
    confidence = counts.ensembl_mouse / validTotal
    suggestedSpecies = 'mouse'
    needsConversion = true
  } else if (counts.entrez / validTotal >= threshold) {
    detectedFormat = 'entrez'
    confidence = counts.entrez / validTotal
    suggestedSpecies = null  // User needs to specify species for Entrez
    needsConversion = true
  } else if (counts.symbol / validTotal >= threshold) {
    detectedFormat = 'symbol'
    confidence = counts.symbol / validTotal
    needsConversion = false  // Already in symbol format
  } else {
    // Mixed format - no clear majority
    detectedFormat = 'mixed'
    const maxCount = Math.max(
      counts.ensembl_human,
      counts.ensembl_mouse,
      counts.entrez,
      counts.symbol
    )
    confidence = maxCount / validTotal

    // Still offer conversion if there's a significant portion of convertible IDs
    const convertibleCount = counts.ensembl_human + counts.ensembl_mouse + counts.entrez
    needsConversion = convertibleCount / validTotal > 0.3  // >30% convertible

    // Suggest species based on majority
    if (counts.ensembl_human > counts.ensembl_mouse) {
      suggestedSpecies = 'human'
    } else if (counts.ensembl_mouse > counts.ensembl_human) {
      suggestedSpecies = 'mouse'
    }
  }

  return {
    detectedFormat,
    confidence,
    sampleSize,
    matchCounts: counts,
    needsConversion,
    suggestedSpecies,
    sampleGenes: sample.slice(0, 10)  // First 10 for preview
  }
}

/**
 * Get human-readable label for a gene ID format
 */
export function getFormatLabel(format: GeneIdFormat): string {
  switch (format) {
    case 'ensembl_human':
      return 'Human Ensembl IDs'
    case 'ensembl_mouse':
      return 'Mouse Ensembl IDs'
    case 'entrez':
      return 'Entrez Gene IDs'
    case 'symbol':
      return 'Gene Symbols'
    case 'mixed':
      return 'Mixed Formats'
    case 'unknown':
    default:
      return 'Unknown Format'
  }
}

/**
 * Check if a single gene ID needs conversion
 */
export function geneNeedsConversion(gene: string): boolean {
  const classification = classifyGeneId(gene)
  return classification === 'ENSEMBL_HUMAN' ||
         classification === 'ENSEMBL_MOUSE' ||
         classification === 'ENTREZ'
}
