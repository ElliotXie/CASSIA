import * as XLSX from 'xlsx'
import Papa from 'papaparse'
import { detectGeneIdFormat, type GeneIdDetectionResult } from './gene-id-detection'

export interface FileData {
  headers: string[]
  data: any[]
  rowCount: number
  clusterCount: number
  geneCount: number
  // Gene ID detection fields
  geneColumn?: string | null
  uniqueGenes?: string[]
  geneIdDetection?: GeneIdDetectionResult
}

export async function processFile(file: File): Promise<FileData> {
  try {
    if (!file) {
      throw new Error('No file provided')
    }

    if (file.size === 0) {
      throw new Error('File is empty')
    }

    if (file.size > 100 * 1024 * 1024) { // 100MB limit
      throw new Error('File is too large. Please upload a file smaller than 100MB.')
    }

    const extension = file.name.split('.').pop()?.toLowerCase()
    
    if (!extension) {
      throw new Error('Could not determine file extension')
    }
    
    if (extension === 'csv') {
      return await processCsvFile(file)
    } else if (extension === 'xlsx' || extension === 'xls') {
      return await processExcelFile(file)
    } else {
      throw new Error(`Unsupported file format: ${extension}. Please upload a CSV or XLSX file.`)
    }
  } catch (error) {
    console.error('File processing error:', error)
    throw error instanceof Error ? error : new Error('Unknown error processing file')
  }
}

async function processCsvFile(file: File): Promise<FileData> {
  return new Promise((resolve, reject) => {
    try {
      Papa.parse(file, {
        header: true,
        skipEmptyLines: 'greedy',
        transformHeader: (header: string) => header.trim(),
        complete: (results) => {
          try {
            console.log('CSV parsing results:', results)
            
            if (results.errors && results.errors.length > 0) {
              console.warn('CSV parsing warnings:', results.errors)
            }
            
            if (!results.data || !Array.isArray(results.data)) {
              reject(new Error('Invalid CSV format: No data found'))
              return
            }
            
            if (results.data.length === 0) {
              reject(new Error('CSV file is empty or contains no valid data'))
              return
            }
            
            // Filter out completely empty rows
            const data = results.data.filter((row: any) => {
              if (!row || typeof row !== 'object') return false
              return Object.values(row).some(val => 
                val !== null && val !== undefined && val !== ''
              )
            })
            
            if (data.length === 0) {
              reject(new Error('No valid data rows found in CSV file'))
              return
            }
            
            console.log(`Filtered data: ${data.length} rows`)
            
            const processed = analyzeData(data)
            resolve(processed)
          } catch (error) {
            console.error('Error in CSV processing:', error)
            reject(new Error(`Error processing CSV data: ${error instanceof Error ? error.message : 'Unknown error'}`))
          }
        },
        error: (error) => {
          console.error('Papa Parse error:', error)
          reject(new Error(`Error parsing CSV file: ${error.message || 'Unknown parsing error'}`))
        }
      })
    } catch (error) {
      reject(new Error(`Failed to start CSV parsing: ${error instanceof Error ? error.message : 'Unknown error'}`))
    }
  })
}

async function processExcelFile(file: File): Promise<FileData> {
  return new Promise((resolve, reject) => {
    try {
      const reader = new FileReader()
      
      reader.onload = (e) => {
        try {
          const data = e.target?.result as ArrayBuffer
          if (!data) {
            reject(new Error('Failed to read Excel file data'))
            return
          }
          
          const workbook = XLSX.read(data, { type: 'array' })
          
          if (!workbook.SheetNames || workbook.SheetNames.length === 0) {
            reject(new Error('Excel file contains no worksheets'))
            return
          }
          
          // Use the first sheet
          const sheetName = workbook.SheetNames[0]
          const worksheet = workbook.Sheets[sheetName]
          
          if (!worksheet) {
            reject(new Error(`Worksheet '${sheetName}' not found`))
            return
          }
          
          // Convert to JSON
          const jsonData = XLSX.utils.sheet_to_json(worksheet)
          
          if (!Array.isArray(jsonData)) {
            reject(new Error('Failed to parse Excel data as JSON'))
            return
          }
          
          console.log(`Excel parsing results: ${jsonData.length} rows`)
          
          const processed = analyzeData(jsonData)
          resolve(processed)
        } catch (error) {
          console.error('Excel processing error:', error)
          reject(new Error(`Error processing Excel file: ${error instanceof Error ? error.message : 'Unknown error'}`))
        }
      }
      
      reader.onerror = (error) => {
        console.error('FileReader error:', error)
        reject(new Error('Error reading Excel file'))
      }
      
      reader.readAsArrayBuffer(file)
    } catch (error) {
      reject(new Error(`Failed to start Excel file reading: ${error instanceof Error ? error.message : 'Unknown error'}`))
    }
  })
}

function analyzeData(data: any[]): FileData {
  try {
    if (!data || !Array.isArray(data) || data.length === 0) {
      throw new Error('File appears to be empty or contains no valid data')
    }
    
    // Find the first non-empty row to get headers
    const firstValidRow = data.find(row => row && typeof row === 'object' && Object.keys(row).length > 0)
    if (!firstValidRow) {
      throw new Error('No valid data rows found in file')
    }
    
    const headers = Object.keys(firstValidRow)
    if (headers.length === 0) {
      throw new Error('No column headers found in file')
    }
    
    // Detect important columns
    const clusterColumn = detectColumn(headers, ['cluster', 'seurat_clusters', 'true cell type'])
    const geneColumn = detectColumn(headers, ['gene', 'genes', 'gene_name', 'symbol'])
    
    // Count unique clusters and genes
    const clusters = new Set()
    const genes = new Set()
    
    data.forEach((row, index) => {
      try {
        if (!row || typeof row !== 'object') return
        
        if (clusterColumn && row[clusterColumn] !== undefined && row[clusterColumn] !== null && row[clusterColumn] !== '') {
          clusters.add(String(row[clusterColumn]).trim())
        }
        if (geneColumn && row[geneColumn] !== undefined && row[geneColumn] !== null && row[geneColumn] !== '') {
          genes.add(String(row[geneColumn]).trim())
        }
      } catch (err) {
        console.warn(`Error processing row ${index}:`, err)
      }
    })
    
    // Convert genes Set to array for detection
    const uniqueGenes = Array.from(genes) as string[]

    // Detect gene ID format
    const geneIdDetection = uniqueGenes.length > 0
      ? detectGeneIdFormat(uniqueGenes)
      : undefined

    if (geneIdDetection) {
      console.log('Gene ID detection result:', {
        format: geneIdDetection.detectedFormat,
        confidence: Math.round(geneIdDetection.confidence * 100) + '%',
        needsConversion: geneIdDetection.needsConversion,
        sampleGenes: geneIdDetection.sampleGenes.slice(0, 3)
      })
    }

    const result: FileData = {
      headers: headers || [],
      data: data || [],
      rowCount: data?.length || 0,
      clusterCount: clusters.size || 0,
      geneCount: genes.size || 0,
      // Gene ID detection fields
      geneColumn: geneColumn,
      uniqueGenes: uniqueGenes,
      geneIdDetection: geneIdDetection
    }

    console.log('Data analysis completed:', result)
    return result
    
  } catch (error) {
    console.error('Error in analyzeData:', error)
    throw new Error(`Failed to analyze file data: ${error instanceof Error ? error.message : 'Unknown error'}`)
  }
}

function detectColumn(headers: string[], possibleNames: string[]): string | null {
  const lowerHeaders = headers.map(h => h.toLowerCase())
  
  for (const name of possibleNames) {
    const index = lowerHeaders.indexOf(name.toLowerCase())
    if (index !== -1) {
      return headers[index]
    }
  }
  
  // Look for partial matches
  for (const name of possibleNames) {
    const match = lowerHeaders.find(h => h.includes(name.toLowerCase()))
    if (match) {
      return headers[lowerHeaders.indexOf(match)]
    }
  }
  
  return null
}

export function validateMarkerData(data: any[]): { isValid: boolean; errors: string[] } {
  const errors: string[] = []
  
  if (!data || data.length === 0) {
    errors.push('File is empty')
    return { isValid: false, errors }
  }
  
  const headers = Object.keys(data[0])
  
  // Check for required columns
  const hasCluster = detectColumn(headers, ['cluster', 'seurat_clusters', 'true cell type'])
  const hasGene = detectColumn(headers, ['gene', 'genes', 'gene_name', 'symbol'])
  
  if (!hasCluster) {
    errors.push('Could not find cluster column (expected: cluster, seurat_clusters, or true cell type)')
  }
  
  if (!hasGene) {
    errors.push('Could not find gene column (expected: gene, genes, gene_name, or symbol)')
  }
  
  // Check for statistical columns (recommended but not required)
  const hasStats = detectColumn(headers, ['p_val', 'pvalue', 'padj', 'p_val_adj', 'avg_log2fc', 'logfc'])
  if (!hasStats) {
    errors.push('Warning: No statistical columns found (p_val, avg_log2FC, etc.). Results may be less accurate.')
  }
  
  return {
    isValid: errors.filter(e => !e.startsWith('Warning:')).length === 0,
    errors
  }
}