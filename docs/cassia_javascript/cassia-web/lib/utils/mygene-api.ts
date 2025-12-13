/**
 * MyGene.info API Integration
 *
 * Fallback API for gene IDs not found in local mappings.
 * https://mygene.info/
 *
 * The API supports CORS and is free to use with reasonable rate limits.
 */

const MYGENE_API = 'https://mygene.info/v3/query'
const TIMEOUT_MS = 10000  // 10 second timeout per batch
const MAX_BATCH_SIZE = 1000  // API limit per request
const MAX_RETRIES = 1  // Retry once on failure

export type Species = 'human' | 'mouse'
export type IdType = 'ensembl' | 'entrez'

interface MyGeneResult {
  query: string
  symbol?: string
  notfound?: boolean
  _id?: string
  _score?: number
}

/**
 * Helper to chunk an array into smaller batches
 */
function chunkArray<T>(array: T[], size: number): T[][] {
  const chunks: T[][] = []
  for (let i = 0; i < array.length; i += size) {
    chunks.push(array.slice(i, i + size))
  }
  return chunks
}

/**
 * Fetch with timeout support
 */
async function fetchWithTimeout(
  url: string,
  options: RequestInit,
  timeout: number
): Promise<Response> {
  const controller = new AbortController()
  const timeoutId = setTimeout(() => controller.abort(), timeout)

  try {
    const response = await fetch(url, {
      ...options,
      signal: controller.signal
    })
    return response
  } finally {
    clearTimeout(timeoutId)
  }
}

/**
 * Get species ID for MyGene API
 */
function getSpeciesTaxId(species: Species): string {
  return species === 'human' ? '9606' : '10090'
}

/**
 * Get the appropriate field scope for the ID type
 */
function getScopes(idType: IdType): string {
  return idType === 'ensembl' ? 'ensembl.gene' : 'entrezgene'
}

/**
 * Query MyGene.info API for gene symbol mappings
 *
 * @param genes - Array of gene IDs to query
 * @param species - 'human' or 'mouse'
 * @param idType - 'ensembl' or 'entrez'
 * @param onProgress - Optional callback for progress updates
 * @returns Map of gene ID to symbol (null if not found)
 */
export async function queryMyGene(
  genes: string[],
  species: Species,
  idType: IdType,
  onProgress?: (completed: number, total: number) => void
): Promise<Map<string, string | null>> {
  const results = new Map<string, string | null>()

  if (genes.length === 0) {
    return results
  }

  const batches = chunkArray(genes, MAX_BATCH_SIZE)
  const speciesTaxId = getSpeciesTaxId(species)
  const scopes = getScopes(idType)

  let completedGenes = 0

  for (let batchIndex = 0; batchIndex < batches.length; batchIndex++) {
    const batch = batches[batchIndex]
    let retries = 0
    let success = false

    while (retries <= MAX_RETRIES && !success) {
      try {
        // Build query parameters
        const params = new URLSearchParams({
          q: batch.join(','),
          scopes: scopes,
          fields: 'symbol',
          species: speciesTaxId,
          size: String(batch.length)
        })

        const response = await fetchWithTimeout(
          `${MYGENE_API}?${params}`,
          {
            method: 'GET',
            headers: {
              'Accept': 'application/json'
            }
          },
          TIMEOUT_MS
        )

        if (!response.ok) {
          throw new Error(`API error: ${response.status} ${response.statusText}`)
        }

        const data = await response.json()

        // Handle both single result and array response
        const dataArray: MyGeneResult[] = Array.isArray(data) ? data : [data]

        // Process results
        for (const result of dataArray) {
          if (result.notfound) {
            results.set(result.query, null)
          } else if (result.symbol) {
            results.set(result.query, result.symbol)
          } else {
            // No symbol found
            results.set(result.query, null)
          }
        }

        // Mark any queries not in response as not found
        for (const gene of batch) {
          if (!results.has(gene)) {
            results.set(gene, null)
          }
        }

        success = true
        completedGenes += batch.length

        if (onProgress) {
          onProgress(completedGenes, genes.length)
        }

      } catch (error) {
        retries++

        if (retries > MAX_RETRIES) {
          console.warn(`MyGene API batch ${batchIndex + 1}/${batches.length} failed after ${MAX_RETRIES} retries:`, error)

          // Mark entire batch as not found on final failure
          for (const gene of batch) {
            results.set(gene, null)
          }

          completedGenes += batch.length
          if (onProgress) {
            onProgress(completedGenes, genes.length)
          }
        } else {
          console.log(`Retrying batch ${batchIndex + 1}/${batches.length} (attempt ${retries + 1})...`)
          // Brief delay before retry
          await new Promise(resolve => setTimeout(resolve, 500))
        }
      }
    }
  }

  return results
}

/**
 * Query a single gene ID (convenience wrapper)
 */
export async function queryMyGeneSingle(
  gene: string,
  species: Species,
  idType: IdType
): Promise<string | null> {
  const results = await queryMyGene([gene], species, idType)
  return results.get(gene) ?? null
}

/**
 * Check if MyGene API is reachable
 */
export async function isMyGeneAvailable(): Promise<boolean> {
  try {
    const response = await fetchWithTimeout(
      'https://mygene.info/v3/metadata',
      { method: 'GET' },
      5000
    )
    return response.ok
  } catch {
    return false
  }
}
