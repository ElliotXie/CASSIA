/**
 * Gene Mapping Loader
 *
 * Handles lazy loading and caching of gene ID mapping files.
 * Files are gzip compressed and decompressed in the browser using pako.
 */

import pako from 'pako'

export interface MappingFile {
  version: string
  source: string
  species: string
  idType: string
  mappings: Record<string, string>
}

export type Species = 'human' | 'mouse'
export type IdType = 'ensembl' | 'entrez'

// In-memory cache for loaded mappings
const mappingCache = new Map<string, Record<string, string>>()

// Track loading promises to prevent duplicate fetches
const loadingPromises = new Map<string, Promise<Record<string, string>>>()

/**
 * Generate cache key for a mapping file
 */
function getCacheKey(species: Species, idType: IdType): string {
  return `${species}_${idType}`
}

/**
 * Get the URL for a mapping file
 */
function getMappingUrl(species: Species, idType: IdType): string {
  return `/data/gene-mappings/${species}_${idType}_to_symbol.json.gz`
}

/**
 * Load and decompress a mapping file
 *
 * @param species - 'human' or 'mouse'
 * @param idType - 'ensembl' or 'entrez'
 * @returns Record mapping gene IDs to symbols
 */
export async function loadMappings(
  species: Species,
  idType: IdType
): Promise<Record<string, string>> {
  const cacheKey = getCacheKey(species, idType)

  // Return from cache if available
  if (mappingCache.has(cacheKey)) {
    return mappingCache.get(cacheKey)!
  }

  // Wait for existing load if in progress
  if (loadingPromises.has(cacheKey)) {
    return loadingPromises.get(cacheKey)!
  }

  // Start loading
  const loadPromise = (async (): Promise<Record<string, string>> => {
    try {
      const url = getMappingUrl(species, idType)
      console.log(`Loading gene mappings from: ${url}`)

      const response = await fetch(url)

      if (!response.ok) {
        throw new Error(`Failed to load mapping file: ${response.status} ${response.statusText}`)
      }

      const compressed = await response.arrayBuffer()
      console.log(`Downloaded ${compressed.byteLength} bytes (compressed)`)

      // Decompress using pako
      const decompressed = pako.ungzip(new Uint8Array(compressed), { to: 'string' })
      console.log(`Decompressed to ${decompressed.length} characters`)

      const data: MappingFile = JSON.parse(decompressed)

      if (!data.mappings || typeof data.mappings !== 'object') {
        throw new Error('Invalid mapping file format: missing mappings object')
      }

      const mappingCount = Object.keys(data.mappings).length
      console.log(`Loaded ${mappingCount} gene mappings for ${species} ${idType}`)

      // Cache the mappings
      mappingCache.set(cacheKey, data.mappings)

      return data.mappings
    } catch (error) {
      console.error(`Error loading mappings for ${species} ${idType}:`, error)
      throw error
    } finally {
      // Clean up loading promise
      loadingPromises.delete(cacheKey)
    }
  })()

  loadingPromises.set(cacheKey, loadPromise)
  return loadPromise
}

/**
 * Check if mappings are already cached
 */
export function isMappingCached(species: Species, idType: IdType): boolean {
  return mappingCache.has(getCacheKey(species, idType))
}

/**
 * Check if mappings are currently being loaded
 */
export function isMappingLoading(species: Species, idType: IdType): boolean {
  return loadingPromises.has(getCacheKey(species, idType))
}

/**
 * Clear the mapping cache (useful for testing or memory management)
 */
export function clearMappingCache(): void {
  mappingCache.clear()
}

/**
 * Get cache statistics
 */
export function getCacheStats(): {
  cachedCount: number
  loadingCount: number
  cachedMappings: string[]
} {
  return {
    cachedCount: mappingCache.size,
    loadingCount: loadingPromises.size,
    cachedMappings: Array.from(mappingCache.keys())
  }
}

/**
 * Preload mappings for a species (both Ensembl and Entrez)
 * Useful to start loading in background before user needs them
 */
export async function preloadMappingsForSpecies(species: Species): Promise<void> {
  await Promise.all([
    loadMappings(species, 'ensembl'),
    loadMappings(species, 'entrez')
  ])
}

/**
 * Try to load mappings with fallback to empty object
 * Useful when mappings are optional
 */
export async function loadMappingsSafe(
  species: Species,
  idType: IdType
): Promise<Record<string, string>> {
  try {
    return await loadMappings(species, idType)
  } catch (error) {
    console.warn(`Failed to load mappings for ${species} ${idType}, using empty mapping:`, error)
    return {}
  }
}
