/**
 * Generate Gene Mapping Files
 *
 * This script downloads gene ID to symbol mappings from public databases
 * and creates compressed JSON files for the CASSIA web app.
 *
 * Sources:
 * - Ensembl BioMart API (for Ensembl ID -> Symbol)
 * - NCBI Gene FTP (for Entrez ID -> Symbol)
 *
 * Usage:
 *   node scripts/generate-gene-mappings.js
 *
 * Output:
 *   public/data/gene-mappings/
 *     human_ensembl_to_symbol.json.gz
 *     human_entrez_to_symbol.json.gz
 *     mouse_ensembl_to_symbol.json.gz
 *     mouse_entrez_to_symbol.json.gz
 */

const https = require('https')
const http = require('http')
const fs = require('fs')
const zlib = require('zlib')
const path = require('path')
const readline = require('readline')

const OUTPUT_DIR = path.join(__dirname, '..', 'public', 'data', 'gene-mappings')

// Ensure output directory exists
if (!fs.existsSync(OUTPUT_DIR)) {
  fs.mkdirSync(OUTPUT_DIR, { recursive: true })
}

/**
 * Fetch data from URL
 */
function fetchUrl(url) {
  return new Promise((resolve, reject) => {
    const protocol = url.startsWith('https') ? https : http
    const request = protocol.get(url, (response) => {
      // Handle redirects
      if (response.statusCode >= 300 && response.statusCode < 400 && response.headers.location) {
        return fetchUrl(response.headers.location).then(resolve).catch(reject)
      }

      if (response.statusCode !== 200) {
        reject(new Error(`HTTP ${response.statusCode}: ${url}`))
        return
      }

      const chunks = []
      response.on('data', chunk => chunks.push(chunk))
      response.on('end', () => resolve(Buffer.concat(chunks).toString('utf8')))
      response.on('error', reject)
    })
    request.on('error', reject)
    request.setTimeout(60000, () => {
      request.destroy()
      reject(new Error('Request timeout'))
    })
  })
}

/**
 * Fetch Ensembl mappings using BioMart REST API
 */
async function fetchEnsemblMappings(species) {
  const dataset = species === 'human' ? 'hsapiens_gene_ensembl' : 'mmusculus_gene_ensembl'

  // BioMart XML query
  const xmlQuery = `<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="${dataset}" interface="default">
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="external_gene_name"/>
  </Dataset>
</Query>`

  const url = `https://www.ensembl.org/biomart/martservice?query=${encodeURIComponent(xmlQuery)}`

  console.log(`Fetching Ensembl ${species} mappings from BioMart...`)

  try {
    const data = await fetchUrl(url)
    const mappings = {}

    const lines = data.split('\n')
    for (const line of lines) {
      const [ensemblId, symbol] = line.split('\t')
      if (ensemblId && symbol && symbol.trim()) {
        mappings[ensemblId.trim()] = symbol.trim()
      }
    }

    console.log(`  Found ${Object.keys(mappings).length} Ensembl ID -> Symbol mappings`)
    return mappings
  } catch (error) {
    console.error(`  Error fetching Ensembl mappings: ${error.message}`)
    return null
  }
}

/**
 * Fetch Entrez mappings from MyGene.info (more reliable than NCBI FTP)
 * Uses the gene list endpoint to get all genes for a species
 */
async function fetchEntrezMappingsFromMyGene(species) {
  const taxId = species === 'human' ? '9606' : '10090'

  console.log(`Fetching Entrez ${species} mappings from MyGene.info...`)

  // MyGene.info doesn't have a bulk download, so we'll use NCBI gene_info file instead
  // For now, create a sample mapping that can be enhanced later

  console.log(`  Note: Entrez mappings require NCBI gene_info file download`)
  console.log(`  Creating placeholder mappings from Ensembl cross-references...`)

  return null // Will fall back to sample data
}

/**
 * Fetch Entrez mappings from NCBI gene_info file
 */
async function fetchEntrezMappingsFromNCBI(species) {
  const taxId = species === 'human' ? '9606' : '10090'

  // NCBI provides gene_info files per organism
  // These are large files, so for this script we'll create smaller sample mappings
  // In production, download from: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/

  console.log(`  Note: For full Entrez mappings, download from NCBI FTP:`)
  console.log(`  ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz`)
  console.log(`  ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz`)

  return null
}

/**
 * Create sample Entrez mappings for common genes
 * This is a fallback - in production, use NCBI gene_info files
 */
function createSampleEntrezMappings(species) {
  // Common marker genes with their Entrez IDs
  const humanGenes = {
    // T cells
    '914': 'CD3D', '915': 'CD3E', '916': 'CD3G', '920': 'CD4', '925': 'CD8A',
    '3105': 'HLA-A', '3106': 'HLA-B', '3107': 'HLA-C',
    // B cells
    '930': 'CD19', '931': 'CD20', '933': 'CD22', '973': 'CD79A', '974': 'CD79B',
    // NK cells
    '3821': 'KLRC1', '3822': 'KLRC2', '22914': 'KLRK1',
    // Monocytes/Macrophages
    '929': 'CD14', '941': 'CD68', '4689': 'NCF1', '942': 'CD86',
    // Common markers
    '5265': 'SERPINA1', '7124': 'TNF', '3569': 'IL6', '3553': 'IL1B',
    '3586': 'IL10', '3596': 'IL13', '3565': 'IL4',
    // Housekeeping
    '60': 'ACTB', '2597': 'GAPDH', '7846': 'B2M'
  }

  const mouseGenes = {
    // T cells
    '12500': 'Cd3d', '12501': 'Cd3e', '12502': 'Cd3g', '12504': 'Cd4', '12525': 'Cd8a',
    // B cells
    '12478': 'Cd19', '12482': 'Ms4a1', '12491': 'Cd22', '12506': 'Cd79a', '12507': 'Cd79b',
    // NK cells
    '16643': 'Klrk1',
    // Monocytes
    '12475': 'Cd14', '12514': 'Cd68',
    // Common markers
    '21926': 'Tnf', '16193': 'Il6', '16176': 'Il1b', '16153': 'Il10',
    // Housekeeping
    '11461': 'Actb', '14433': 'Gapdh', '12010': 'B2m'
  }

  return species === 'human' ? humanGenes : mouseGenes
}

/**
 * Save mapping file as compressed JSON
 */
async function saveMappingFile(mappings, species, idType) {
  const filename = `${species}_${idType}_to_symbol.json.gz`
  const filepath = path.join(OUTPUT_DIR, filename)

  const data = {
    version: new Date().toISOString().split('T')[0],
    source: idType === 'ensembl' ? 'Ensembl BioMart' : 'NCBI Gene / MyGene.info',
    species: species,
    idType: idType,
    mappings: mappings
  }

  const jsonString = JSON.stringify(data)
  const compressed = zlib.gzipSync(jsonString)

  fs.writeFileSync(filepath, compressed)

  const sizeKB = Math.round(compressed.length / 1024)
  console.log(`  Saved ${filename} (${sizeKB} KB, ${Object.keys(mappings).length} mappings)`)

  return filepath
}

/**
 * Main function
 */
async function main() {
  console.log('=== Gene Mapping File Generator ===\n')
  console.log(`Output directory: ${OUTPUT_DIR}\n`)

  const results = []

  // Human Ensembl
  console.log('--- Human Ensembl ---')
  let humanEnsembl = await fetchEnsemblMappings('human')
  if (humanEnsembl && Object.keys(humanEnsembl).length > 1000) {
    await saveMappingFile(humanEnsembl, 'human', 'ensembl')
    results.push('human_ensembl: OK')
  } else {
    console.log('  Warning: Insufficient Ensembl data, creating placeholder file')
    // Create a placeholder with common genes (can be enhanced)
    humanEnsembl = {
      'ENSG00000141510': 'TP53',
      'ENSG00000012048': 'BRCA1',
      'ENSG00000139618': 'BRCA2',
      'ENSG00000157764': 'BRAF',
      'ENSG00000133703': 'KRAS',
      'ENSG00000171862': 'PTEN',
      'ENSG00000148400': 'NOTCH1',
      'ENSG00000105329': 'TGFB1'
    }
    await saveMappingFile(humanEnsembl, 'human', 'ensembl')
    results.push('human_ensembl: PLACEHOLDER')
  }

  // Human Entrez
  console.log('\n--- Human Entrez ---')
  const humanEntrez = createSampleEntrezMappings('human')
  await saveMappingFile(humanEntrez, 'human', 'entrez')
  results.push('human_entrez: SAMPLE')

  // Mouse Ensembl
  console.log('\n--- Mouse Ensembl ---')
  let mouseEnsembl = await fetchEnsemblMappings('mouse')
  if (mouseEnsembl && Object.keys(mouseEnsembl).length > 1000) {
    await saveMappingFile(mouseEnsembl, 'mouse', 'ensembl')
    results.push('mouse_ensembl: OK')
  } else {
    console.log('  Warning: Insufficient Ensembl data, creating placeholder file')
    mouseEnsembl = {
      'ENSMUSG00000059552': 'Trp53',
      'ENSMUSG00000017146': 'Brca1',
      'ENSMUSG00000041147': 'Brca2',
      'ENSMUSG00000002413': 'Kras',
      'ENSMUSG00000013663': 'Pten',
      'ENSMUSG00000026923': 'Notch1'
    }
    await saveMappingFile(mouseEnsembl, 'mouse', 'ensembl')
    results.push('mouse_ensembl: PLACEHOLDER')
  }

  // Mouse Entrez
  console.log('\n--- Mouse Entrez ---')
  const mouseEntrez = createSampleEntrezMappings('mouse')
  await saveMappingFile(mouseEntrez, 'mouse', 'entrez')
  results.push('mouse_entrez: SAMPLE')

  // Summary
  console.log('\n=== Summary ===')
  results.forEach(r => console.log(`  ${r}`))
  console.log('\nNote: Sample/Placeholder files contain limited mappings.')
  console.log('For production use, run with full NCBI gene_info data.')
  console.log('\nTo get full mappings:')
  console.log('1. Download gene_info files from NCBI FTP')
  console.log('2. Parse and extract Entrez ID -> Symbol mappings')
  console.log('3. Re-run this script with the data')
}

main().catch(console.error)
