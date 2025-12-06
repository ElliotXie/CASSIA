/**
 * CASSIA Report HTML Export Utilities
 * Generates standalone HTML files for all report types
 */

import type { BatchReportData } from './BatchReport'
import type { ClusterHypothesis } from './HypothesisReport'
import type { EvaluationResult } from './EvaluationReport'
import type { SubclusterResult } from './SubclusteringReport'
import type { ClusterUncertainty } from './UncertaintyReport'

/**
 * Common HTML template wrapper with styling
 */
function htmlTemplate(title: string, content: string): string {
  return `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>${escapeHtml(title)}</title>
    <style>
        * { box-sizing: border-box; margin: 0; padding: 0; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            line-height: 1.6;
            color: #333;
            background: linear-gradient(135deg, #f0fdfa 0%, #ecfdf5 100%);
            min-height: 100vh;
            padding: 2rem;
        }
        .container { max-width: 1400px; margin: 0 auto; }
        h1 { color: #0d9488; font-size: 2rem; margin-bottom: 0.5rem; }
        h2 { color: #0f766e; font-size: 1.5rem; margin: 1.5rem 0 1rem; border-bottom: 2px solid #99f6e4; padding-bottom: 0.5rem; }
        h3 { color: #115e59; font-size: 1.2rem; margin: 1rem 0 0.5rem; }
        .header {
            background: linear-gradient(135deg, #0d9488 0%, #059669 100%);
            color: white;
            padding: 2rem;
            border-radius: 1rem;
            margin-bottom: 2rem;
            box-shadow: 0 4px 6px -1px rgba(0,0,0,0.1);
        }
        .header h1 { color: white; }
        .header p { color: #a7f3d0; margin-top: 0.5rem; }
        .meta { display: flex; gap: 1rem; margin-top: 1rem; flex-wrap: wrap; }
        .meta-item { background: rgba(255,255,255,0.2); padding: 0.5rem 1rem; border-radius: 0.5rem; font-size: 0.9rem; }
        .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1rem; margin: 1.5rem 0; }
        .stat-card {
            background: white;
            border: 1px solid #99f6e4;
            border-radius: 1rem;
            padding: 1.5rem;
            text-align: center;
            box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        }
        .stat-value { font-size: 2rem; font-weight: bold; color: #0d9488; }
        .stat-label { font-size: 0.85rem; color: #6b7280; text-transform: uppercase; letter-spacing: 0.05em; }
        .card {
            background: white;
            border: 1px solid #d1fae5;
            border-radius: 1rem;
            padding: 1.5rem;
            margin-bottom: 1.5rem;
            box-shadow: 0 2px 8px rgba(0,0,0,0.05);
            transition: box-shadow 0.3s;
        }
        .card:hover { box-shadow: 0 4px 12px rgba(0,0,0,0.1); }
        .card-header { display: flex; justify-content: space-between; align-items: start; margin-bottom: 1rem; }
        .card-title { font-weight: bold; font-size: 1.1rem; color: #1f2937; }
        .badge {
            display: inline-block;
            padding: 0.25rem 0.75rem;
            border-radius: 9999px;
            font-size: 0.75rem;
            font-weight: 600;
        }
        .badge-teal { background: linear-gradient(135deg, #14b8a6 0%, #10b981 100%); color: white; }
        .badge-amber { background: linear-gradient(135deg, #f59e0b 0%, #ea580c 100%); color: white; }
        .badge-green { background: #dcfce7; color: #166534; }
        .badge-yellow { background: #fef9c3; color: #854d0e; }
        .badge-red { background: #fee2e2; color: #991b1b; }
        .predicted-type {
            background: linear-gradient(135deg, #f0fdfa 0%, #d1fae5 100%);
            border-left: 4px solid #14b8a6;
            padding: 1rem;
            border-radius: 0.75rem;
            margin: 1rem 0;
        }
        .predicted-label { font-size: 0.75rem; color: #0d9488; text-transform: uppercase; font-weight: 600; }
        .predicted-value { font-size: 1.25rem; font-weight: bold; color: #115e59; margin-top: 0.25rem; }
        .sub-types { margin: 1rem 0; }
        .sub-type-item {
            padding: 0.5rem 1rem;
            border-radius: 0.5rem;
            margin-bottom: 0.5rem;
            font-size: 0.9rem;
        }
        .sub-type-1 { background: #ccfbf1; color: #0f766e; border-left: 3px solid #14b8a6; }
        .sub-type-2 { background: #d1fae5; color: #047857; border-left: 3px solid #10b981; }
        .sub-type-3 { background: #ecfdf5; color: #059669; border-left: 3px solid #34d399; }
        .metadata {
            background: #f9fafb;
            border-radius: 0.5rem;
            padding: 0.75rem;
            display: flex;
            flex-wrap: wrap;
            gap: 1rem;
            font-size: 0.85rem;
            color: #6b7280;
        }
        .metadata-item { display: flex; align-items: center; gap: 0.25rem; }
        .markers-preview { font-size: 0.85rem; color: #6b7280; margin-top: 0.75rem; }
        table { width: 100%; border-collapse: collapse; margin: 1rem 0; }
        th, td { padding: 0.75rem 1rem; text-align: left; border-bottom: 1px solid #e5e7eb; }
        th { background: #0d9488; color: white; font-weight: 600; }
        tr:nth-child(even) { background: #f9fafb; }
        tr:hover { background: #f0fdfa; }
        .section { margin: 2rem 0; }
        .collapsible { cursor: pointer; user-select: none; }
        .collapsible-content { display: none; padding: 1rem; background: #1f2937; color: #f3f4f6; border-radius: 0.5rem; margin-top: 0.5rem; white-space: pre-wrap; font-family: monospace; font-size: 0.85rem; overflow-x: auto; }
        .collapsible.active + .collapsible-content { display: block; }
        .grid-2 { display: grid; grid-template-columns: repeat(auto-fit, minmax(400px, 1fr)); gap: 1.5rem; }
        .footer { text-align: center; color: #9ca3af; font-size: 0.9rem; margin-top: 3rem; padding-top: 2rem; border-top: 1px solid #e5e7eb; }
        @media print {
            body { background: white; padding: 0; }
            .card { break-inside: avoid; box-shadow: none; border: 1px solid #ccc; }
            .header { background: #0d9488; -webkit-print-color-adjust: exact; print-color-adjust: exact; }
        }
    </style>
</head>
<body>
    <div class="container">
        ${content}
        <div class="footer">
            <p>Generated by CASSIA - ${new Date().toLocaleString()}</p>
        </div>
    </div>
    <script>
        document.querySelectorAll('.collapsible').forEach(el => {
            el.addEventListener('click', () => el.classList.toggle('active'));
        });
    </script>
</body>
</html>`
}

/**
 * Escape HTML special characters
 */
function escapeHtml(text: string): string {
  if (!text) return ''
  return text
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#039;')
}

/**
 * Get score badge class based on score value
 */
function getScoreBadgeClass(score: number | string | undefined): string {
  if (score === undefined || score === null) return 'badge'
  const numScore = typeof score === 'string' ? parseFloat(score) : score
  if (isNaN(numScore)) return 'badge'
  if (numScore >= 90) return 'badge badge-teal'
  if (numScore >= 75) return 'badge badge-green'
  if (numScore >= 60) return 'badge badge-yellow'
  return 'badge badge-red'
}

/**
 * Truncate text with ellipsis
 */
function truncate(text: string, maxLength: number = 100): string {
  if (!text) return ''
  if (text.length <= maxLength) return text
  return text.slice(0, maxLength) + '...'
}

/**
 * Export Batch Report to HTML
 */
export function exportBatchReportHTML(
  data: BatchReportData[],
  title: string = 'CASSIA Batch Analysis Report'
): string {
  // Input validation
  if (!data || !Array.isArray(data) || data.length === 0) {
    return htmlTemplate(title, '<div style="text-align: center; padding: 3rem; color: #6b7280;"><p>No data available to export.</p></div>')
  }

  const timestamp = new Date().toLocaleString()
  const tissues = [...new Set(data.map(d => d.tissue).filter(Boolean))]
  const models = [...new Set(data.map(d => d.model).filter(Boolean))]

  const statsHtml = `
    <div class="stats-grid">
      <div class="stat-card">
        <div class="stat-value">${data.length}</div>
        <div class="stat-label">Total Clusters</div>
      </div>
      <div class="stat-card">
        <div class="stat-value">${tissues.length}</div>
        <div class="stat-label">Tissues</div>
      </div>
      <div class="stat-card">
        <div class="stat-value">${models.length}</div>
        <div class="stat-label">Models</div>
      </div>
    </div>
  `

  const cardsHtml = data.map(row => {
    const subTypes = Array.isArray(row.subTypes)
      ? row.subTypes
      : typeof row.subTypes === 'string'
        ? row.subTypes.split(',').map(s => s.trim()).filter(Boolean)
        : []

    const subTypesHtml = subTypes.length > 0
      ? `<div class="sub-types">
          <div style="font-size: 0.75rem; color: #6b7280; text-transform: uppercase; margin-bottom: 0.5rem;">Sub Cell Types</div>
          ${subTypes.slice(0, 3).map((type, i) => `<div class="sub-type-item sub-type-${i + 1}">${escapeHtml(type)}</div>`).join('')}
          ${subTypes.length > 3 ? `<div style="color: #9ca3af; font-size: 0.85rem; padding-left: 0.5rem;">+${subTypes.length - 3} more</div>` : ''}
        </div>`
      : ''

    return `
      <div class="card">
        <div class="card-header">
          <div class="card-title">${escapeHtml(row.clusterId)}</div>
          <div style="display: flex; gap: 0.5rem;">
            <span class="badge badge-teal">${row.markerCount || 0} markers</span>
            ${row.iterations ? `<span class="badge badge-amber">${row.iterations} iter</span>` : ''}
            ${row.score !== undefined ? `<span class="${getScoreBadgeClass(row.score)}">${escapeHtml(String(row.score))}</span>` : ''}
          </div>
        </div>
        <div class="predicted-type">
          <div class="predicted-label">Predicted</div>
          <div class="predicted-value">${escapeHtml(row.mainType || 'Unknown')}</div>
        </div>
        ${subTypesHtml}
        <div class="metadata">
          <div class="metadata-item">Model: ${escapeHtml(row.model || 'N/A')}</div>
          <div class="metadata-item">Provider: ${escapeHtml(row.provider || 'N/A')}</div>
          <div class="metadata-item">Tissue: ${escapeHtml(row.tissue || 'N/A')}</div>
          <div class="metadata-item">Species: ${escapeHtml(row.species || 'N/A')}</div>
        </div>
        <div class="markers-preview">
          <strong>Markers:</strong> ${escapeHtml(truncate(row.markers || '', 150))}
        </div>
        ${row.conversationHistory ? `
          <div style="margin-top: 1rem;">
            <div class="collapsible" style="color: #0d9488; font-weight: 600; font-size: 0.9rem;">
              ▶ View Full Analysis
            </div>
            <div class="collapsible-content">${escapeHtml(row.conversationHistory)}</div>
          </div>
        ` : ''}
      </div>
    `
  }).join('')

  const content = `
    <div class="header">
      <h1>${escapeHtml(title)}</h1>
      <p>Comprehensive Cell Type Annotation Analysis</p>
      <div class="meta">
        <div class="meta-item">Generated: ${timestamp}</div>
        ${data[0]?.tissue ? `<div class="meta-item">Tissue: ${escapeHtml(data[0].tissue)}</div>` : ''}
        ${data[0]?.species ? `<div class="meta-item">Species: ${escapeHtml(data[0].species)}</div>` : ''}
      </div>
    </div>
    ${statsHtml}
    <div class="grid-2">
      ${cardsHtml}
    </div>
  `

  return htmlTemplate(title, content)
}

/**
 * Export Hypothesis Report to HTML
 */
export function exportHypothesisReportHTML(
  data: ClusterHypothesis[],
  title: string = 'CASSIA Hypothesis Generation Report'
): string {
  // Input validation
  if (!data || !Array.isArray(data) || data.length === 0) {
    return htmlTemplate(title, '<div style="text-align: center; padding: 3rem; color: #6b7280;"><p>No hypothesis data available to export.</p></div>')
  }

  const clustersHtml = data.map(cluster => {
    const hypothesesHtml = cluster.hypotheses.length > 0
      ? `
        <table>
          <thead>
            <tr>
              <th style="width: 80px;">Rank</th>
              <th style="width: 30%;">Predicted Cell Type</th>
              <th>Reasoning</th>
            </tr>
          </thead>
          <tbody>
            ${cluster.hypotheses.map((hyp, i) => `
              <tr>
                <td><span class="${i < 3 ? 'badge badge-' + ['teal', 'green', 'yellow'][i] : 'badge'}">${hyp.rank}</span></td>
                <td style="font-weight: 600;">${escapeHtml(hyp.cellType)}</td>
                <td>${escapeHtml(hyp.reasoning)}</td>
              </tr>
            `).join('')}
          </tbody>
        </table>
      `
      : '<p style="color: #9ca3af;">No hypotheses generated for this cluster.</p>'

    return `
      <div class="card">
        <h3 style="color: #0d9488; margin-bottom: 1rem;">Cluster: ${escapeHtml(cluster.clusterName)}</h3>
        ${cluster.error ? `<div class="badge badge-red">${escapeHtml(cluster.error)}</div>` : hypothesesHtml}
        ${cluster.rawResponse ? `
          <div style="margin-top: 1rem;">
            <div class="collapsible" style="color: #0d9488; font-weight: 600; font-size: 0.9rem;">
              ▶ View Raw LLM Output
            </div>
            <div class="collapsible-content">${escapeHtml(cluster.rawResponse)}</div>
          </div>
        ` : ''}
      </div>
    `
  }).join('')

  const content = `
    <div class="header">
      <h1>${escapeHtml(title)}</h1>
      <p>Ranked cell type hypotheses for each cluster</p>
      <div class="meta">
        <div class="meta-item">Generated: ${new Date().toLocaleString()}</div>
        <div class="meta-item">Clusters: ${data.length}</div>
      </div>
    </div>
    ${clustersHtml}
  `

  return htmlTemplate(title, content)
}

/**
 * Export Evaluation Report to HTML
 */
export function exportEvaluationReportHTML(
  data: EvaluationResult[],
  modelName?: string,
  title: string = 'LLM Cell Type Annotation Evaluation Report'
): string {
  // Input validation
  if (!data || !Array.isArray(data) || data.length === 0) {
    return htmlTemplate(title, '<div style="text-align: center; padding: 3rem; color: #6b7280;"><p>No evaluation data available to export.</p></div>')
  }

  const scores = data.map(d => d.score).filter(s => typeof s === 'number' && !isNaN(s))

  // Handle empty scores array
  if (scores.length === 0) {
    return htmlTemplate(title, '<div style="text-align: center; padding: 3rem; color: #6b7280;"><p>No valid scores available to export.</p></div>')
  }

  const mean = scores.reduce((a, b) => a + b, 0) / scores.length
  const sorted = [...scores].sort((a, b) => a - b)
  const median = scores.length % 2 === 0
    ? (sorted[scores.length / 2 - 1] + sorted[scores.length / 2]) / 2
    : sorted[Math.floor(scores.length / 2)]

  const isSimilarityScale = Math.max(...scores) > 10

  const statsHtml = `
    <div class="stats-grid">
      <div class="stat-card">
        <div class="stat-value">${mean.toFixed(2)}</div>
        <div class="stat-label">Mean Score</div>
      </div>
      <div class="stat-card">
        <div class="stat-value">${median.toFixed(2)}</div>
        <div class="stat-label">Median Score</div>
      </div>
      <div class="stat-card">
        <div class="stat-value">${Math.min(...scores).toFixed(2)}</div>
        <div class="stat-label">Min Score</div>
      </div>
      <div class="stat-card">
        <div class="stat-value">${Math.max(...scores).toFixed(2)}</div>
        <div class="stat-label">Max Score</div>
      </div>
    </div>
  `

  const sampleHtml = `
    <div class="section">
      <h2>Sample Results</h2>
      <table>
        <thead>
          <tr>
            <th>Gold Standard</th>
            <th>Prediction</th>
            <th>Score</th>
            <th>Explanation</th>
          </tr>
        </thead>
        <tbody>
          ${data.slice(0, 20).map(row => `
            <tr>
              <td style="font-weight: 500;">${escapeHtml(row.goldStandard)}</td>
              <td>${escapeHtml(row.prediction)}</td>
              <td><span class="${getScoreBadgeClass(row.score)}">${escapeHtml(String(row.score))}</span></td>
              <td style="max-width: 300px;">${escapeHtml(truncate(row.reasoning || '', 200))}</td>
            </tr>
          `).join('')}
        </tbody>
      </table>
      <p style="color: #9ca3af; font-size: 0.9rem; margin-top: 0.5rem;">Showing first 20 results.</p>
    </div>
  `

  const content = `
    <div class="header">
      <h1>${escapeHtml(title)}${modelName ? ` - ${escapeHtml(modelName)}` : ''}</h1>
      <p>Evaluation of LLM predictions against gold standard</p>
      <div class="meta">
        <div class="meta-item">Generated: ${new Date().toLocaleString()}</div>
        <div class="meta-item">Total Samples: ${data.length}</div>
        <div class="meta-item">Score Type: ${isSimilarityScale ? 'Similarity (0-100)' : 'Discrete (0-5)'}</div>
      </div>
    </div>
    ${statsHtml}
    ${sampleHtml}
  `

  return htmlTemplate(title, content)
}

/**
 * Export Subclustering Report to HTML
 */
export function exportSubclusteringReportHTML(
  data: SubclusterResult[],
  modelName?: string,
  title: string = 'Subclustering Annotation Report'
): string {
  // Input validation
  if (!data || !Array.isArray(data) || data.length === 0) {
    return htmlTemplate(title, '<div style="text-align: center; padding: 3rem; color: #6b7280;"><p>No subclustering data available to export.</p></div>')
  }

  const tableHtml = `
    <table>
      <thead>
        <tr>
          <th style="width: 100px;">Cluster</th>
          <th style="width: 25%;">Annotation</th>
          <th style="width: 20%;">Key Markers</th>
          <th>Reasoning</th>
        </tr>
      </thead>
      <tbody>
        ${data.map(row => `
          <tr>
            <td style="font-weight: bold; color: #0d9488; font-size: 1.1rem;">${escapeHtml(row.clusterId)}</td>
            <td>
              <div style="font-weight: 600; color: #1f2937;">${escapeHtml(row.mainCellType)}</div>
              ${row.subCellType ? `<div style="color: #6b7280; font-size: 0.9rem; margin-top: 0.25rem;">${escapeHtml(row.subCellType)}</div>` : ''}
            </td>
            <td style="font-size: 0.9rem; color: #6b7280;">${escapeHtml(truncate(row.keyMarkers || '', 100))}</td>
            <td>
              <div style="background: #f9fafb; border-radius: 0.5rem; padding: 0.75rem; font-size: 0.9rem; color: #4b5563;">
                ${escapeHtml(row.reason || 'No reasoning available')}
              </div>
            </td>
          </tr>
        `).join('')}
      </tbody>
    </table>
  `

  const uniqueMain = [...new Set(data.map(d => d.mainCellType))].length
  const uniqueSub = [...new Set(data.map(d => d.subCellType).filter(Boolean))].length

  const content = `
    <div class="header">
      <h1>${escapeHtml(title)}${modelName ? ` - ${escapeHtml(modelName)}` : ''}</h1>
      <p>Hierarchical annotation results</p>
      <div class="meta">
        <div class="meta-item">Generated: ${new Date().toLocaleString()}</div>
        <div class="meta-item">Total Subclusters: ${data.length}</div>
      </div>
    </div>
    <div class="stats-grid">
      <div class="stat-card">
        <div class="stat-value">${data.length}</div>
        <div class="stat-label">Total Subclusters</div>
      </div>
      <div class="stat-card">
        <div class="stat-value">${uniqueMain}</div>
        <div class="stat-label">Unique Main Types</div>
      </div>
      <div class="stat-card">
        <div class="stat-value">${uniqueSub}</div>
        <div class="stat-label">Unique Sub Types</div>
      </div>
    </div>
    ${tableHtml}
  `

  return htmlTemplate(title, content)
}

/**
 * Export Uncertainty Report to HTML
 */
export function exportUncertaintyReportHTML(
  data: ClusterUncertainty[],
  title: string = 'CASSIA Uncertainty Quantification Report'
): string {
  // Input validation
  if (!data || !Array.isArray(data) || data.length === 0) {
    return htmlTemplate(title, '<div style="text-align: center; padding: 3rem; color: #6b7280;"><p>No uncertainty data available to export.</p></div>')
  }

  const validScores = data.map(d => d.consensusScore).filter(s => typeof s === 'number' && !isNaN(s))
  const avgScore = validScores.length > 0 ? validScores.reduce((a, b) => a + b, 0) / validScores.length : 0
  const highConfidence = data.filter(d => d.consensusScore >= 80).length
  const lowConfidence = data.filter(d => d.consensusScore < 50).length

  const clustersHtml = data.map(cluster => {
    const scoreColor = cluster.consensusScore >= 80 ? '#22c55e' : cluster.consensusScore >= 50 ? '#eab308' : '#ef4444'

    const roundsHtml = cluster.roundResults.length > 0
      ? `
        <table style="margin-top: 1rem;">
          <thead>
            <tr>
              <th>Round</th>
              <th>Main Type</th>
              <th>Sub Type</th>
              <th>Agreement</th>
            </tr>
          </thead>
          <tbody>
            ${cluster.roundResults.map(round => {
              const mainMatch = round.mainType.toLowerCase() === cluster.consensusMainType.toLowerCase()
              const subMatch = round.subType.toLowerCase() === cluster.consensusSubType.toLowerCase()
              const agreement = mainMatch && subMatch ? 'Full' : mainMatch ? 'Partial' : 'None'
              const agreementClass = agreement === 'Full' ? 'badge-green' : agreement === 'Partial' ? 'badge-yellow' : 'badge-red'
              return `
                <tr>
                  <td>${round.iteration}</td>
                  <td style="${mainMatch ? 'color: #16a34a; font-weight: 500;' : ''}">${escapeHtml(round.mainType)}</td>
                  <td style="${subMatch ? 'color: #16a34a; font-weight: 500;' : ''}">${escapeHtml(round.subType)}</td>
                  <td><span class="badge ${agreementClass}">${agreement} Match</span></td>
                </tr>
              `
            }).join('')}
          </tbody>
        </table>
      `
      : ''

    return `
      <div class="card">
        <div class="card-header">
          <div class="card-title">${escapeHtml(cluster.clusterId)}</div>
          <div style="display: flex; align-items: center; gap: 0.5rem;">
            <span style="font-size: 1.5rem; font-weight: bold; color: ${scoreColor};">${cluster.consensusScore.toFixed(0)}%</span>
          </div>
        </div>
        <div class="predicted-type">
          <div class="predicted-label">Consensus</div>
          <div class="predicted-value">${escapeHtml(cluster.consensusMainType)}</div>
          <div style="color: #6b7280; margin-top: 0.25rem;">${escapeHtml(cluster.consensusSubType)}</div>
        </div>
        ${cluster.possibleMixedTypes && cluster.possibleMixedTypes.length > 0 ? `
          <div style="margin: 1rem 0;">
            <span style="font-size: 0.75rem; color: #6b7280; text-transform: uppercase;">Possible Mixed Types:</span>
            <div style="display: flex; flex-wrap: wrap; gap: 0.5rem; margin-top: 0.5rem;">
              ${cluster.possibleMixedTypes.map(t => `<span class="badge badge-yellow">${escapeHtml(t)}</span>`).join('')}
            </div>
          </div>
        ` : ''}
        ${cluster.llmReasoning ? `
          <div style="background: #f0fdfa; border-radius: 0.5rem; padding: 1rem; margin: 1rem 0; font-size: 0.9rem; color: #4b5563;">
            <strong style="color: #0d9488;">LLM Reasoning:</strong> ${escapeHtml(cluster.llmReasoning)}
          </div>
        ` : ''}
        ${roundsHtml}
      </div>
    `
  }).join('')

  const content = `
    <div class="header">
      <h1>${escapeHtml(title)}</h1>
      <p>Consensus scores and per-round analysis</p>
      <div class="meta">
        <div class="meta-item">Generated: ${new Date().toLocaleString()}</div>
        <div class="meta-item">Total Clusters: ${data.length}</div>
      </div>
    </div>
    <div class="stats-grid">
      <div class="stat-card">
        <div class="stat-value">${data.length}</div>
        <div class="stat-label">Total Clusters</div>
      </div>
      <div class="stat-card">
        <div class="stat-value">${avgScore.toFixed(1)}%</div>
        <div class="stat-label">Avg Score</div>
      </div>
      <div class="stat-card">
        <div class="stat-value" style="color: #22c55e;">${highConfidence}</div>
        <div class="stat-label">High Confidence</div>
      </div>
      <div class="stat-card">
        <div class="stat-value" style="color: #ef4444;">${lowConfidence}</div>
        <div class="stat-label">Low Confidence</div>
      </div>
    </div>
    ${clustersHtml}
  `

  return htmlTemplate(title, content)
}

/**
 * Download HTML string as file
 */
export function downloadHTML(html: string, filename: string): void {
  const blob = new Blob([html], { type: 'text/html;charset=utf-8' })
  const url = URL.createObjectURL(blob)
  const link = document.createElement('a')
  link.href = url
  link.download = filename
  document.body.appendChild(link)
  link.click()
  document.body.removeChild(link)
  URL.revokeObjectURL(url)
}
