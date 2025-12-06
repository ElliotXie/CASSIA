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
 * Generates a standalone HTML file matching the Python version exactly
 */
export function exportSubclusteringReportHTML(
  data: SubclusterResult[],
  modelName?: string,
  title: string = 'Subclustering Annotation Report'
): string {
  // Input validation - use simple standalone HTML for empty data
  if (!data || !Array.isArray(data) || data.length === 0) {
    return `<!DOCTYPE html>
<html>
<head>
    <title>${escapeHtml(title)}</title>
    <style>
        body { font-family: 'Segoe UI', Arial, sans-serif; background: #f7f9fa; margin: 0; padding: 40px; }
        .container { max-width: 1200px; margin: 0 auto; background: #fff; border-radius: 16px; box-shadow: 0 4px 24px rgba(0,0,0,0.08); padding: 32px; text-align: center; color: #6b7280; }
    </style>
</head>
<body>
    <div class="container">
        <p>No subclustering data available to export.</p>
    </div>
</body>
</html>`
  }

  // Format timestamp to match Python's strftime('%Y-%m-%d %H:%M:%S')
  const now = new Date()
  const timestamp = `${now.getFullYear()}-${String(now.getMonth() + 1).padStart(2, '0')}-${String(now.getDate()).padStart(2, '0')} ${String(now.getHours()).padStart(2, '0')}:${String(now.getMinutes()).padStart(2, '0')}:${String(now.getSeconds()).padStart(2, '0')}`

  // Build table rows with popup markers (matching Python exactly)
  const rows = data.map(row => {
    const cluster = row.clusterId || ''
    const mainType = row.mainCellType || ''
    const subType = row.subCellType || ''
    const keyMarkers = row.keyMarkers || ''
    const reason = row.reason || ''

    // Escape quotes for JavaScript (matching Python's escaping)
    const escapedMarkers = String(keyMarkers).replace(/\\/g, '\\\\').replace(/'/g, "\\'").replace(/"/g, '&quot;').replace(/\n/g, '<br>')
    const escapedCluster = String(cluster).replace(/'/g, "\\'")

    return `
        <tr>
            <td class="cluster-col">${escapeHtml(String(cluster))}</td>
            <td class="annotation-col">
                <div class="main-type">${escapeHtml(mainType)}</div>
                <div class="sub-type">${escapeHtml(subType)}</div>
            </td>
            <td class="marker-col">
                <button class="marker-toggle" onclick="showMarkerPopup('${escapedCluster}', '${escapedMarkers}')">Show Markers</button>
            </td>
            <td class="reasoning-col"><div class="reasoning-box">${escapeHtml(reason)}</div></td>
        </tr>`
  }).join('')

  const displayModelName = modelName || title

  // Generate complete standalone HTML matching Python version exactly
  return `<!DOCTYPE html>
<html>
<head>
    <title>Subclustering Annotation Report - ${escapeHtml(displayModelName)}</title>
    <style>
        body { font-family: 'Segoe UI', Arial, sans-serif; background: #f7f9fa; margin: 0; }
        .container { max-width: 1200px; margin: 40px auto; background: #fff; border-radius: 16px; box-shadow: 0 4px 24px rgba(0,0,0,0.08); padding: 32px 32px 24px 32px; }
        h1 { color: #1976d2; margin-bottom: 8px; }
        .meta { color: #555; font-size: 15px; margin-bottom: 24px; }
        table { width: 100%; border-collapse: separate; border-spacing: 0; }
        th, td { padding: 12px 10px; text-align: left; }
        th { background: #1976d2; color: #fff; font-size: 16px; font-weight: 600; border-top-left-radius: 8px; border-top-right-radius: 8px; }
        tr { background: #fff; transition: background 0.2s; }
        tr:nth-child(even) { background: #f3f6fa; }
        tr:hover { background: #e3f2fd; }
        .cluster-col { width: 5%; font-weight: bold; color: #1976d2; font-size: 18px; }
        .annotation-col { width: 20%; }
        .main-type { font-weight: bold; color: #222; font-size: 16px; }
        .sub-type { color: #888; font-size: 14px; margin-top: 2px; }
        .marker-col { width: 20%; font-size: 13px; color: #888; }
        .marker-toggle { background: #e3f2fd; color: #1976d2; border: none; border-radius: 5px; padding: 3px 10px; font-size: 13px; cursor: pointer; }
        .marker-toggle:hover { background: #bbdefb; }
        .reasoning-col { width: 55%; }
        /* Modal popup styles */
        .modal-overlay { display: none; position: fixed; top: 0; left: 0; width: 100%; height: 100%; background: rgba(0,0,0,0.5); z-index: 1000; justify-content: center; align-items: center; }
        .modal-overlay.active { display: flex; }
        .modal-content { background: #fff; border-radius: 12px; padding: 24px; max-width: 600px; max-height: 80vh; overflow-y: auto; box-shadow: 0 8px 32px rgba(0,0,0,0.2); position: relative; }
        .modal-header { display: flex; justify-content: space-between; align-items: center; margin-bottom: 16px; border-bottom: 1px solid #e0e0e0; padding-bottom: 12px; }
        .modal-title { font-size: 18px; font-weight: 600; color: #1976d2; margin: 0; }
        .modal-close { background: #f44336; color: #fff; border: none; border-radius: 50%; width: 28px; height: 28px; font-size: 18px; cursor: pointer; display: flex; align-items: center; justify-content: center; }
        .modal-close:hover { background: #d32f2f; }
        .modal-body { font-size: 14px; color: #444; line-height: 1.6; }
        .reasoning-box { background: #fffde7; border-left: 5px solid #ffe082; border-radius: 7px; padding: 12px 16px; font-size: 15px; color: #444; box-shadow: 0 1px 4px rgba(255,193,7,0.07); }
    </style>
    <script>
    function showMarkerPopup(cluster, markers) {
        document.getElementById('modal-cluster').textContent = 'Cluster ' + cluster + ' - Key Markers';
        document.getElementById('modal-markers').innerHTML = markers;
        document.getElementById('marker-modal').classList.add('active');
    }
    function closeMarkerPopup() {
        document.getElementById('marker-modal').classList.remove('active');
    }
    // Close modal when clicking outside
    document.addEventListener('click', function(e) {
        var modal = document.getElementById('marker-modal');
        if (e.target === modal) {
            closeMarkerPopup();
        }
    });
    // Close modal with Escape key
    document.addEventListener('keydown', function(e) {
        if (e.key === 'Escape') {
            closeMarkerPopup();
        }
    });
    </script>
</head>
<body>
    <div class="container">
        <h1>Subclustering Annotation Report</h1>
        <div class="meta"><b>Model:</b> ${escapeHtml(displayModelName)} &nbsp; | &nbsp; <b>Generated:</b> ${timestamp}</div>
        <table>
            <tr>
                <th>Cluster</th>
                <th>Annotation<br><span style="font-weight:normal;font-size:12px">(Main / Subtype)</span></th>
                <th>Top Markers</th>
                <th>Reasoning</th>
            </tr>
            ${rows}
        </table>
    </div>
    <!-- Modal Popup -->
    <div id="marker-modal" class="modal-overlay">
        <div class="modal-content">
            <div class="modal-header">
                <h3 id="modal-cluster" class="modal-title">Key Markers</h3>
                <button class="modal-close" onclick="closeMarkerPopup()">&times;</button>
            </div>
            <div id="modal-markers" class="modal-body"></div>
        </div>
    </div>
</body>
</html>`
}

/**
 * Generate SVG pie chart for type distribution
 * Matches Python's matplotlib pie chart functionality
 */
function generateSVGPieChart(
  data: { name: string; value: number }[],
  title: string
): string {
  if (!data || data.length === 0) return ''

  const colors = ['#22c55e', '#3b82f6', '#f59e0b', '#ef4444', '#8b5cf6', '#ec4899', '#06b6d4', '#84cc16']
  const total = data.reduce((sum, d) => sum + d.value, 0)
  if (total === 0) return ''

  const size = 200
  const radius = 70
  const centerX = size / 2
  const centerY = size / 2

  let currentAngle = -Math.PI / 2 // Start from top

  const slices = data.map((d, i) => {
    const sliceAngle = (d.value / total) * 2 * Math.PI
    const startAngle = currentAngle
    const endAngle = currentAngle + sliceAngle

    const x1 = centerX + radius * Math.cos(startAngle)
    const y1 = centerY + radius * Math.sin(startAngle)
    const x2 = centerX + radius * Math.cos(endAngle)
    const y2 = centerY + radius * Math.sin(endAngle)

    const largeArcFlag = sliceAngle > Math.PI ? 1 : 0
    const pathData = `M ${centerX} ${centerY} L ${x1} ${y1} A ${radius} ${radius} 0 ${largeArcFlag} 1 ${x2} ${y2} Z`

    currentAngle = endAngle

    return `<path d="${pathData}" fill="${colors[i % colors.length]}" stroke="white" stroke-width="2"/>`
  }).join('')

  // Legend
  const legendItems = data.map((d, i) => {
    const percent = Math.round((d.value / total) * 100)
    return `
      <div style="display: flex; align-items: center; gap: 6px; margin-bottom: 4px;">
        <div style="width: 12px; height: 12px; border-radius: 2px; background: ${colors[i % colors.length]};"></div>
        <span style="font-size: 11px; color: #4b5563;">${escapeHtml(d.name)} (${percent}%)</span>
      </div>
    `
  }).join('')

  return `
    <div style="text-align: center;">
      <div style="font-weight: 600; font-size: 13px; color: #374151; margin-bottom: 8px;">${escapeHtml(title)}</div>
      <div style="display: flex; align-items: center; justify-content: center; gap: 16px;">
        <svg width="${size}" height="${size}" viewBox="0 0 ${size} ${size}">
          ${slices}
        </svg>
        <div style="text-align: left;">
          ${legendItems}
        </div>
      </div>
    </div>
  `
}

/**
 * Calculate type distribution from results
 */
function calculateDistribution(results: Array<{ mainType: string; subType: string }>, field: 'mainType' | 'subType'): { name: string; value: number }[] {
  const counts: Record<string, number> = {}
  results.forEach(r => {
    const value = r[field] || 'Unknown'
    counts[value] = (counts[value] || 0) + 1
  })
  return Object.entries(counts)
    .map(([name, value]) => ({ name, value }))
    .sort((a, b) => b.value - a.value)
}

/**
 * Export Uncertainty Report to HTML
 * Matches Python's generate_uq_batch_html_report() exactly
 */
export function exportUncertaintyReportHTML(
  data: ClusterUncertainty[],
  title: string = 'Batch Uncertainty Quantification Report',
  model?: string,
  provider?: string
): string {
  // Format timestamp to match Python's strftime('%Y-%m-%d %H:%M:%S')
  const now = new Date()
  const timestamp = `${now.getFullYear()}-${String(now.getMonth() + 1).padStart(2, '0')}-${String(now.getDate()).padStart(2, '0')} ${String(now.getHours()).padStart(2, '0')}:${String(now.getMinutes()).padStart(2, '0')}:${String(now.getSeconds()).padStart(2, '0')}`

  // Input validation
  if (!data || !Array.isArray(data) || data.length === 0) {
    return `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>${escapeHtml(title)}</title>
    <style>
        body { font-family: 'Segoe UI', Roboto, -apple-system, sans-serif; background: linear-gradient(135deg, #f5f7fa 0%, #e4e8ec 100%); min-height: 100vh; padding: 20px; }
        .container { max-width: 1200px; margin: 0 auto; background: white; border-radius: 16px; box-shadow: 0 4px 24px rgba(0,0,0,0.1); padding: 40px; text-align: center; color: #6b7280; }
    </style>
</head>
<body>
    <div class="container">
        <p>No uncertainty data available to export.</p>
    </div>
</body>
</html>`
  }

  // Calculate summary statistics
  const totalClusters = data.length
  const validScores = data.map(d => d.consensusScore).filter(s => typeof s === 'number' && !isNaN(s))
  const avgScore = validScores.length > 0 ? validScores.reduce((a, b) => a + b, 0) / validScores.length : 0

  // Get score class for average
  const getScoreClass = (score: number): string => {
    if (score >= 80) return 'score-high'
    if (score >= 50) return 'score-medium'
    return 'score-low'
  }

  const avgScoreClass = getScoreClass(avgScore)

  // Build cluster rows
  const clusterRows = data.map(cluster => {
    const mainType = cluster.consensusMainType || 'Unknown'
    const subType = cluster.consensusSubType || 'Unknown'
    const score = cluster.consensusScore || 0
    const scoreClass = getScoreClass(score)
    const scorePct = Math.round(score)

    // Get results (use originalResults with fallback)
    const originalResults = cluster.originalResults || cluster.roundResults || []
    const unifiedLlmResults = cluster.unifiedResultsLlm || []

    // Build original iterations table
    const originalIterRows = originalResults.map((r, i) =>
      `<tr><td>${r.iteration || i + 1}</td><td>${escapeHtml(r.mainType || 'N/A')}</td><td>${escapeHtml(r.subType || 'N/A')}</td></tr>`
    ).join('')
    const originalTable = originalIterRows || '<tr><td colspan="3">No iteration data</td></tr>'

    // Build unified LLM iterations table
    const unifiedLlmRows = unifiedLlmResults.map((r, i) =>
      `<tr><td>${r.iteration || i + 1}</td><td>${escapeHtml(r.mainType || 'N/A')}</td><td>${escapeHtml(r.subType || 'N/A')}</td></tr>`
    ).join('')
    const unifiedLlmTable = unifiedLlmRows || '<tr><td colspan="3">No data</td></tr>'

    // Generate pie charts from unified results or original results
    const resultsForCharts = unifiedLlmResults.length > 0 ? unifiedLlmResults : originalResults
    const mainTypeDist = calculateDistribution(resultsForCharts, 'mainType')
    const subTypeDist = calculateDistribution(resultsForCharts, 'subType')
    const mainPieChart = generateSVGPieChart(mainTypeDist, 'Main Cell Type Distribution')
    const subPieChart = generateSVGPieChart(subTypeDist, 'Sub Cell Type Distribution')

    // Safe ID for JavaScript
    const detailId = `cluster_detail_${cluster.clusterId.replace(/[^a-zA-Z0-9]/g, '_')}`

    // Score color for inline display
    const scoreColor = score >= 80 ? '#22c55e' : score >= 50 ? '#eab308' : '#ef4444'

    return `
            <tr class="cluster-row">
                <td style="font-weight: 600;">${escapeHtml(cluster.clusterId)}</td>
                <td>${escapeHtml(mainType)}</td>
                <td>${escapeHtml(subType)}</td>
                <td style="text-align: center;">
                    <span class="score-pill" style="background: ${scoreColor};">${scorePct}%</span>
                </td>
                <td style="text-align: center;">
                    <button class="toggle-btn" onclick="toggleClusterDetail('${detailId}')">Show</button>
                </td>
            </tr>
            <tr id="${detailId}" class="detail-row" style="display: none;">
                <td colspan="5">
                    <div style="display: flex; flex-wrap: wrap; gap: 20px; margin-bottom: 15px;">
                        ${mainPieChart}
                        ${subPieChart}
                    </div>
                    <div class="iteration-box">
                        <strong>Per-Iteration Results (Original):</strong>
                        <table class="iteration-table">
                            <thead><tr><th>Iter</th><th>Main Type</th><th>Sub Type</th></tr></thead>
                            <tbody>${originalTable}</tbody>
                        </table>
                    </div>
                    <div class="iteration-box" style="margin-top: 15px;">
                        <strong>Unified Per-Iteration Results (LLM):</strong>
                        <table class="iteration-table">
                            <thead><tr><th>Iter</th><th>Main Type</th><th>Sub Type</th></tr></thead>
                            <tbody>${unifiedLlmTable}</tbody>
                        </table>
                    </div>
                </td>
            </tr>`
  }).join('')

  // Generate complete standalone HTML matching Python version
  return `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>${escapeHtml(title)}</title>
    <style>
        * {
            box-sizing: border-box;
            margin: 0;
            padding: 0;
        }
        body {
            font-family: 'Segoe UI', Roboto, -apple-system, sans-serif;
            background: linear-gradient(135deg, #f5f7fa 0%, #e4e8ec 100%);
            min-height: 100vh;
            padding: 20px;
            line-height: 1.6;
        }
        .container {
            max-width: 1200px;
            margin: 0 auto;
            background: white;
            border-radius: 16px;
            box-shadow: 0 4px 24px rgba(0,0,0,0.1);
            overflow: hidden;
        }
        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            text-align: center;
        }
        .header h1 {
            font-size: 28px;
            margin-bottom: 10px;
            font-weight: 700;
        }
        .header-meta {
            font-size: 14px;
            opacity: 0.9;
        }
        .header-meta span {
            margin: 0 10px;
        }
        .section {
            padding: 25px 30px;
            border-bottom: 1px solid #e5e7eb;
        }
        .section:last-child {
            border-bottom: none;
        }
        .section-title {
            font-size: 18px;
            font-weight: 600;
            color: #1f2937;
            margin-bottom: 15px;
            display: flex;
            align-items: center;
            gap: 8px;
        }
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 15px;
        }
        .summary-card {
            background: linear-gradient(145deg, #f8fafc, #e2e8f0);
            border-radius: 12px;
            padding: 20px;
            text-align: center;
        }
        .summary-value {
            font-size: 32px;
            font-weight: 700;
            color: #1f2937;
        }
        .summary-label {
            font-size: 12px;
            color: #6b7280;
            margin-top: 5px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        .score-badge {
            width: 80px;
            height: 80px;
            border-radius: 50%;
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            margin: 0 auto 10px;
            color: white;
            font-weight: bold;
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
        }
        .score-high {
            background: linear-gradient(135deg, #22c55e, #16a34a);
        }
        .score-medium {
            background: linear-gradient(135deg, #eab308, #ca8a04);
        }
        .score-low {
            background: linear-gradient(135deg, #ef4444, #dc2626);
        }
        .score-pill {
            display: inline-block;
            padding: 4px 12px;
            border-radius: 12px;
            color: white;
            font-size: 13px;
            font-weight: 600;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 10px;
        }
        th {
            background: #f1f5f9;
            padding: 12px;
            text-align: left;
            font-weight: 600;
            color: #374151;
            font-size: 13px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        td {
            padding: 12px;
            border-bottom: 1px solid #e5e7eb;
            color: #4b5563;
        }
        .cluster-row:hover {
            background: #f8fafc;
        }
        .row-error {
            background-color: #fef2f2 !important;
        }
        .toggle-btn {
            background: #e0e7ff;
            color: #4338ca;
            border: none;
            padding: 6px 12px;
            border-radius: 6px;
            cursor: pointer;
            font-size: 12px;
            font-weight: 500;
            transition: all 0.2s;
        }
        .toggle-btn:hover {
            background: #c7d2fe;
        }
        .detail-row {
            background: #f8fafc;
        }
        .iteration-box {
            background: white;
            padding: 15px;
            border-radius: 8px;
            border: 1px solid #e5e7eb;
        }
        .iteration-table {
            margin-top: 10px;
            font-size: 13px;
        }
        .iteration-table th {
            background: #f8fafc;
            padding: 8px;
        }
        .iteration-table td {
            padding: 8px;
        }
        .footer {
            text-align: center;
            padding: 20px;
            color: #9ca3af;
            font-size: 12px;
        }
    </style>
    <script>
        function toggleClusterDetail(id) {
            var row = document.getElementById(id);
            var btn = row.previousElementSibling.querySelector('.toggle-btn');
            if (row.style.display === 'none') {
                row.style.display = 'table-row';
                btn.textContent = 'Hide';
            } else {
                row.style.display = 'none';
                btn.textContent = 'Show';
            }
        }
    </script>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>${escapeHtml(title)}</h1>
            <div class="header-meta">
                <span><strong>Model:</strong> ${escapeHtml(model || 'N/A')}</span>
                <span>|</span>
                <span><strong>Provider:</strong> ${escapeHtml(provider || 'N/A')}</span>
                <span>|</span>
                <span><strong>Total Clusters:</strong> ${totalClusters}</span>
            </div>
        </div>

        <div class="section">
            <h3 class="section-title">Summary Statistics</h3>
            <div class="summary-grid">
                <div class="summary-card">
                    <div class="summary-value">${totalClusters}</div>
                    <div class="summary-label">Total Clusters</div>
                </div>
                <div class="summary-card">
                    <div class="score-badge ${avgScoreClass}">
                        <span style="font-size: 24px;">${Math.round(avgScore)}%</span>
                    </div>
                    <div class="summary-label">Avg LLM Consensus Score</div>
                </div>
            </div>
        </div>

        <div class="section">
            <h3 class="section-title">Cluster Results</h3>
            <table>
                <thead>
                    <tr>
                        <th>Cluster ID</th>
                        <th>Main Cell Type</th>
                        <th>Sub Cell Type</th>
                        <th style="text-align: center;">LLM Consensus Score</th>
                        <th style="text-align: center;">Details</th>
                    </tr>
                </thead>
                <tbody>
                    ${clusterRows}
                </tbody>
            </table>
        </div>

        <div class="footer">
            Generated by CASSIA on ${timestamp}
        </div>
    </div>
</body>
</html>`
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
