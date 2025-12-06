/**
 * CASSIA Batch Analysis HTML Report Generator
 *
 * This module generates beautiful, interactive HTML reports from CASSIA batch analysis results.
 * Features: Modal popups for conversation history, search/filter, export options.
 * Theme: Elegant high-tech nature design with teal/emerald colors and glass-morphism effects.
 *
 * Ported from Python: generate_batch_report.py
 */

// ----------------- Helper Functions -----------------

/**
 * Escape HTML special characters to prevent XSS
 * @param {string} text - Text to escape
 * @returns {string} HTML-escaped text
 */
function escapeHtml(text) {
    if (text === null || text === undefined) return '';
    const str = String(text);
    const htmlEscapes = {
        '&': '&amp;',
        '<': '&lt;',
        '>': '&gt;',
        '"': '&quot;',
        "'": '&#39;'
    };
    return str.replace(/[&<>"']/g, char => htmlEscapes[char]);
}

/**
 * Truncate text with ellipsis if too long
 * @param {string} text - Text to truncate
 * @param {number} maxLength - Maximum length (default: 100)
 * @returns {string} Truncated text
 */
function truncateText(text, maxLength = 100) {
    if (!text) return '';
    if (text.length <= maxLength) return text;
    return text.substring(0, maxLength) + '...';
}

/**
 * Parse conversation history into structured sections
 * @param {string} historyText - Raw conversation history text with sections separated by " | "
 * @returns {Object} Dictionary with keys: 'annotation', 'validator', 'formatting', 'scoring'
 */
function parseConversationHistory(historyText) {
    const result = {
        annotation: '',
        validator: '',
        formatting: '',
        scoring: ''
    };

    if (!historyText) return result;

    const sections = historyText.split(' | ');

    for (const section of sections) {
        const trimmed = section.trim();
        if (trimmed.startsWith('Final Annotation Agent:')) {
            result.annotation = trimmed.replace('Final Annotation Agent:', '').trim();
        } else if (trimmed.startsWith('Coupling Validator:')) {
            result.validator = trimmed.replace('Coupling Validator:', '').trim();
        } else if (trimmed.startsWith('Formatting Agent:')) {
            result.formatting = trimmed.replace('Formatting Agent:', '').trim();
        } else if (trimmed.startsWith('Scoring Agent:')) {
            result.scoring = trimmed.replace('Scoring Agent:', '').trim();
        }
    }

    return result;
}

/**
 * Format analysis text with proper HTML structure
 * Converts markdown-style formatting to HTML and preserves line breaks
 * @param {string} text - Raw analysis text (may contain newlines)
 * @returns {string} HTML formatted text
 */
function formatAnalysisText(text) {
    if (!text) return '<p>No content available.</p>';

    // Escape HTML first
    let formatted = escapeHtml(text);

    // Convert **bold** to <strong>
    formatted = formatted.replace(/\*\*(.+?)\*\*/g, '<strong>$1</strong>');

    // Convert *italic* to <em>
    formatted = formatted.replace(/\*([^*]+)\*/g, '<em>$1</em>');

    // Convert ### headers
    formatted = formatted.replace(/^### (.+)$/gm, '<h4>$1</h4>');
    formatted = formatted.replace(/^## (.+)$/gm, '<h3>$1</h3>');
    formatted = formatted.replace(/^# (.+)$/gm, '<h2>$1</h2>');

    // Convert Step N. patterns to headers
    formatted = formatted.replace(/^(Step \d+[\.\:].+?)$/gm, '<h4 class="step-header">$1</h4>');

    // Convert bullet points (* item) at start of line
    formatted = formatted.replace(/^\* (.+)$/gm, '<li>$1</li>');

    // Convert numbered lists (1. item) at start of line
    formatted = formatted.replace(/^(\d+)\. (.+)$/gm, '<li>$2</li>');

    // Process line by line to handle structure
    const lines = formatted.split('\n');
    const resultLines = [];
    let inList = false;

    for (let line of lines) {
        line = line.trim();
        if (!line) {
            if (inList) {
                resultLines.push('</ul>');
                inList = false;
            }
            resultLines.push('<br>');
            continue;
        }

        if (line.startsWith('<li>')) {
            if (!inList) {
                resultLines.push('<ul>');
                inList = true;
            }
            resultLines.push(line);
        } else if (line.startsWith('<h')) {
            if (inList) {
                resultLines.push('</ul>');
                inList = false;
            }
            resultLines.push(line);
        } else {
            if (inList) {
                resultLines.push('</ul>');
                inList = false;
            }
            resultLines.push(`${line}<br>`);
        }
    }

    if (inList) {
        resultLines.push('</ul>');
    }

    return resultLines.join('\n');
}

/**
 * Parse and format JSON summary from Formatting Agent
 * @param {string} jsonText - Raw JSON text
 * @returns {string} HTML formatted summary
 */
function formatJsonSummary(jsonText) {
    try {
        // Find JSON content
        const jsonMatch = jsonText.match(/\{[\s\S]*\}/);
        if (jsonMatch) {
            const data = JSON.parse(jsonMatch[0]);

            const mainType = data.main_cell_type || 'Not specified';
            const subTypes = data.sub_cell_types || [];
            const mixedTypes = data.possible_mixed_cell_types || [];

            // Generate ranked subtype items with color coding
            const rankClasses = ['rank-1', 'rank-2', 'rank-3'];
            let subtypeItems;
            if (subTypes.length > 0) {
                subtypeItems = subTypes.map((t, i) =>
                    `<li class="${i < 3 ? rankClasses[i] : ''}">${escapeHtml(String(t))}</li>`
                ).join('');
            } else {
                subtypeItems = '<li class="empty">None identified</li>';
            }

            const mixedItems = mixedTypes.length > 0
                ? mixedTypes.map(t => `<li>${escapeHtml(String(t))}</li>`).join('')
                : '<li class="empty">None identified</li>';

            return `
            <div class="json-summary">
                <div class="summary-field">
                    <span class="field-label">Main Cell Type:</span>
                    <span class="field-value main-type">${escapeHtml(String(mainType))}</span>
                </div>
                <div class="summary-field">
                    <span class="field-label">Sub Cell Types:</span>
                    <ul class="type-list">
                        ${subtypeItems}
                    </ul>
                </div>
                <div class="summary-field">
                    <span class="field-label">Possible Mixed Types:</span>
                    <ul class="type-list">
                        ${mixedItems}
                    </ul>
                </div>
            </div>
            `;
        }
    } catch (e) {
        console.warn('Could not parse JSON in report:', e);
    }

    return `<pre class="json-raw">${escapeHtml(jsonText)}</pre>`;
}

// ----------------- Card Generation Functions -----------------

/**
 * Generate HTML for a single cluster card
 * @param {Object} row - Dictionary containing cluster data
 * @param {number} index - Cluster index for unique IDs
 * @returns {string} HTML string for the cluster card
 */
function generateClusterCard(row, index) {
    const trueType = escapeHtml(String(row['Cluster ID'] || 'Unknown'));
    const mainType = escapeHtml(String(row['Predicted General Cell Type'] || 'Unknown'));
    const subTypes = String(row['Predicted Detailed Cell Type'] || '');
    const mixedTypes = String(row['Possible Mixed Cell Types'] || '');
    const markerNum = String(row['Marker Number'] || 'N/A');
    const markerList = String(row['Marker List'] || '');
    const iterations = String(row['Iterations'] || '1');
    const model = escapeHtml(String(row['Model'] || 'N/A'));
    const provider = escapeHtml(String(row['Provider'] || 'N/A'));
    const tissue = escapeHtml(String(row['Tissue'] || 'N/A'));
    const species = escapeHtml(String(row['Species'] || 'N/A'));
    const additionalInfo = String(row['Additional Info'] || 'N/A');

    // Extract quality score if available
    const scoreValue = row['Score'];
    let scoreHtml = '';
    if (scoreValue !== null && scoreValue !== undefined && scoreValue !== '') {
        const scoreNum = parseFloat(scoreValue);
        if (!isNaN(scoreNum)) {
            let scoreClass;
            if (scoreNum > 90) {
                scoreClass = 'score-high';
            } else if (scoreNum >= 75) {
                scoreClass = 'score-medium';
            } else {
                scoreClass = 'score-low';
            }
            scoreHtml = `<span class="badge badge-score ${scoreClass}" title="Quality Score">Quality Score: ${Math.round(scoreNum)}</span>`;
        }
    }

    // Format sub types as list with ranked colors
    let subTypesHtml = '';
    if (subTypes && subTypes.trim()) {
        const subList = subTypes.split(',').map(s => s.trim());
        const rankClasses = ['rank-1', 'rank-2', 'rank-3'];
        subTypesHtml = '<ul class="sub-types-list">' +
            subList.slice(0, 3).map((s, i) =>
                `<li class="${rankClasses[i]}">${escapeHtml(s)}</li>`
            ).join('');
        if (subList.length > 3) {
            subTypesHtml += `<li class="more">+${subList.length - 3} more</li>`;
        }
        subTypesHtml += '</ul>';
    } else {
        subTypesHtml = '<span class="empty-field">None identified</span>';
    }

    // Format mixed types
    let mixedHtml = '';
    if (mixedTypes && mixedTypes.trim()) {
        mixedHtml = `
        <div class="card-field mixed-types">
            <span class="field-icon">
                <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <circle cx="12" cy="12" r="10"></circle>
                    <path d="M8 12h8"></path>
                    <path d="M12 8v8"></path>
                </svg>
            </span>
            <span class="field-label">Mixed Types:</span>
            <span class="field-value">${escapeHtml(truncateText(mixedTypes, 60))}</span>
        </div>
        `;
    }

    // Format additional info
    let additionalHtml = '';
    if (additionalInfo && additionalInfo.trim() && additionalInfo !== 'N/A' && additionalInfo.toLowerCase() !== 'nan') {
        additionalHtml = `
        <div class="card-field additional-info">
            <span class="field-icon">
                <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <circle cx="12" cy="12" r="10"></circle>
                    <path d="M12 16v-4"></path>
                    <path d="M12 8h.01"></path>
                </svg>
            </span>
            <span class="field-label">Additional Info:</span>
            <span class="field-value">${escapeHtml(truncateText(additionalInfo, 80))}</span>
        </div>
        `;
    }

    // Extract and format merged groupings if available
    const merged1 = String(row['Merged_Grouping_1'] || '');
    const merged2 = String(row['Merged_Grouping_2'] || '');
    const merged3 = String(row['Merged_Grouping_3'] || '');

    let mergedHtml = '';
    if (merged1.trim() || merged2.trim() || merged3.trim()) {
        const mergedItems = [];
        if (merged1.trim()) {
            mergedItems.push(`<div class="merged-item"><span class="merged-label">Broad:</span><span class="merged-tag merged-broad">${escapeHtml(truncateText(merged1, 40))}</span></div>`);
        }
        if (merged2.trim()) {
            mergedItems.push(`<div class="merged-item"><span class="merged-label">Detailed:</span><span class="merged-tag merged-detailed">${escapeHtml(truncateText(merged2, 40))}</span></div>`);
        }
        if (merged3.trim()) {
            mergedItems.push(`<div class="merged-item"><span class="merged-label">Specific:</span><span class="merged-tag merged-specific">${escapeHtml(truncateText(merged3, 40))}</span></div>`);
        }

        mergedHtml = `
        <div class="card-section merged-section">
            <div class="section-label">Merged Cell Type Groups (Broad -> Specific):</div>
            <div class="merged-groupings">
                ${mergedItems.join('')}
            </div>
        </div>
        `;
    }

    // Truncate marker list for preview
    const markerPreview = truncateText(markerList, 80);

    return `
    <div class="cluster-card" data-index="${index}" data-cluster="${trueType}"
         data-maintype="${mainType}" data-tissue="${tissue}" data-species="${species}"
         data-model="${model}" data-provider="${provider}">
        <div class="card-header">
            <h3 class="cluster-name">${trueType}</h3>
            <div class="card-badges">
                <span class="badge badge-markers" title="Number of markers">${markerNum} markers</span>
                <span class="badge badge-iterations" title="Iterations">${iterations} iter</span>
                ${scoreHtml}
            </div>
        </div>

        <div class="card-main-type">
            <span class="main-type-label">Predicted:</span>
            <span class="main-type-value">${mainType}</span>
        </div>

        <div class="card-section">
            <div class="section-label">Sub Cell Types:</div>
            ${subTypesHtml}
        </div>

        ${mergedHtml}

        ${mixedHtml}

        <div class="card-meta">
            <div class="meta-row">
                <span class="meta-item" title="Model">
                    <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                        <path d="M12 2L2 7l10 5 10-5-10-5z"></path>
                        <path d="M2 17l10 5 10-5"></path>
                        <path d="M2 12l10 5 10-5"></path>
                    </svg>
                    ${model}
                </span>
                <span class="meta-item" title="Provider">
                    <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                        <path d="M22 12h-4l-3 9L9 3l-3 9H2"></path>
                    </svg>
                    ${provider}
                </span>
            </div>
            <div class="meta-row">
                <span class="meta-item" title="Tissue">
                    <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                        <path d="M20.84 4.61a5.5 5.5 0 0 0-7.78 0L12 5.67l-1.06-1.06a5.5 5.5 0 0 0-7.78 7.78l1.06 1.06L12 21.23l7.78-7.78 1.06-1.06a5.5 5.5 0 0 0 0-7.78z"></path>
                    </svg>
                    ${tissue}
                </span>
                <span class="meta-item" title="Species">
                    <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                        <circle cx="12" cy="12" r="10"></circle>
                        <path d="M2 12h20"></path>
                        <path d="M12 2a15.3 15.3 0 0 1 4 10 15.3 15.3 0 0 1-4 10 15.3 15.3 0 0 1-4-10 15.3 15.3 0 0 1 4-10z"></path>
                    </svg>
                    ${species}
                </span>
            </div>
        </div>

        <div class="card-markers">
            <span class="markers-label">Markers:</span>
            <span class="markers-preview">${escapeHtml(markerPreview)}</span>
        </div>

        ${additionalHtml}

        <button class="view-analysis-btn" onclick="openModal(${index})">
            <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                <path d="M1 12s4-8 11-8 11 8 11 8-4 8-11 8-11-8-11-8z"></path>
                <circle cx="12" cy="12" r="3"></circle>
            </svg>
            View Full Analysis
        </button>
    </div>
    `;
}

/**
 * Generate HTML for modal content (conversation history)
 * @param {Object} row - Dictionary containing cluster data
 * @param {number} index - Cluster index
 * @returns {string} HTML string for modal content
 */
function generateModalContent(row, index) {
    const trueType = escapeHtml(String(row['Cluster ID'] || 'Unknown'));
    const mainType = escapeHtml(String(row['Predicted General Cell Type'] || 'Unknown'));
    const conversation = String(row['Conversation History'] || '');
    const markerList = escapeHtml(String(row['Marker List'] || ''));

    // Check for pre-formatted annotation HTML
    const preFormattedAnnotation = row['_formatted_annotation_html'] || '';

    // Parse conversation history
    const sections = parseConversationHistory(conversation);

    // Format each section
    let annotationHtml;
    if (preFormattedAnnotation) {
        annotationHtml = preFormattedAnnotation;
    } else {
        annotationHtml = sections.annotation
            ? formatAnalysisText(sections.annotation)
            : '<p>No annotation data available.</p>';
    }

    let validatorHtml = '';
    if (sections.validator) {
        const isPassed = sections.validator.toUpperCase().includes('VALIDATION PASSED');
        const statusClass = isPassed ? 'passed' : 'failed';
        const statusText = isPassed ? 'PASSED' : 'REVIEW NEEDED';
        validatorHtml = `
        <div class="validation-status ${statusClass}">
            <span class="status-icon">${isPassed ? '&#10004;' : '&#10008;'}</span>
            <span class="status-text">Validation ${statusText}</span>
        </div>
        <div class="validator-content">${escapeHtml(sections.validator)}</div>
        `;
    } else {
        validatorHtml = '<p>No validation data available.</p>';
    }

    const formattingHtml = sections.formatting
        ? formatJsonSummary(sections.formatting)
        : '<p>No formatting data available.</p>';

    // Scoring section if present
    let scoringHtml = '';
    if (sections.scoring) {
        scoringHtml = `
        <div class="modal-section scoring-section">
            <h3 class="section-title">
                <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <polygon points="12 2 15.09 8.26 22 9.27 17 14.14 18.18 21.02 12 17.77 5.82 21.02 7 14.14 2 9.27 8.91 8.26 12 2"></polygon>
                </svg>
                Quality Score
            </h3>
            <div class="section-content">
                ${formatAnalysisText(sections.scoring)}
            </div>
        </div>
        `;
    }

    // Merged groupings section if available
    const merged1 = String(row['Merged_Grouping_1'] || '');
    const merged2 = String(row['Merged_Grouping_2'] || '');
    const merged3 = String(row['Merged_Grouping_3'] || '');

    let mergedSectionHtml = '';
    if (merged1.trim() || merged2.trim() || merged3.trim()) {
        mergedSectionHtml = `
        <div class="modal-section merged-groupings-section">
            <h3 class="section-title">
                <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <path d="M17 21v-2a4 4 0 0 0-4-4H5a4 4 0 0 0-4 4v2"></path>
                    <circle cx="9" cy="7" r="4"></circle>
                    <path d="M23 21v-2a4 4 0 0 0-3-3.87"></path>
                    <path d="M16 3.13a4 4 0 0 1 0 7.75"></path>
                </svg>
                Merged Cell Type Groupings
            </h3>
            <div class="section-content">
                <div class="merged-detail-grid">
                    <div class="merged-detail-item">
                        <span class="merged-level">Broad:</span>
                        <span class="merged-value">${merged1.trim() ? escapeHtml(merged1) : 'N/A'}</span>
                    </div>
                    <div class="merged-detail-item">
                        <span class="merged-level">Detailed:</span>
                        <span class="merged-value">${merged2.trim() ? escapeHtml(merged2) : 'N/A'}</span>
                    </div>
                    <div class="merged-detail-item">
                        <span class="merged-level">Specific:</span>
                        <span class="merged-value">${merged3.trim() ? escapeHtml(merged3) : 'N/A'}</span>
                    </div>
                </div>
            </div>
        </div>
        `;
    }

    return `
    <div class="modal-content" id="modal-content-${index}">
        <div class="modal-header">
            <div class="modal-title-section">
                <h2 class="modal-title">${trueType}</h2>
                <p class="modal-subtitle">Predicted: <strong>${mainType}</strong></p>
            </div>
            <button class="modal-close" onclick="closeModal()">&times;</button>
        </div>

        <div class="modal-body">
            <div class="markers-full">
                <h4>Full Marker List:</h4>
                <p class="markers-text">${markerList}</p>
            </div>

            <div class="modal-section annotation-section">
                <h3 class="section-title">
                    <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                        <circle cx="11" cy="11" r="8"></circle>
                        <path d="M21 21l-4.35-4.35"></path>
                    </svg>
                    Annotation Analysis
                </h3>
                <div class="section-content">
                    ${annotationHtml}
                </div>
            </div>

            <div class="modal-section validator-section">
                <h3 class="section-title">
                    <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                        <path d="M22 11.08V12a10 10 0 1 1-5.93-9.14"></path>
                        <polyline points="22 4 12 14.01 9 11.01"></polyline>
                    </svg>
                    Validation Check
                </h3>
                <div class="section-content">
                    ${validatorHtml}
                </div>
            </div>

            <div class="modal-section formatting-section">
                <h3 class="section-title">
                    <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                        <path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"></path>
                        <polyline points="14 2 14 8 20 8"></polyline>
                        <line x1="16" y1="13" x2="8" y2="13"></line>
                        <line x1="16" y1="17" x2="8" y2="17"></line>
                        <polyline points="10 9 9 9 8 9"></polyline>
                    </svg>
                    Structured Summary
                </h3>
                <div class="section-content">
                    ${formattingHtml}
                </div>
            </div>

            ${scoringHtml}

            ${mergedSectionHtml}
        </div>
    </div>
    `;
}

// ----------------- CSS Styles -----------------

/**
 * Return the CSS styles for the report (exact copy from Python)
 * @returns {string} CSS styles
 */
function getCssStyles() {
    return `
    :root {
        --primary: #0d9488;
        --primary-light: #14b8a6;
        --primary-dark: #0f766e;
        --secondary: #059669;
        --accent: #f59e0b;
        --accent-light: #fbbf24;
        --background: linear-gradient(135deg, #f0fdf4 0%, #ecfeff 50%, #f0fdf4 100%);
        --surface: rgba(255, 255, 255, 0.85);
        --surface-hover: rgba(255, 255, 255, 0.95);
        --text-primary: #134e4a;
        --text-secondary: #64748b;
        --text-muted: #94a3b8;
        --border: rgba(13, 148, 136, 0.15);
        --shadow: 0 8px 32px rgba(13, 148, 136, 0.12);
        --shadow-hover: 0 12px 40px rgba(13, 148, 136, 0.2);
    }

    * {
        margin: 0;
        padding: 0;
        box-sizing: border-box;
    }

    body {
        font-family: 'Inter', 'Segoe UI', system-ui, -apple-system, sans-serif;
        background: var(--background);
        min-height: 100vh;
        color: var(--text-primary);
        line-height: 1.6;
    }

    .container {
        max-width: 1600px;
        margin: 0 auto;
        padding: 20px;
    }

    /* Header */
    .report-header {
        background: linear-gradient(135deg, var(--primary) 0%, var(--secondary) 100%);
        color: white;
        padding: 40px;
        border-radius: 24px;
        margin-bottom: 30px;
        box-shadow: var(--shadow);
        position: relative;
        overflow: hidden;
    }

    .report-header::before {
        content: '';
        position: absolute;
        top: -50%;
        right: -50%;
        width: 100%;
        height: 200%;
        background: radial-gradient(circle, rgba(255,255,255,0.1) 0%, transparent 60%);
        pointer-events: none;
    }

    .header-content {
        position: relative;
        z-index: 1;
    }

    .report-title {
        font-size: 2.5rem;
        font-weight: 800;
        margin-bottom: 8px;
        letter-spacing: -0.5px;
    }

    .report-subtitle {
        font-size: 1.1rem;
        opacity: 0.9;
        font-weight: 400;
    }

    .report-meta {
        display: flex;
        gap: 24px;
        margin-top: 20px;
        flex-wrap: wrap;
    }

    .meta-badge {
        display: flex;
        align-items: center;
        gap: 8px;
        background: rgba(255, 255, 255, 0.2);
        padding: 8px 16px;
        border-radius: 12px;
        font-size: 0.9rem;
        backdrop-filter: blur(10px);
    }

    /* Stats Bar */
    .stats-bar {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
        gap: 16px;
        margin-bottom: 30px;
    }

    .stat-card {
        background: var(--surface);
        backdrop-filter: blur(10px);
        border: 1px solid var(--border);
        border-radius: 16px;
        padding: 20px;
        text-align: center;
        box-shadow: var(--shadow);
    }

    .stat-value {
        font-size: 2rem;
        font-weight: 700;
        color: var(--primary);
        line-height: 1;
    }

    .stat-label {
        font-size: 0.85rem;
        color: var(--text-secondary);
        margin-top: 8px;
        text-transform: uppercase;
        letter-spacing: 0.5px;
    }

    /* Controls */
    .controls-bar {
        background: var(--surface);
        backdrop-filter: blur(10px);
        border: 1px solid var(--border);
        border-radius: 16px;
        padding: 20px;
        margin-bottom: 30px;
        display: flex;
        flex-wrap: wrap;
        gap: 16px;
        align-items: center;
        box-shadow: var(--shadow);
    }

    .search-box {
        flex: 1;
        min-width: 250px;
        position: relative;
    }

    .search-box input {
        width: 100%;
        padding: 12px 16px 12px 44px;
        border: 2px solid var(--border);
        border-radius: 12px;
        font-size: 1rem;
        background: white;
        color: var(--text-primary);
        transition: all 0.3s ease;
    }

    .search-box input:focus {
        outline: none;
        border-color: var(--primary);
        box-shadow: 0 0 0 4px rgba(13, 148, 136, 0.1);
    }

    .search-box svg {
        position: absolute;
        left: 14px;
        top: 50%;
        transform: translateY(-50%);
        color: var(--text-muted);
    }

    .filter-select {
        padding: 12px 16px;
        border: 2px solid var(--border);
        border-radius: 12px;
        font-size: 0.9rem;
        background: white;
        color: var(--text-primary);
        cursor: pointer;
        min-width: 140px;
        transition: all 0.3s ease;
    }

    .filter-select:focus {
        outline: none;
        border-color: var(--primary);
    }

    .btn {
        display: inline-flex;
        align-items: center;
        gap: 8px;
        padding: 12px 20px;
        border: none;
        border-radius: 12px;
        font-size: 0.9rem;
        font-weight: 600;
        cursor: pointer;
        transition: all 0.3s ease;
    }

    .btn-primary {
        background: linear-gradient(135deg, var(--primary) 0%, var(--secondary) 100%);
        color: white;
    }

    .btn-primary:hover {
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(13, 148, 136, 0.3);
    }

    .btn-secondary {
        background: white;
        color: var(--text-primary);
        border: 2px solid var(--border);
    }

    .btn-secondary:hover {
        border-color: var(--primary);
        color: var(--primary);
    }

    /* Cluster Grid */
    .cluster-grid {
        display: grid;
        grid-template-columns: repeat(auto-fill, minmax(360px, 1fr));
        gap: 24px;
    }

    .cluster-card {
        background: var(--surface);
        backdrop-filter: blur(10px);
        border: 1px solid var(--border);
        border-radius: 20px;
        padding: 24px;
        transition: all 0.3s ease;
        display: flex;
        flex-direction: column;
        gap: 16px;
    }

    .cluster-card:hover {
        transform: translateY(-4px);
        box-shadow: var(--shadow-hover);
        border-color: var(--primary-light);
        background: var(--surface-hover);
    }

    .card-header {
        display: flex;
        justify-content: space-between;
        align-items: flex-start;
        gap: 12px;
    }

    .cluster-name {
        font-size: 1.1rem;
        font-weight: 700;
        color: var(--text-primary);
        line-height: 1.3;
    }

    .card-badges {
        display: flex;
        gap: 6px;
        flex-wrap: wrap;
        justify-content: flex-end;
    }

    .badge {
        padding: 4px 10px;
        border-radius: 8px;
        font-size: 0.75rem;
        font-weight: 600;
        white-space: nowrap;
    }

    .badge-markers {
        background: linear-gradient(135deg, var(--primary-light), var(--primary));
        color: white;
    }

    .badge-iterations {
        background: linear-gradient(135deg, var(--accent-light), var(--accent));
        color: white;
    }

    .badge-score {
        font-weight: 600;
    }

    .score-high {
        background: linear-gradient(135deg, #10b981, #059669) !important;
        color: white;
    }

    .score-medium {
        background: linear-gradient(135deg, #86efac, #4ade80) !important;
        color: #166534;
    }

    .score-low {
        background: linear-gradient(135deg, #ef4444, #dc2626) !important;
        color: white;
    }

    .card-main-type {
        background: linear-gradient(135deg, rgba(13, 148, 136, 0.08), rgba(5, 150, 105, 0.08));
        border-radius: 12px;
        padding: 16px;
        border-left: 4px solid var(--primary);
    }

    .main-type-label {
        font-size: 0.8rem;
        color: var(--text-secondary);
        text-transform: uppercase;
        letter-spacing: 0.5px;
        display: block;
        margin-bottom: 4px;
    }

    .main-type-value {
        font-size: 1.15rem;
        font-weight: 700;
        color: var(--primary-dark);
    }

    .card-section {
        padding: 0;
    }

    .section-label {
        font-size: 0.8rem;
        color: var(--text-secondary);
        text-transform: uppercase;
        letter-spacing: 0.5px;
        margin-bottom: 8px;
    }

    .sub-types-list {
        list-style: none;
        display: flex;
        flex-direction: column;
        gap: 6px;
    }

    .sub-types-list li {
        font-size: 0.9rem;
        color: var(--text-primary);
        padding-left: 16px;
        position: relative;
    }

    .sub-types-list li::before {
        content: '';
        position: absolute;
        left: 0;
        top: 8px;
        width: 6px;
        height: 6px;
        background: var(--secondary);
        border-radius: 50%;
    }

    .sub-types-list li.more {
        color: var(--text-muted);
        font-style: italic;
    }

    .sub-types-list li.more::before {
        background: var(--text-muted);
    }

    /* Ranked subtype colors with BACKGROUND */
    .sub-types-list li.rank-1 {
        background: rgba(13, 148, 136, 0.25);
        color: #0d9488;
        font-weight: 600;
        padding: 8px 12px;
        border-radius: 6px;
        border-left: 3px solid #0d9488;
        margin-bottom: 4px;
    }
    .sub-types-list li.rank-1::before {
        display: none;
    }

    .sub-types-list li.rank-2 {
        background: rgba(20, 184, 166, 0.15);
        color: #134e4a;
        padding: 8px 12px;
        border-radius: 6px;
        border-left: 3px solid #14b8a6;
        margin-bottom: 4px;
    }
    .sub-types-list li.rank-2::before {
        display: none;
    }

    .sub-types-list li.rank-3 {
        background: rgba(94, 234, 212, 0.12);
        color: #64748b;
        padding: 8px 12px;
        border-radius: 6px;
        border-left: 3px solid #5eead4;
        margin-bottom: 4px;
    }
    .sub-types-list li.rank-3::before {
        display: none;
    }

    .empty-field {
        color: var(--text-muted);
        font-style: italic;
        font-size: 0.9rem;
    }

    .card-field {
        display: flex;
        align-items: flex-start;
        gap: 8px;
        font-size: 0.85rem;
        color: var(--text-secondary);
    }

    .card-field .field-icon {
        color: var(--primary);
        flex-shrink: 0;
        margin-top: 2px;
    }

    .card-field .field-label {
        font-weight: 600;
        white-space: nowrap;
    }

    .card-meta {
        display: flex;
        flex-direction: column;
        gap: 8px;
        padding: 12px;
        background: rgba(100, 116, 139, 0.05);
        border-radius: 10px;
    }

    .meta-row {
        display: flex;
        gap: 16px;
        flex-wrap: wrap;
    }

    .meta-item {
        display: flex;
        align-items: center;
        gap: 6px;
        font-size: 0.8rem;
        color: var(--text-secondary);
    }

    .meta-item svg {
        color: var(--primary);
    }

    .card-markers {
        font-size: 0.8rem;
        color: var(--text-secondary);
    }

    .markers-label {
        font-weight: 600;
    }

    .markers-preview {
        color: var(--text-muted);
        word-break: break-word;
    }

    .view-analysis-btn {
        display: flex;
        align-items: center;
        justify-content: center;
        gap: 8px;
        width: 100%;
        padding: 14px 20px;
        background: linear-gradient(135deg, var(--primary) 0%, var(--secondary) 100%);
        color: white;
        border: none;
        border-radius: 12px;
        font-size: 0.95rem;
        font-weight: 600;
        cursor: pointer;
        transition: all 0.3s ease;
        margin-top: auto;
    }

    .view-analysis-btn:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 20px rgba(13, 148, 136, 0.3);
    }

    /* Modal */
    .modal-overlay {
        display: none;
        position: fixed;
        top: 0;
        left: 0;
        width: 100%;
        height: 100%;
        background: rgba(15, 23, 42, 0.7);
        backdrop-filter: blur(8px);
        z-index: 1000;
        justify-content: center;
        align-items: center;
        padding: 20px;
        opacity: 0;
        transition: opacity 0.3s ease;
    }

    .modal-overlay.active {
        display: flex;
        opacity: 1;
    }

    .modal-container {
        background: white;
        border-radius: 24px;
        max-width: 900px;
        width: 100%;
        max-height: 90vh;
        overflow: hidden;
        box-shadow: 0 25px 50px rgba(0, 0, 0, 0.25);
        animation: modalSlideIn 0.3s ease;
    }

    @keyframes modalSlideIn {
        from {
            opacity: 0;
            transform: translateY(20px) scale(0.95);
        }
        to {
            opacity: 1;
            transform: translateY(0) scale(1);
        }
    }

    .modal-content {
        display: none;
    }

    .modal-content.active {
        display: block;
    }

    .modal-header {
        display: flex;
        justify-content: space-between;
        align-items: flex-start;
        padding: 28px 32px;
        background: linear-gradient(135deg, var(--primary) 0%, var(--secondary) 100%);
        color: white;
    }

    .modal-title {
        font-size: 1.5rem;
        font-weight: 700;
        margin-bottom: 4px;
    }

    .modal-subtitle {
        font-size: 1rem;
        opacity: 0.9;
    }

    .modal-close {
        background: rgba(255, 255, 255, 0.2);
        border: none;
        color: white;
        font-size: 1.8rem;
        width: 40px;
        height: 40px;
        border-radius: 12px;
        cursor: pointer;
        transition: all 0.2s ease;
        display: flex;
        align-items: center;
        justify-content: center;
        line-height: 1;
    }

    .modal-close:hover {
        background: rgba(255, 255, 255, 0.3);
        transform: rotate(90deg);
    }

    .modal-body {
        padding: 28px 32px;
        max-height: calc(90vh - 120px);
        overflow-y: auto;
    }

    .markers-full {
        background: rgba(13, 148, 136, 0.05);
        border-radius: 12px;
        padding: 16px;
        margin-bottom: 24px;
        border-left: 4px solid var(--primary);
    }

    .markers-full h4 {
        font-size: 0.85rem;
        color: var(--text-secondary);
        text-transform: uppercase;
        letter-spacing: 0.5px;
        margin-bottom: 8px;
    }

    .markers-text {
        font-size: 0.9rem;
        color: var(--text-primary);
        word-break: break-word;
        line-height: 1.7;
    }

    .modal-section {
        margin-bottom: 24px;
        border-radius: 16px;
        overflow: hidden;
        border: 1px solid var(--border);
    }

    .section-title {
        display: flex;
        align-items: center;
        gap: 10px;
        padding: 16px 20px;
        font-size: 1.05rem;
        font-weight: 600;
        margin: 0;
    }

    .section-content {
        padding: 20px;
        background: white;
        font-size: 0.95rem;
        line-height: 1.8;
    }

    .annotation-section .section-title {
        background: linear-gradient(135deg, rgba(13, 148, 136, 0.1), rgba(13, 148, 136, 0.05));
        color: var(--primary-dark);
    }

    .validator-section .section-title {
        background: linear-gradient(135deg, rgba(5, 150, 105, 0.1), rgba(5, 150, 105, 0.05));
        color: #047857;
    }

    .formatting-section .section-title {
        background: linear-gradient(135deg, rgba(245, 158, 11, 0.1), rgba(245, 158, 11, 0.05));
        color: #b45309;
    }

    .scoring-section .section-title {
        background: linear-gradient(135deg, rgba(99, 102, 241, 0.1), rgba(99, 102, 241, 0.05));
        color: #4338ca;
    }

    .validation-status {
        display: inline-flex;
        align-items: center;
        gap: 8px;
        padding: 10px 16px;
        border-radius: 10px;
        font-weight: 600;
        margin-bottom: 12px;
    }

    .validation-status.passed {
        background: rgba(5, 150, 105, 0.1);
        color: #047857;
    }

    .validation-status.failed {
        background: rgba(239, 68, 68, 0.1);
        color: #dc2626;
    }

    .validator-content {
        color: var(--text-secondary);
        font-size: 0.9rem;
    }

    .json-summary {
        display: flex;
        flex-direction: column;
        gap: 16px;
    }

    .summary-field {
        display: flex;
        flex-direction: column;
        gap: 8px;
    }

    .summary-field .field-label {
        font-size: 0.8rem;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        color: var(--accent);
        font-weight: 600;
    }

    .summary-field .field-value {
        color: var(--text-primary);
    }

    .summary-field .main-type {
        font-size: 1.2rem;
        font-weight: 700;
        color: var(--primary-dark);
    }

    .type-list {
        list-style: none;
        display: flex;
        flex-direction: column;
        gap: 8px;
    }

    .type-list li {
        padding: 8px 12px;
        background: rgba(13, 148, 136, 0.05);
        border-radius: 8px;
        font-size: 0.9rem;
    }

    .type-list li.empty {
        color: var(--text-muted);
        font-style: italic;
        background: rgba(100, 116, 139, 0.05);
    }

    /* Ranked subtypes in modal type-list */
    .type-list li.rank-1 {
        background: rgba(13, 148, 136, 0.15);
        color: #0d9488;
        font-weight: 600;
        border-left: 3px solid #0d9488;
    }

    .type-list li.rank-2 {
        background: rgba(20, 184, 166, 0.12);
        color: #134e4a;
        border-left: 3px solid #14b8a6;
    }

    .type-list li.rank-3 {
        background: rgba(94, 234, 212, 0.12);
        color: #64748b;
        border-left: 3px solid #5eead4;
    }

    /* Merged cell type groupings styles */
    .merged-section {
        margin-top: 8px;
    }

    .merged-groupings {
        display: flex;
        flex-direction: column;
        gap: 6px;
    }

    .merged-item {
        display: flex;
        align-items: center;
        gap: 8px;
    }

    .merged-label {
        font-size: 0.7rem;
        font-weight: 600;
        color: var(--text-secondary);
        text-transform: uppercase;
        min-width: 55px;
    }

    .merged-tag {
        padding: 4px 10px;
        border-radius: 6px;
        font-size: 0.8rem;
        font-weight: 500;
        max-width: 100%;
        overflow: hidden;
        text-overflow: ellipsis;
        white-space: nowrap;
    }

    .merged-broad {
        background: rgba(13, 148, 136, 0.15);
        color: #0d9488;
        border: 1px solid rgba(13, 148, 136, 0.3);
    }

    .merged-detailed {
        background: rgba(5, 150, 105, 0.15);
        color: #059669;
        border: 1px solid rgba(5, 150, 105, 0.3);
    }

    .merged-specific {
        background: rgba(245, 158, 11, 0.15);
        color: #b45309;
        border: 1px solid rgba(245, 158, 11, 0.3);
    }

    /* Merged groupings section in modal */
    .merged-groupings-section .section-title {
        background: linear-gradient(135deg, rgba(139, 92, 246, 0.1), rgba(139, 92, 246, 0.05));
        color: #7c3aed;
    }

    .merged-detail-grid {
        display: flex;
        flex-direction: column;
        gap: 12px;
    }

    .merged-detail-item {
        display: flex;
        flex-direction: column;
        gap: 4px;
        padding: 12px;
        background: rgba(139, 92, 246, 0.05);
        border-radius: 8px;
        border-left: 3px solid #7c3aed;
    }

    .merged-level {
        font-size: 0.75rem;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        color: #7c3aed;
        font-weight: 600;
    }

    .merged-value {
        font-size: 0.95rem;
        color: var(--text-primary);
    }

    .json-raw {
        background: #1e293b;
        color: #e2e8f0;
        padding: 16px;
        border-radius: 8px;
        font-size: 0.85rem;
        overflow-x: auto;
    }

    .step-header {
        color: var(--primary);
        font-size: 1rem;
        margin-top: 20px;
        margin-bottom: 10px;
        padding-bottom: 6px;
        border-bottom: 2px solid rgba(13, 148, 136, 0.2);
    }

    .section-content h2, .section-content h3, .section-content h4 {
        color: var(--text-primary);
        margin-top: 16px;
        margin-bottom: 8px;
    }

    .section-content ul {
        margin: 12px 0;
        padding-left: 24px;
    }

    .section-content li {
        margin-bottom: 6px;
    }

    .section-content p {
        margin-bottom: 12px;
    }

    .section-content strong {
        color: var(--text-primary);
    }

    /* Print Styles */
    @media print {
        .controls-bar,
        .view-analysis-btn,
        .modal-overlay {
            display: none !important;
        }

        .cluster-card {
            break-inside: avoid;
            box-shadow: none;
            border: 1px solid #ddd;
        }

        .report-header {
            background: var(--primary) !important;
            -webkit-print-color-adjust: exact;
            print-color-adjust: exact;
        }
    }

    /* Responsive */
    @media (max-width: 768px) {
        .container {
            padding: 12px;
        }

        .report-header {
            padding: 24px;
            border-radius: 16px;
        }

        .report-title {
            font-size: 1.8rem;
        }

        .cluster-grid {
            grid-template-columns: 1fr;
        }

        .controls-bar {
            flex-direction: column;
        }

        .search-box {
            width: 100%;
        }

        .filter-select {
            width: 100%;
        }

        .modal-container {
            max-height: 95vh;
            border-radius: 16px;
        }

        .modal-header {
            padding: 20px;
        }

        .modal-body {
            padding: 20px;
        }
    }

    /* No results message */
    .no-results {
        text-align: center;
        padding: 60px 20px;
        color: var(--text-secondary);
    }

    .no-results svg {
        width: 64px;
        height: 64px;
        color: var(--text-muted);
        margin-bottom: 16px;
    }

    .no-results h3 {
        font-size: 1.2rem;
        margin-bottom: 8px;
        color: var(--text-primary);
    }
    `;
}

// ----------------- JavaScript Code -----------------

/**
 * Return the JavaScript code for interactivity
 * @returns {string} JavaScript code
 */
function getJavaScript() {
    return `
    let currentModalIndex = null;

    function openModal(index) {
        const overlay = document.getElementById('modal-overlay');
        const allContents = document.querySelectorAll('.modal-content');

        allContents.forEach(c => c.classList.remove('active'));

        const content = document.getElementById('modal-content-' + index);
        if (content) {
            content.classList.add('active');
            overlay.classList.add('active');
            currentModalIndex = index;
            document.body.style.overflow = 'hidden';
        }
    }

    function closeModal() {
        const overlay = document.getElementById('modal-overlay');
        overlay.classList.remove('active');
        currentModalIndex = null;
        document.body.style.overflow = '';
    }

    // Close modal on overlay click
    document.getElementById('modal-overlay').addEventListener('click', function(e) {
        if (e.target === this) {
            closeModal();
        }
    });

    // Close modal on Escape key
    document.addEventListener('keydown', function(e) {
        if (e.key === 'Escape' && currentModalIndex !== null) {
            closeModal();
        }
    });

    // Search functionality
    function filterCards() {
        const searchTerm = document.getElementById('search-input').value.toLowerCase();
        const tissueFilter = document.getElementById('filter-tissue').value;
        const speciesFilter = document.getElementById('filter-species').value;
        const modelFilter = document.getElementById('filter-model').value;
        const providerFilter = document.getElementById('filter-provider').value;

        const cards = document.querySelectorAll('.cluster-card');
        let visibleCount = 0;

        cards.forEach(card => {
            const text = card.textContent.toLowerCase();
            const tissue = card.dataset.tissue;
            const species = card.dataset.species;
            const model = card.dataset.model;
            const provider = card.dataset.provider;

            const matchesSearch = searchTerm === '' || text.includes(searchTerm);
            const matchesTissue = tissueFilter === '' || tissue === tissueFilter;
            const matchesSpecies = speciesFilter === '' || species === speciesFilter;
            const matchesModel = modelFilter === '' || model === modelFilter;
            const matchesProvider = providerFilter === '' || provider === providerFilter;

            if (matchesSearch && matchesTissue && matchesSpecies && matchesModel && matchesProvider) {
                card.style.display = '';
                visibleCount++;
            } else {
                card.style.display = 'none';
            }
        });

        // Show/hide no results message
        const noResults = document.getElementById('no-results');
        if (noResults) {
            noResults.style.display = visibleCount === 0 ? 'block' : 'none';
        }

        // Update visible count
        const countEl = document.getElementById('visible-count');
        if (countEl) {
            countEl.textContent = visibleCount;
        }
    }

    // Clear filters
    function clearFilters() {
        document.getElementById('search-input').value = '';
        document.getElementById('filter-tissue').value = '';
        document.getElementById('filter-species').value = '';
        document.getElementById('filter-model').value = '';
        document.getElementById('filter-provider').value = '';
        filterCards();
    }

    // Export to PDF (print)
    function exportToPDF() {
        window.print();
    }

    // Initialize search input listener
    document.getElementById('search-input').addEventListener('input', filterCards);
    `;
}

// ----------------- Main Report Generation Functions -----------------

/**
 * Generate an interactive HTML report directly from data (list of dictionaries)
 * This preserves newlines in conversation history for better formatting.
 *
 * @param {Array} rows - List of dictionaries containing analysis data for each cluster
 * @param {string|null} outputPath - Output path for HTML file (null for browser - returns string)
 * @param {string} reportTitle - Title displayed in the report header
 * @returns {string} HTML content string
 */
export function generateBatchHtmlReportFromData(rows, outputPath = null, reportTitle = 'CASSIA Batch Analysis Report') {
    if (!rows || rows.length === 0) {
        throw new Error('No data provided');
    }

    // Extract unique values for filters
    const tissues = [...new Set(rows.map(r => String(r['Tissue'] || '')).filter(t => t))].sort();
    const speciesList = [...new Set(rows.map(r => String(r['Species'] || '')).filter(s => s))].sort();
    const models = [...new Set(rows.map(r => String(r['Model'] || '')).filter(m => m))].sort();
    const providers = [...new Set(rows.map(r => String(r['Provider'] || '')).filter(p => p))].sort();

    // Generate timestamp
    const timestamp = new Date().toLocaleString();

    // Build filter options HTML
    const tissueOptions = '<option value="">All Tissues</option>' +
        tissues.map(t => `<option value="${escapeHtml(t)}">${escapeHtml(t)}</option>`).join('');
    const speciesOptions = '<option value="">All Species</option>' +
        speciesList.map(s => `<option value="${escapeHtml(s)}">${escapeHtml(s)}</option>`).join('');
    const modelOptions = '<option value="">All Models</option>' +
        models.map(m => `<option value="${escapeHtml(m)}">${escapeHtml(m)}</option>`).join('');
    const providerOptions = '<option value="">All Providers</option>' +
        providers.map(p => `<option value="${escapeHtml(p)}">${escapeHtml(p)}</option>`).join('');

    // Generate cluster cards
    const cardsHtml = rows.map((row, i) => generateClusterCard(row, i)).join('\n');

    // Generate modal contents
    const modalsHtml = rows.map((row, i) => generateModalContent(row, i)).join('\n');

    // Get primary tissue/species/model for header
    const primaryTissue = tissues[0] || 'N/A';
    const primarySpecies = speciesList[0] || 'N/A';

    // Build complete HTML
    const htmlContent = `<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>${escapeHtml(reportTitle)}</title>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap" rel="stylesheet">
    <style>
${getCssStyles()}
    </style>
</head>
<body>
    <div class="container">
        <header class="report-header">
            <div class="header-content">
                <h1 class="report-title">${escapeHtml(reportTitle)}</h1>
                <p class="report-subtitle">Comprehensive Cell Type Annotation Analysis</p>
                <div class="report-meta">
                    <span class="meta-badge">
                        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                            <rect x="3" y="4" width="18" height="18" rx="2" ry="2"></rect>
                            <line x1="16" y1="2" x2="16" y2="6"></line>
                            <line x1="8" y1="2" x2="8" y2="6"></line>
                            <line x1="3" y1="10" x2="21" y2="10"></line>
                        </svg>
                        ${timestamp}
                    </span>
                    <span class="meta-badge">
                        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                            <path d="M20.84 4.61a5.5 5.5 0 0 0-7.78 0L12 5.67l-1.06-1.06a5.5 5.5 0 0 0-7.78 7.78l1.06 1.06L12 21.23l7.78-7.78 1.06-1.06a5.5 5.5 0 0 0 0-7.78z"></path>
                        </svg>
                        ${escapeHtml(primaryTissue)}
                    </span>
                    <span class="meta-badge">
                        <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                            <circle cx="12" cy="12" r="10"></circle>
                            <path d="M2 12h20"></path>
                            <path d="M12 2a15.3 15.3 0 0 1 4 10 15.3 15.3 0 0 1-4 10 15.3 15.3 0 0 1-4-10 15.3 15.3 0 0 1 4-10z"></path>
                        </svg>
                        ${escapeHtml(primarySpecies)}
                    </span>
                </div>
            </div>
        </header>

        <div class="stats-bar">
            <div class="stat-card">
                <div class="stat-value">${rows.length}</div>
                <div class="stat-label">Total Clusters</div>
            </div>
            <div class="stat-card">
                <div class="stat-value" id="visible-count">${rows.length}</div>
                <div class="stat-label">Showing</div>
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

        <div class="controls-bar">
            <div class="search-box">
                <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <circle cx="11" cy="11" r="8"></circle>
                    <path d="M21 21l-4.35-4.35"></path>
                </svg>
                <input type="text" id="search-input" placeholder="Search clusters, cell types, markers...">
            </div>
            <select class="filter-select" id="filter-tissue" onchange="filterCards()">
                ${tissueOptions}
            </select>
            <select class="filter-select" id="filter-species" onchange="filterCards()">
                ${speciesOptions}
            </select>
            <select class="filter-select" id="filter-model" onchange="filterCards()">
                ${modelOptions}
            </select>
            <select class="filter-select" id="filter-provider" onchange="filterCards()">
                ${providerOptions}
            </select>
            <button class="btn btn-secondary" onclick="clearFilters()">
                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <path d="M18 6L6 18"></path>
                    <path d="M6 6l12 12"></path>
                </svg>
                Clear
            </button>
            <button class="btn btn-primary" onclick="exportToPDF()">
                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                    <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path>
                    <polyline points="7 10 12 15 17 10"></polyline>
                    <line x1="12" y1="15" x2="12" y2="3"></line>
                </svg>
                Export PDF
            </button>
        </div>

        <div class="cluster-grid">
            ${cardsHtml}
        </div>

        <div id="no-results" class="no-results" style="display: none;">
            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                <circle cx="11" cy="11" r="8"></circle>
                <path d="M21 21l-4.35-4.35"></path>
            </svg>
            <h3>No clusters found</h3>
            <p>Try adjusting your search or filters</p>
        </div>
    </div>

    <div class="modal-overlay" id="modal-overlay">
        <div class="modal-container">
            ${modalsHtml}
        </div>
    </div>

    <script>
${getJavaScript()}
    </script>
</body>
</html>`;

    console.log(`HTML report generated: ${rows.length} clusters`);
    return htmlContent;
}

export default { generateBatchHtmlReportFromData };
