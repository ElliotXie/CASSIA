# CASSIA Annotation Boost - JavaScript Implementation

🚀 **Complete 100% Python-compatible implementation of CASSIA Annotation Boost functionality**

## Overview

The Annotation Boost module provides iterative marker analysis with hypothesis generation and gene checking for enhanced cell type annotation. This is a sophisticated system that allows for deep, methodical analysis of cell types through multiple rounds of gene expression evaluation.

## 🎯 Key Features

### ✅ **100% Python Compatibility**
- All functions maintain identical signatures and behavior to the Python implementation
- Same prompts, search strategies, and analysis workflows
- Compatible with all existing Python CASSIA workflows

### 🔬 **Iterative Analysis Methods**
- **Breadth-First Search**: Test multiple hypotheses simultaneously 
- **Depth-First Search**: Focus on one hypothesis at a time for deeper analysis
- **Conversation Summarization**: Intelligent summarization of previous analysis
- **Gene Expression Lookup**: Dynamic marker gene validation

### 📊 **Advanced Reporting**
- **HTML Summary Reports**: Professional scientific reports with styling
- **Raw Conversation Logs**: Complete analysis history in text format
- **Multiple Report Styles**: Per-iteration or gene-focused summaries
- **Interactive Visualizations**: Gene badges and structured layouts

### 🧬 **Flexible Data Handling**
- **CSV File Input**: Compatible with standard marker expression files
- **JavaScript Array Input**: Direct data input without file I/O
- **Multiple LLM Providers**: OpenAI, Anthropic, OpenRouter, custom endpoints
- **Conversation History Modes**: Full, summarized, or none

## 📁 Implementation Files

```
backend/src/
├── annotationBoost.js           # Main implementation (1,950+ lines)
├── llm_utils.js                 # LLM provider integration
└── ...                          # Other CASSIA modules

backend/test/
├── test_annotation_boost.js     # Comprehensive test suite
├── test_annotation_boost_import.js  # Quick import verification
└── ...                          # Other tests

backend/examples/
└── annotation_boost_example.js  # Usage examples and demos
```

## 🔧 Core Functions

### Main Functions
- `runCASSIAAnnotationboost()` - Standard iterative analysis
- `runCASSIAAnnotationboostAdditionalTask()` - Analysis with specific biological questions

### Analysis Functions
- `iterativeMarkerAnalysis()` - Core iterative workflow
- `extractGenesFromConversation()` - Parse gene requests from LLM responses
- `getMarkerInfo()` - Retrieve expression data for requested genes

### Prompt Generators
- `promptHypothesisGenerator()` - Breadth-first analysis prompts
- `promptHypothesisGeneratorDepthFirst()` - Focused analysis prompts
- `promptHypothesisGeneratorAdditionalTask()` - Task-specific prompts

### Report Generation
- `generateSummaryReport()` - Create structured analysis reports
- `formatSummaryToHtml()` - Convert to professional HTML format
- `saveRawConversationText()` - Export complete conversation logs

### Utility Functions
- `summarizeConversationHistory()` - Intelligent conversation summarization
- `prepareAnalysisData()` - Data loading and preparation
- `formatSummaryToHtml()` - HTML report formatting

## 🚀 Quick Start

### Basic Usage

```javascript
import { runCASSIAAnnotationboost } from './src/annotationBoost.js';

const result = await runCASSIAAnnotationboost({
    fullResultPath: "path/to/full_results.csv",
    marker: "path/to/markers.csv", 
    clusterName: "T cell",
    majorClusterInfo: "Human lung tissue scRNA-seq",
    outputName: "t_cell_analysis",
    numIterations: 3,
    provider: "openrouter",
    searchStrategy: "breadth"
});
```

### Advanced Usage with Custom Data

```javascript
// Custom marker data
const markerData = [
    { gene: "CD3D", avg_log2FC: 2.5, p_val_adj: 1.2e-15 },
    { gene: "CD3E", avg_log2FC: 2.1, p_val_adj: 3.4e-12 },
    // ... more genes
];

const result = await runCASSIAAnnotationboost({
    fullResultPath: "results.csv",
    marker: markerData, // Use array instead of file
    clusterName: "T cell",
    majorClusterInfo: "Human PBMC",
    outputName: "custom_analysis",
    numIterations: 2,
    searchStrategy: "depth", // Focus on one hypothesis
    reportStyle: "total_summary", // Gene-focused report
    conversationHistoryMode: "final" // Use summarization
});
```

### Analysis with Additional Tasks

```javascript
const result = await runCASSIAAnnotationboostAdditionalTask({
    fullResultPath: "results.csv",
    marker: "markers.csv",
    clusterName: "T cell", 
    majorClusterInfo: "Human tissue",
    outputName: "regulatory_analysis",
    additionalTask: "determine if these are regulatory T cells and assess activation state",
    searchStrategy: "depth",
    numIterations: 3
});
```

## 📊 Configuration Options

### Search Strategies
- **`"breadth"`**: Test multiple hypotheses per iteration (default)
- **`"depth"`**: Focus on one hypothesis, go deeper into subtypes

### Report Styles  
- **`"per_iteration"`**: Structured by analysis rounds (default)
- **`"total_summary"`**: Organized by gene groups and findings

### Conversation History Modes
- **`"final"`**: Use AI summarization of previous analysis (default)
- **`"full"`**: Include complete conversation history
- **`"none"`**: Start fresh without previous context

### LLM Providers
- **`"openrouter"`**: OpenRouter API (default)
- **`"openai"`**: OpenAI API
- **`"anthropic"`**: Anthropic API
- **Custom endpoints**: HTTP URLs

## 🧪 Testing

```bash
# Run comprehensive test suite
npm run test-annotation-boost

# Quick import verification
node test/test_annotation_boost_import.js

# View usage examples
node examples/annotation_boost_example.js
```

## 📈 Expected Output

### Success Response
```javascript
{
    status: 'success',
    raw_text_path: 'analysis_raw_conversation.txt',
    summary_report_path: 'analysis_summary.html', 
    execution_time: 45.2,
    analysis_text: 'FINAL ANNOTATION COMPLETED...'
}
```

### Generated Files
- **Raw Conversation**: Complete analysis dialogue with LLM
- **HTML Report**: Professional scientific report with styling
- **Gene Expression Data**: Retrieved marker information
- **Error Logs**: Diagnostic information if issues occur

## 🔬 Analysis Workflow

1. **Data Preparation**: Load CSV files or process custom data
2. **Context Setup**: Prepare cluster information and conversation history
3. **Initial Hypothesis**: Generate starting analysis prompt
4. **Iterative Analysis**:
   - LLM generates hypotheses and requests specific genes
   - System retrieves expression data for requested genes
   - LLM analyzes results and refines hypotheses
   - Process repeats until conclusion reached
5. **Report Generation**: Create structured HTML and text reports

## 💡 Best Practices

### Iteration Count
- **Quick Analysis**: 1-2 iterations for speed
- **Standard Analysis**: 3-5 iterations for thoroughness  
- **Deep Analysis**: 5-10 iterations for complex cases

### Search Strategy Selection
- **Breadth-First**: Use when you want to explore multiple possibilities
- **Depth-First**: Use when you have a strong hypothesis to validate

### Performance Tips
- Use `conversationHistoryMode: "none"` for faster execution
- Set lower `numIterations` for initial testing
- Use `temperature: 0` for consistent results

## 🔧 Integration

The annotation boost module integrates seamlessly with other CASSIA components:

```javascript
// Complete CASSIA workflow
import { 
    runCASSIABatch,           // Initial batch analysis
    runCASSIAScoreBatch,      // Quality scoring
    runCASSIAAnnotationboost, // Enhanced analysis
    runCASSIAGenerateScoreReport // Final reports
} from './index.js';

// 1. Initial analysis
const batchResults = await runCASSIABatch({...});

// 2. Score results  
const scoredResults = await runCASSIAScoreBatch({...});

// 3. Enhanced analysis for specific clusters
const enhancedResults = await runCASSIAAnnotationboost({...});

// 4. Generate comprehensive reports
const reports = await runCASSIAGenerateScoreReport({...});
```

## 🎉 Implementation Status

### ✅ Completed Features
- [x] All prompt generation functions (4 variants)
- [x] Iterative marker analysis engine
- [x] Gene expression data retrieval
- [x] Conversation summarization with LLM
- [x] HTML report generation with professional styling
- [x] Raw conversation text export
- [x] Error handling and diagnostics
- [x] Comprehensive test suite
- [x] Usage examples and documentation
- [x] Integration with main CASSIA module

### 🧬 Biological Accuracy
- Maintains all scientific prompts from Python implementation
- Preserves exact analysis workflows and decision logic
- Compatible with all cell type annotation use cases
- Supports advanced biological questions and hypotheses

## 📚 Documentation

- **Main README**: Core CASSIA functionality overview
- **Test README**: Comprehensive testing documentation  
- **Examples**: Real-world usage patterns
- **API Reference**: Complete function documentation

---

🎯 **Ready for Production**: The annotation boost module is fully implemented, tested, and ready for integration with any JavaScript-based single-cell analysis pipeline!