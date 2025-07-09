# Enhanced Cell Type Comparison Documentation

## Overview

The CASSIA cell type comparison functionality has been significantly enhanced with structured output, intelligent extraction, and beautiful HTML visualization capabilities.

## 🌟 New Features

### 1. **Structured Output with Tags**

The AI models now provide responses in a structured format using XML-like tags:

```xml
<celltype>T cell</celltype>
<reasoning>
CD3 is the definitive T cell marker, present on all T cells. CD4 and CD8 
identify T helper and cytotoxic T cell subsets respectively.
</reasoning>
<score>95</score>
```

### 2. **Intelligent Score & Reasoning Extraction**

- Automatically parses structured responses from any AI model
- Extracts both numerical scores (0-100) and detailed reasoning
- Robust fallback mechanisms for malformed responses
- Handles multiple response formats gracefully

### 3. **Interactive HTML Reports**

- Beautiful, professional-looking HTML reports
- Color-coded scoring (high/medium/low)
- Expandable sections for raw model responses
- Cross-model comparison tables
- Mobile-friendly responsive design

### 4. **Enhanced CSV Output**

- Structured columns for each cell type's score and reasoning
- Separate columns: `[CellType]_score`, `[CellType]_reasoning`
- Raw model responses included for transparency
- Easy analysis with data science tools

## 🚀 Usage

### Basic Usage

```python
import CASSIA

# Enhanced comparison with HTML visualization
results = CASSIA.compareCelltypes(
    tissue="blood",
    celltypes=["T cell", "B cell", "NK cell"],
    marker_set="CD3, CD4, CD8, CD19, CD20, CD56",
    species="human",
    generate_html_report=True  # NEW: Creates HTML visualization
)

# Access structured results
for result in results['results']:
    model = result['model']
    scores = result['extracted_scores']
    for celltype, data in scores.items():
        print(f"{model} - {celltype}: Score {data['score']}")
        print(f"Reasoning: {data['reasoning']}")
```

### Function Parameters

```python
compareCelltypes(
    tissue,                    # str: Tissue type
    celltypes,                 # list: 2-4 cell types to compare  
    marker_set,                # str: Comma-separated marker genes
    species="human",           # str: Species being analyzed
    model_list=None,           # list: Models to use (optional)
    output_file=None,          # str: Output filename (optional)
    generate_html_report=True  # bool: Generate HTML report (NEW!)
)
```

### Return Value

```python
{
    'results': [
        {
            'model': 'anthropic/claude-3.5-sonnet',
            'extracted_scores': {
                'T cell': {
                    'score': '95',
                    'reasoning': 'CD3 is the definitive T cell marker...'
                },
                'B cell': {
                    'score': '70', 
                    'reasoning': 'CD19 and CD20 are B cell markers...'
                }
            },
            'response': 'Full raw response...',
            'status': 'success'
        }
    ],
    'csv_file': 'comparison_results.csv',
    'html_content': '<html>...</html>'
}
```

## 📁 Output Files

### CSV File Structure

| Column | Description |
|--------|-------------|
| `model` | AI model used |
| `tissue` | Tissue type analyzed |
| `species` | Species analyzed |
| `cell_types` | Cell types compared |
| `status` | Success/error status |
| `[CellType]_score` | Score for each cell type |
| `[CellType]_reasoning` | Reasoning for each cell type |
| `raw_response` | Complete model response |

### HTML Report Features

- **Header Section**: Professional title and analysis info
- **Model Sections**: Each AI model gets its own section with:
  - Cell type cards showing scores and reasoning
  - Color-coded score badges
  - Expandable raw response sections
- **Summary Table**: Cross-model comparison with color coding:
  - 🟢 High scores (70+): Green highlighting
  - 🟡 Medium scores (40-69): Yellow highlighting  
  - 🔴 Low scores (<40): Red highlighting
- **Interactive Elements**: Toggle buttons to show/hide raw responses

## 🧪 Testing

### Comprehensive Test Suite

```bash
python test_enhanced_comparison.py
```

This runs:
- ✅ Score extraction testing with mock data
- ✅ HTML report generation testing  
- ✅ Edge case and error handling testing
- ✅ Live API comparison testing (requires API key)

### Simple Example

```bash
python example_enhanced_comparison.py
```

Demonstrates:
- Basic usage patterns
- Extraction functionality
- HTML report generation
- Error handling

## 🔧 Technical Details

### Extraction Function

```python
from CASSIA.cell_type_comparison import extract_celltype_scores

# Extract scores from any response format
scores = extract_celltype_scores(response_text, celltypes)
```

### HTML Generation Function

```python
from CASSIA.cell_type_comparison import generate_comparison_html_report

# Generate standalone HTML report
html_content = generate_comparison_html_report(results, "report.html")
```

## 🎯 Use Cases

1. **Research Publications**: Generate professional reports for papers
2. **Collaborative Analysis**: Share interactive HTML reports with colleagues  
3. **Method Comparison**: Compare multiple AI models on the same data
4. **Quality Control**: Detailed reasoning helps validate results
5. **Educational**: Clear explanations help students understand cell type identification

## 🔄 Backward Compatibility

- All existing code using `CASSIA.compareCelltypes()` continues to work
- New features are opt-in via the `generate_html_report` parameter
- Enhanced output structure while maintaining the same interface
- Previous CSV format is still supported (with additional columns)

## 📊 Example Results

### Command Line Output
```
Model: anthropic/claude-3.5-haiku
Extracted scores: {
    'T cell': {
        'score': '95', 
        'reasoning': 'CD3 is the definitive T cell marker...'
    },
    'B cell': {
        'score': '70',
        'reasoning': 'CD19 and CD20 are B cell markers...'
    }
}
```

### HTML Report Preview
The generated HTML reports feature:
- Professional styling with gradients and shadows
- Responsive grid layout for cell type cards
- Interactive toggle buttons for detailed views
- Color-coded summary tables
- Mobile-friendly design

## 🚀 Future Enhancements

Potential improvements:
- Integration with additional AI providers
- Advanced visualization options (charts, graphs)
- Batch processing for multiple marker sets
- Export to additional formats (PDF, Word)
- Custom scoring rubrics and thresholds

---

**Version**: Enhanced in CASSIA v0.2.18+  
**Dependencies**: pandas, requests, re, json  
**License**: Same as CASSIA main package 