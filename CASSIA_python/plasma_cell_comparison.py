#!/usr/bin/env python3
"""
Plasma cell subtype comparison using CASSIA enhanced functionality.
"""

import CASSIA

# The marker here are copy from CASSIA's previous results.
marker = "IGLL5, IGLV6-57, JCHAIN, FAM92B, IGLC3, IGLC2, IGHV3-7, IGKC, TNFRSF17, IGHG1, AC026369.3, IGHV3-23, IGKV4-1, IGKV1-5, IGHA1, IGLV3-1, IGLV2-11, MYL2, MZB1, IGHG3, IGHV3-74, IGHM, ANKRD36BP2, AMPD1, IGKV3-20, IGHA2, DERL3, AC104699.1, LINC02362, AL391056.1, LILRB4, CCL3, BMP6, UBE2QL1, LINC00309, AL133467.1, GPRC5D, FCRL5, DNAAF1, AP002852.1, AC007569.1, CXorf21, RNU1-85P, U62317.4, TXNDC5, LINC02384, CCR10, BFSP2, APOBEC3A, AC106897.1"

print("🧬 Running Plasma Cell Subtype Comparison")
print("=" * 50)

results = CASSIA.compareCelltypes(
    tissue="large intestine",
    celltypes=["Plasma Cells", "IgA-secreting Plasma Cells", "IgG-secreting Plasma Cells", "IgM-secreting Plasma Cells"],
    marker_set=marker,
    species="human",
    output_file="plasama_cell_subtype.csv"
)

# Extract all cell types and models for summary table
celltypes = ["Plasma Cells", "IgA-secreting Plasma Cells", "IgG-secreting Plasma Cells", "IgM-secreting Plasma Cells"]
models = []
scores_matrix = {}

for result in results['results']:
    if result['status'] == 'success':
        model = result['model']
        models.append(model)
        scores_matrix[model] = {}
        
        extracted_scores = result.get('extracted_scores', {})
        for celltype in celltypes:
            if celltype in extracted_scores:
                scores_matrix[model][celltype] = extracted_scores[celltype]['score']
            else:
                scores_matrix[model][celltype] = 'N/A'

print("\n📊 SUMMARY SCORE TABLE")
print("=" * 80)

# Header
print(f"{'Cell Type':<30}", end="")
for model in models:
    model_short = model.split('/')[-1] if '/' in model else model
    print(f"{model_short:<20}", end="")
print()

print("-" * 80)

# Scores for each cell type
for celltype in celltypes:
    print(f"{celltype:<30}", end="")
    for model in models:
        score = scores_matrix.get(model, {}).get(celltype, 'N/A')
        print(f"{score:<20}", end="")
    print()

print("\n" + "=" * 80)
print("📝 DETAILED REASONING BY MODEL")
print("=" * 80)

for result in results['results']:
    model = result['model']
    status = result['status']
    extracted_scores = result.get('extracted_scores', {})
    
    print(f"\n🤖 MODEL: {model}")
    print("-" * 60)
    
    if extracted_scores:
        for celltype in celltypes:
            if celltype in extracted_scores:
                score = extracted_scores[celltype]['score']
                reasoning = extracted_scores[celltype]['reasoning']
                
                print(f"\n📈 {celltype}: Score {score}/100")
                print("💭 Reasoning:")
                
                # Format reasoning with proper line breaks and indentation
                import textwrap
                wrapped_text = textwrap.fill(reasoning, width=70, initial_indent="   ", subsequent_indent="   ")
                print(wrapped_text)
                print()
            else:
                print(f"\n📈 {celltype}: N/A")
                print("💭 Reasoning: No data available")
                print()
    else:
        print("  ⚠️  No structured scores extracted")

print("✅ Comparison completed!")
print(f"\n📁 Generated Files:")
print(f"  📊 CSV Report: plasama_cell_subtype.csv")
print(f"  🌐 HTML Report: plasama_cell_subtype_report.html") 