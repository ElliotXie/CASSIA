##final version
import pandas as pd
import scanpy as sc

adata = sc.read_h5ad("C:/Users/ellio/Downloads/GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad")
# Subset the AnnData object to only include lung tissue
lung_adata = adata[adata.obs['tissue'] == 'lung'].copy()


#find the grouping column
print(lung_adata.obs['Broad cell type'].unique())
print(lung_adata.obs['Granular cell type'].unique())
print(lung_adata.obs['Cell types level 2'].unique())
print(lung_adata.obs['Cell types level 3'].unique())

#generate UMAP
sc.pl.umap(adata, color=['Broad cell type', 'Granular cell type'])

# Normalize the data by total counts per cell and scale to 10,000 reads per cell
sc.pp.normalize_total(lung_adata, target_sum=1e4)
# Log-transform the data after adding a pseudocount of 1
sc.pp.log1p(lung_adata)
#remove batch effect
sc.pp.combat(lung_adata, key='batch') 




def analyze_markers(adata, groupby, n_genes=10):
    # Perform differential expression analysis
    sc.tl.rank_genes_groups(adata, groupby=groupby, method='t-test', use_raw=False)
    
    # Get the default Scanpy results
    default_markers = sc.get.rank_genes_groups_df(adata, group=None)
    
    # Modified version
    modified_markers = default_markers.copy()
    modified_markers = modified_markers.sort_values(['group', 'logfoldchanges'], ascending=[True, False])
    
    top_markers = modified_markers.groupby('group').apply(
        lambda x: ', '.join(x['names'].head(n_genes))
    ).reset_index()
    top_markers.columns = [groupby, f'Top {n_genes} Markers']
    
    return default_markers, top_markers

# List of annotation levels
annotation_levels = ['Broad cell type', 'Granular cell type', 'Cell types level 2', 'Cell types level 3']

# Dictionaries to store results for each annotation level
default_results = {}
modified_results = {}

# Perform analysis for each annotation level
for level in annotation_levels:
    print(f"Analyzing {level}...")
    default_df, modified_df = analyze_markers(lung_adata, level)
    default_results[level] = default_df
    modified_results[level] = modified_df
    print(f"Analysis for {level} completed.")

# Print and export results for each level
for level in annotation_levels:
    print(f"\n--- Results for {level} ---")
    
    print("\nDefault Scanpy output (first 10 rows):")
    print(default_results[level].head(10))
    
    print("\nModified output:")
    print(modified_results[level].head())
    
    # Export both DataFrames to separate CSV files
    default_results[level].to_csv(f"default_markers_{level.replace(' ', '_').lower()}.csv", index=False)
    modified_results[level].to_csv(f"modified_markers_{level.replace(' ', '_').lower()}.csv", index=False)
    
    print(f"Results for {level} exported to CSV files.")

print("\nAll analyses completed and results exported.")