import os
import pandas as pd
import re
import numpy as np
from typing import List, Dict, Any, Union, Optional

def debug_gene_extraction(
    marker_data: Union[str, pd.DataFrame],
    test_genes: List[str] = ["CD133", "CD9", "ChAT", "DCLK1", "EDNRB"],
    case_sensitive: bool = True
):
    """
    Debug function to test gene extraction from marker data.
    
    Args:
        marker_data: Path to marker data CSV or DataFrame
        test_genes: Sample genes to test for presence
        case_sensitive: Whether to use case-sensitive matching
    """
    print(f"=== Gene Extraction Debug ===")
    
    # Load the marker data
    try:
        if isinstance(marker_data, pd.DataFrame):
            marker_df = marker_data
            print(f"Using provided DataFrame with shape: {marker_df.shape}")
        else:
            print(f"Loading marker data from: {marker_data}")
            try:
                # Try with index_col=0 first
                marker_df = pd.read_csv(marker_data, index_col=0)
                print(f"Loaded with index_col=0, shape: {marker_df.shape}")
            except Exception as e1:
                # If that fails, try without index_col
                try:
                    marker_df = pd.read_csv(marker_data)
                    print(f"Loaded without index_col, shape: {marker_df.shape}")
                except Exception as e2:
                    print(f"Error loading marker data: {str(e1)}\nSecond attempt: {str(e2)}")
                    return
        
        # Check index type
        print(f"Index type: {type(marker_df.index).__name__}")
        print(f"First 5 index values: {list(marker_df.index[:5])}")
        
        # Display column names
        print(f"Columns: {marker_df.columns.tolist()}")
        
        # Check if there's a gene symbol column instead of index
        potential_gene_cols = [col for col in marker_df.columns if 
                              any(keyword in col.lower() for keyword in 
                                 ['gene', 'symbol', 'name', 'id', 'feature'])]
        
        if potential_gene_cols:
            print(f"Potential gene name columns: {potential_gene_cols}")
            
        # Check for the test genes in the index (case-sensitive)
        found_exact = [gene for gene in test_genes if gene in marker_df.index]
        print(f"Found genes (exact match): {found_exact}")
        
        # Check for the test genes in the index (case-insensitive)
        if not case_sensitive:
            # Create a lowercase mapping of index values
            lower_index = {str(idx).lower(): idx for idx in marker_df.index}
            found_lower = [gene for gene in test_genes if str(gene).lower() in lower_index]
            
            if found_lower:
                print(f"Found genes (case-insensitive): {found_lower}")
                print(f"Actual index entries: {[lower_index[gene.lower()] for gene in found_lower if gene.lower() in lower_index]}")
        
        # Check for potential alternate names/cases in the index
        for gene in test_genes:
            pattern = re.compile(f".*{re.escape(gene)}.*", re.IGNORECASE)
            matches = [idx for idx in marker_df.index if pattern.match(str(idx))]
            if matches:
                print(f"Potential matches for {gene}: {matches}")
        
        # If there are potential gene columns, check there too
        for col in potential_gene_cols:
            found_in_col = []
            for gene in test_genes:
                matches = marker_df[marker_df[col].str.contains(gene, case=case_sensitive, na=False)]
                if not matches.empty:
                    found_in_col.append(gene)
            
            if found_in_col:
                print(f"Found in column '{col}': {found_in_col}")
        
        # Print first few rows for inspection
        print("\nSample data (first 5 rows):")
        print(marker_df.head())
        
    except Exception as e:
        print(f"Error during debug: {str(e)}")
        import traceback
        traceback.print_exc()

def analyze_conversation_for_genes(conversation: str):
    """
    Debug function to analyze how genes are extracted from a conversation.
    
    Args:
        conversation: Sample conversation to test
    """
    print(f"\n=== Conversation Gene Extraction Debug ===")
    
    # Print the conversation
    print(f"Conversation text: {conversation[:200]}...")
    
    # Test extraction with different patterns
    patterns = [
        (r'<check_genes>\s*(.*?)\s*</check_genes>', "Standard check_genes tags"),
        (r'<check_genes>(.*?)</check_genes>', "Without whitespace handling"),
        (r'check_genes:(.*?)(?:$|\n)', "Informal format with colon"),
        (r'genes to check:?(.*?)(?:$|\n)', "Natural language mention"),
        (r'(?:^|\n)\s*[\-\*]\s*([^:\n]*?)(?:$|\n)', "Bullet point format")
    ]
    
    for pattern, description in patterns:
        gene_lists = re.findall(pattern, conversation, re.DOTALL | re.IGNORECASE)
        print(f"\nPattern ({description}):")
        print(f"  Raw matches: {gene_lists}")
        
        if gene_lists:
            # Process gene lists (similar to the extraction function)
            all_genes = []
            for gene_list in gene_lists:
                # Clean and normalize the gene list
                cleaned_list = re.sub(r'[\]\[\)\(]', '', gene_list)
                cleaned_list = re.sub(r'\s+', ' ', cleaned_list)
                # Split by comma or space, depending on formatting
                genes = re.split(r',\s*|\s+', cleaned_list)
                genes = [g.strip() for g in genes if g.strip()]
                all_genes.extend(genes)
            
            print(f"  Extracted genes: {all_genes}")

# Function to test the actual filter_marker function logic
def test_filter_marker(marker_data: Union[str, pd.DataFrame], gene_list: List[str]):
    """
    Test the filter_marker function with specific genes.
    
    Args:
        marker_data: Path to marker data CSV or DataFrame
        gene_list: List of genes to test
    """
    print(f"\n=== Testing filter_marker Function ===")
    
    try:
        # Load marker data
        if isinstance(marker_data, pd.DataFrame):
            marker_df = marker_data
        else:
            try:
                marker_df = pd.read_csv(marker_data, index_col=0)
            except:
                marker_df = pd.read_csv(marker_data)
        
        print(f"Loaded marker data with shape: {marker_df.shape}")
        
        # Implement the filter_marker logic directly
        # Remove any 'Unnamed: 0' column if it exists
        if 'Unnamed: 0' in marker_df.columns:
            marker_df = marker_df.drop(columns=['Unnamed: 0'])
            print("Dropped 'Unnamed: 0' column")

        # Identify valid genes and NA genes
        valid_genes = []
        na_genes = []
        for gene in gene_list:
            if gene in marker_df.index:
                # Check if all values for this gene are NA
                gene_data = marker_df.loc[gene]
                print(f"Gene {gene} data: {gene_data}")
                if isinstance(gene_data, pd.Series):
                    if gene_data.isna().all() or (gene_data == 'NA').all():
                        na_genes.append(gene)
                        print(f"Gene {gene} has all NA values")
                    else:
                        valid_genes.append(gene)
                        print(f"Gene {gene} found with valid data")
                else:
                    # Handle DataFrames for multi-index
                    if gene_data.isna().all().all() or (gene_data == 'NA').all().all():
                        na_genes.append(gene)
                        print(f"Gene {gene} has all NA values")
                    else:
                        valid_genes.append(gene)
                        print(f"Gene {gene} found with valid data")
            else:
                na_genes.append(gene)
                print(f"Gene {gene} not found in index")
        
        print(f"\nValid genes: {valid_genes}")
        print(f"NA/Not found genes: {na_genes}")
        
        # Create result DataFrame with only valid genes
        if valid_genes:
            result = marker_df.loc[valid_genes].copy()
            print(f"\nResult DataFrame:\n{result}")
        else:
            # If no valid genes, create an empty dataframe with the same columns
            result = pd.DataFrame(columns=marker_df.columns)
            print("\nNo valid genes found, empty DataFrame created")
            
    except Exception as e:
        print(f"Error in test_filter_marker: {str(e)}")
        import traceback
        traceback.print_exc()

# Function to examine marker data structure
def examine_marker_structure(marker_data: Union[str, pd.DataFrame]):
    """
    Examine the structure of the marker data.
    
    Args:
        marker_data: Path to marker data CSV or DataFrame
    """
    print(f"\n=== Examining Marker Data Structure ===")
    
    try:
        # Load marker data
        if isinstance(marker_data, pd.DataFrame):
            marker_df = marker_data
            print(f"Using provided DataFrame")
        else:
            try:
                marker_df = pd.read_csv(marker_data, index_col=0)
                print(f"Loaded with index_col=0")
            except:
                marker_df = pd.read_csv(marker_data)
                print(f"Loaded without index_col")
        
        print(f"Shape: {marker_df.shape}")
        print(f"Index name: {marker_df.index.name}")
        print(f"Column names: {marker_df.columns.tolist()}")
        
        # Check index type and uniqueness
        print(f"Index type: {marker_df.index.dtype}")
        print(f"Index is unique: {marker_df.index.is_unique}")
        duplicate_index = marker_df.index.duplicated()
        if duplicate_index.any():
            print(f"Warning: {duplicate_index.sum()} duplicate index values found")
            print(f"Examples: {marker_df.index[duplicate_index][:5].tolist()}")
        
        # Check for column uniqueness
        for col in marker_df.columns:
            if pd.api.types.is_string_dtype(marker_df[col]):
                duplicate_values = marker_df[col].duplicated()
                if duplicate_values.any():
                    print(f"Column '{col}' has {duplicate_values.sum()} duplicate values")
                
                # Check if this column could be a gene identifier
                if col.lower() in ['gene', 'gene_id', 'gene_name', 'symbol']:
                    print(f"Potential gene identifier column: '{col}'")
                    print(f"First 5 values: {marker_df[col].head().tolist()}")
        
        # Data types
        print("\nColumn data types:")
        print(marker_df.dtypes)
        
        # Missing values
        print("\nMissing values per column:")
        print(marker_df.isna().sum())
        
        # Sample data
        print("\nFirst 5 rows:")
        print(marker_df.head())
        
    except Exception as e:
        print(f"Error examining marker structure: {str(e)}")
        import traceback
        traceback.print_exc()

# Testing all the debugging functions
def run_gene_diagnostics(
    marker_data: Union[str, pd.DataFrame],
    test_conversation: str = "",
    gene_list: List[str] = ["CD133", "CD9", "ChAT", "DCLK1", "EDNRB", "ERBB3", "FABP7", "GFAP", "KIT", "LGR5", "NGFR", "NKX2-2", "NOS1", "OLIG2", "PGP9.5", "PROM1", "RET", "S100B", "SOX9", "UCHL1", "VIP"]
):
    """Run all diagnostics on gene extraction"""
    print(f"Running comprehensive gene extraction diagnostics")
    
    # First examine the overall structure
    examine_marker_structure(marker_data)
    
    # Test case-sensitive and insensitive lookups
    debug_gene_extraction(marker_data, gene_list, case_sensitive=True)
    debug_gene_extraction(marker_data, gene_list, case_sensitive=False)
    
    # Test actual filter function
    test_filter_marker(marker_data, gene_list)
    
    if test_conversation:
        analyze_conversation_for_genes(test_conversation)
    
    print("\nDiagnostics complete.")

# Example usage:
if __name__ == "__main__":
    # You can run this script directly
    marker_path = "data/processed.csv"  # Path to your marker data
    sample_convo = """
    Based on the expression data, I'd like to check the following genes:
    <check_genes>CD133, CD9, ChAT, DCLK1, EDNRB</check_genes>
    
    Let's also look at these markers:
    <check_genes>ERBB3, FABP7, GFAP, KIT, LGR5</check_genes>
    """
    
    run_gene_diagnostics(marker_path, sample_convo) 