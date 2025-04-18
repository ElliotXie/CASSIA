import pandas as pd
import numpy as np
import json
import os
import requests
from typing import List, Dict, Union, Tuple, Optional, Any

class LLMEvaluator:
    """
    A class to evaluate cell type annotations using LLM-based scoring.
    Compares annotated cell types (from LLMs) against gold standard annotations.
    """
    
    def __init__(self, api_key: str = None, model: str = "deepseek/deepseek-chat-v3-0324"):
        """
        Initialize the LLM evaluator.
        
        Args:
            api_key (str): OpenRouter API key. If None, will try to get from environment.
            model (str): Model to use for evaluation from OpenRouter.
        """
        self.api_key = api_key or os.environ.get("OPENROUTER_API_KEY")
        if not self.api_key:
            raise ValueError("OpenRouter API key must be provided or set as OPENROUTER_API_KEY environment variable")
        
        self.model = model
    
    def get_single_celltype_prompts(self, 
                                   predicted_celltype: str, 
                                   gold_standard: str, 
                                   tissue: str, 
                                   species: str) -> Tuple[str, str]:
        """
        Generate system and user prompts for evaluating a single cell type prediction.
        
        Args:
            predicted_celltype (str): The cell type predicted by the annotation system
            gold_standard (str): The gold standard annotation
            tissue (str): The tissue context
            species (str): The species context
            
        Returns:
            Tuple[str, str]: System prompt and user prompt
        """
        
        system_prompt = """You are an expert cell biologist tasked with evaluating the accuracy of cell type annotations.
You will be given a predicted cell type and a gold standard (correct) cell type annotation.
Evaluate how accurate and specific the prediction is compared to the gold standard, considering cell ontology relationships.

Score the prediction on a scale of 0-5 where:
- 5: 100% exact match to the gold standard
- 4: General cell type is correct AND subtype is mostly correct but missing subtle details (e.g., activation state, exhaustion status, or minor variations like vascular endothelial cell vs. lymphatic endothelial cell)
- 3: General cell type is correct but subtype is not quite right, though still close on the ontology tree (e.g., macrophage vs. dendritic cell - both are myeloid cells)
- 2: Only the general cell type is correct but the subtype is very far in the ontology tree (e.g., T cell vs. myeloid cell - both are immune cells but far apart)
- 1: Incorrect general cell type, but prediction is somewhat related to the actual cell type in the broader classification
- 0: Completely irrelevant prediction that makes no sense in relation to the gold standard

Your response must include:
1. A numeric score (0-5)
2. A brief explanation of why you assigned this score, specifically referencing the cell ontology relationship
3. A JSON-formatted result with the format:
{"score": X, "explanation": "Your explanation here"}

Note that if a predicted cell type is more specific than the gold standard, give it a 4 and point that out in your explanation. For example, if true type is epithelicall cell and predicted is basal cell, that is a 4.

Note that when the orginal celltype is not specific, such as proliferating cell  cells, if the predicted celltype such as keratinocyte is a celltype known to be sharing that feature, give it a 4.

Note that if the predicted celltype is functionaly similar to the gold standard consider boost up the score.


"""
        
        user_prompt = f"""
Evaluate the following cell type annotation prediction against the gold standard:

Context:
- Tissue: {tissue}
- Species: {species}

Predicted cell type: {predicted_celltype}
Gold standard annotation: {gold_standard}

Please provide your score (0-5) and explanation, followed by a JSON in the format:
{{
  "score": X,
  "explanation": "Your explanation"
}}
"""
        
        return system_prompt, user_prompt
    
    def get_multiple_celltypes_prompts(self, 
                                      predicted_celltypes: List[str], 
                                      gold_standards: List[str], 
                                      tissue: str, 
                                      species: str) -> Tuple[str, str]:
        """
        Generate system and user prompts for evaluating multiple cell type predictions.
        
        Args:
            predicted_celltypes (List[str]): List of predicted cell types
            gold_standards (List[str]): List of corresponding gold standards
            tissue (str): The tissue context
            species (str): The species context
            
        Returns:
            Tuple[str, str]: System prompt and user prompt
        """
        
        system_prompt = """You are an expert cell biologist tasked with evaluating the accuracy of cell type annotations.
You will be given a set of predicted cell types and their corresponding gold standard (correct) annotations.
Evaluate how accurate and specific each prediction is compared to its gold standard, considering cell ontology relationships.

For each pair, score the prediction on a scale of 0-5 where:
- 5: 100% exact match to the gold standard，or if the predicted cell type is a more specific subtype of the gold standard.
- 4: General cell type is correct AND subtype is mostly correct but missing subtle details (e.g., activation state, exhaustion status, or minor variations like vascular endothelial cell vs. lymphatic endothelial cell)
- 3: The general cell type is correct, but the predicted subtype does not exactly match the gold standard. However, the prediction is still closely related in the cell ontology.
For example, the prediction is a differentiated cell type while the gold standard is an undifferentiated or progenitor cell type (e.g., “muscle cell progenitor cell” vs. “muscle cell”).
Another example: the predicted and gold standard cell types are derived from the same ancestor but belong to different lineages (e.g., “dendritic cell” vs. “macrophage”—both derived from monocytes).
- 2: Only the general cell type is correct but the subtype is very far in the ontology tree. For example, predicted is muscle cell but gold standard is neuron cell.
- 1: Incorrect general cell type, but prediction is somewhat related to the actual cell type in the broader classification
- 0: Completely irrelevant prediction that makes no sense in relation to the gold standard

Your response must include:
1. Individual scores for each prediction (0-5)
2. Brief explanations for your scoring decisions, specifically referencing cell ontology relationships
3. A JSON-formatted result with the format:
{
  "individual_scores": [score1, score2, ...],
  "explanations": ["explanation1", "explanation2", ...]
}

Additional Note:

If the predicted cell type is a more specific subtype of the gold standard, assign a score of 5.
Example: If the gold standard is "epithelial cell" and the prediction is "basal cell," assign a score of 5 and explain that "basal cell" is a specific subtype of "epithelial cell."

If the gold standard is a broad or feature-based label (such as "proliferating cell"), and the prediction is a cell type commonly known to exhibit that feature (for example, "keratinocyte" is known to proliferate), assign a score of 4.

"""
        
        # Construct the pairs for evaluation
        pairs = []
        for i, (pred, gold) in enumerate(zip(predicted_celltypes, gold_standards)):
            pairs.append(f"Pair {i+1}:\nPredicted: {pred}\nGold standard: {gold}\n")
        
        pairs_text = "\n".join(pairs)
        
        user_prompt = f"""
Evaluate the following cell type annotation predictions against their gold standards:

Context:
- Tissue: {tissue}
- Species: {species}

{pairs_text}

Please provide your evaluation with individual scores (0-5) and explanations, followed by a JSON in the specified format.
"""
        
        return system_prompt, user_prompt
    
    def evaluate_single_celltype(self, 
                                predicted_celltype: str, 
                                gold_standard: str, 
                                tissue: str, 
                                species: str) -> Dict[str, Any]:
        """
        Evaluate a single cell type annotation against a gold standard.
        
        Args:
            predicted_celltype (str): The cell type predicted by the annotation system
            gold_standard (str): The gold standard annotation
            tissue (str): The tissue context
            species (str): The species context
            
        Returns:
            Dict containing score and explanation
        """
        
        system_prompt, user_prompt = self.get_single_celltype_prompts(
            predicted_celltype, gold_standard, tissue, species
        )
        
        response = self._call_llm(system_prompt, user_prompt)
        result = self._extract_score_json(response)
        return result
    
    def evaluate_multiple_celltypes(self, 
                                   predicted_celltypes: List[str], 
                                   gold_standards: List[str], 
                                   tissue: str, 
                                   species: str) -> Dict[str, Any]:
        """
        Evaluate multiple cell type annotations against their gold standards.
        
        Args:
            predicted_celltypes (List[str]): List of predicted cell types
            gold_standards (List[str]): List of corresponding gold standards
            tissue (str): The tissue context
            species (str): The species context
            
        Returns:
            Dict containing individual scores and explanations
        """
        
        if len(predicted_celltypes) != len(gold_standards):
            raise ValueError("Length of predicted_celltypes must match length of gold_standards")
        
        system_prompt, user_prompt = self.get_multiple_celltypes_prompts(
            predicted_celltypes, gold_standards, tissue, species
        )
        
        response = self._call_llm(system_prompt, user_prompt)
        result = self._extract_multiple_scores_json(response)
        return result
    
    def batch_evaluate_from_dataframe(self, 
                                     df: pd.DataFrame, 
                                     predicted_col: str, 
                                     gold_col: str,
                                     tissue_col: str = None,
                                     species_col: str = None,
                                     default_tissue: str = "unknown",
                                     default_species: str = "human") -> pd.DataFrame:
        """
        Evaluate a batch of predictions from a dataframe.
        
        Args:
            df (pd.DataFrame): Dataframe containing predictions and gold standards
            predicted_col (str): Column name for predictions
            gold_col (str): Column name for gold standards
            tissue_col (str, optional): Column name for tissue context
            species_col (str, optional): Column name for species context
            default_tissue (str): Default tissue if tissue_col is None
            default_species (str): Default species if species_col is None
            
        Returns:
            pd.DataFrame with original data plus evaluation results
        """
        
        # Create copies of the dataframe columns for modification
        result_df = df.copy()
        result_df['evaluation_score'] = None
        result_df['evaluation_explanation'] = None
        
        for idx, row in df.iterrows():
            # Get tissue and species, using default if not in dataframe
            tissue = row.get(tissue_col, default_tissue) if tissue_col else default_tissue
            species = row.get(species_col, default_species) if species_col else default_species
            
            # Evaluate this row
            eval_result = self.evaluate_single_celltype(
                predicted_celltype=row[predicted_col],
                gold_standard=row[gold_col],
                tissue=tissue,
                species=species
            )
            
            # Store results
            result_df.at[idx, 'evaluation_score'] = eval_result.get('score', 0)
            result_df.at[idx, 'evaluation_explanation'] = eval_result.get('explanation', '')
        
        return result_df
    
    def _call_llm(self, system_prompt: str, user_prompt: str) -> str:
        """
        Call the LLM using OpenRouter API.
        
        Args:
            system_prompt (str): System prompt for the LLM
            user_prompt (str): User prompt for the LLM
            
        Returns:
            str: LLM response text
        """
        try:
            response = requests.post(
                url="https://openrouter.ai/api/v1/chat/completions",
                headers={
                    "Authorization": f"Bearer {self.api_key}",
                    "HTTP-Referer": "https://localhost:5000",  # Required for OpenRouter
                    "Content-Type": "application/json"
                },
                json={
                    "model": self.model,
                    "temperature": 0,  # Use zero temperature for evaluation
                    "max_tokens": 2000,
                    "messages": [
                        {"role": "system", "content": system_prompt},
                        {"role": "user", "content": user_prompt}
                    ]
                }
            )
            
            # Check if request was successful
            if response.status_code == 200:
                response_data = response.json()
                return response_data['choices'][0]['message']['content']
            else:
                print(f"Error: OpenRouter API returned status code {response.status_code}")
                print(f"Response: {response.text}")
                return ''
                
        except Exception as e:
            print(f"Error making OpenRouter API request: {str(e)}")
            return ''
    
    def _extract_score_json(self, response: str) -> Dict[str, Any]:
        """
        Extract the score and explanation from the LLM response.
        Tries to find and parse a JSON object; falls back to regex if needed.
        
        Args:
            response (str): LLM response text
            
        Returns:
            Dict with score and explanation
        """
        try:
            # First try to find JSON content
            import re
            json_match = re.search(r'\{[\s\S]*?\}', response)
            if json_match:
                json_str = json_match.group(0)
                result = json.loads(json_str)
                
                # Validate that we have the expected fields
                if 'score' in result:
                    return result
            
            # Fallback to regex for score if JSON parsing failed
            score_match = re.search(r'score:?\s*(\d+)', response, re.IGNORECASE)
            if score_match:
                score = int(score_match.group(1))
                
                # Extract explanation if possible
                explanation_match = re.search(r'explanation:?\s*(.+?)(?:\n|$)', response, re.IGNORECASE | re.DOTALL)
                explanation = explanation_match.group(1).strip() if explanation_match else "No explanation provided"
                
                return {
                    "score": score,
                    "explanation": explanation
                }
            
            # If all else fails, return a default response
            return {
                "score": 0,
                "explanation": "Failed to extract score and explanation from LLM response"
            }
            
        except Exception as e:
            print(f"Error extracting score from response: {str(e)}")
            return {
                "score": 0,
                "explanation": f"Error extracting score: {str(e)}"
            }
    
    def _extract_multiple_scores_json(self, response: str) -> Dict[str, Any]:
        """
        Extract multiple scores and explanations from the LLM response.
        
        Args:
            response (str): LLM response text
            
        Returns:
            Dict with individual scores and explanations
        """
        try:
            # First try to find JSON content
            import re
            json_match = re.search(r'\{[\s\S]*?\}', response)
            if json_match:
                json_str = json_match.group(0)
                result = json.loads(json_str)
                
                # Validate that we have the expected fields
                if 'individual_scores' in result:
                    return result
            
            # If JSON parsing failed, return a default response
            return {
                "individual_scores": [],
                "explanations": []
            }
            
        except Exception as e:
            print(f"Error extracting scores from response: {str(e)}")
            return {
                "individual_scores": [],
                "explanations": []
            }


# Example functions to generate simulated data for testing

def generate_simulated_data(n_samples: int = 10) -> pd.DataFrame:
    """
    Generate simulated data for testing the evaluator with an ontology-based scoring system.
    
    Args:
        n_samples (int): Number of samples to generate
        
    Returns:
        pd.DataFrame: Dataframe with simulated data
    """
    species_list = ["human", "mouse", "rat"]
    tissue_list = ["brain", "lung", "liver", "kidney", "heart", "spleen"]
    
    # Gold standard cell types with ontology hierarchy
    # Format: [(gold_standard, general_type)]
    cell_types = [
        ("CD8+ cytotoxic T cell", "T cell"),
        ("CD4+ helper T cell", "T cell"),
        ("CD4+ regulatory T cell", "T cell"),
        ("Memory B cell", "B cell"),
        ("Plasma cell", "B cell"),
        ("M1 macrophage", "Macrophage"),
        ("M2 macrophage", "Macrophage"),
        ("Alveolar macrophage", "Macrophage"),
        ("Conventional dendritic cell", "Dendritic cell"),
        ("Plasmacytoid dendritic cell", "Dendritic cell"),
        ("Neutrophil", "Granulocyte"),
        ("Eosinophil", "Granulocyte"),
        ("NK cell", "Lymphoid cell"),
        ("Vascular endothelial cell", "Endothelial cell"),
        ("Lymphatic endothelial cell", "Endothelial cell"),
        ("Type I pneumocyte", "Epithelial cell"),
        ("Type II pneumocyte", "Epithelial cell"),
        ("Ciliated epithelial cell", "Epithelial cell"),
        ("Fibroblast", "Stromal cell"),
        ("Myofibroblast", "Stromal cell")
    ]
    
    # Map cell types to their lineage for score=3 examples
    lineage_map = {
        "T cell": ["NK cell", "B cell"],  # Close immune cell types (same lymphoid lineage)
        "B cell": ["T cell", "NK cell"],
        "Macrophage": ["Dendritic cell", "Monocyte"],  # Close myeloid cells
        "Dendritic cell": ["Macrophage", "Monocyte"],
        "Granulocyte": ["Monocyte", "Macrophage"],
        "Lymphoid cell": ["T cell", "B cell"],
        "Endothelial cell": ["Vascular endothelial cell", "Lymphatic endothelial cell"],
        "Epithelial cell": ["Type I pneumocyte", "Type II pneumocyte", "Ciliated epithelial cell"],
        "Stromal cell": ["Fibroblast", "Myofibroblast"]
    }
    
    # Generate data
    data = []
    for i in range(n_samples):
        # Select random gold standard
        gold_standard, general_type = cell_types[i % len(cell_types)]
        
        # Generate simulated prediction with varying accuracy based on the 0-5 scale
        accuracy_level = np.random.choice([
            "exact",     # Score 5: Exact match
            "subtle",    # Score 4: Subtle differences
            "close",     # Score 3: Close on ontology tree
            "general",   # Score 2: Only general type correct
            "distant",   # Score 1: Incorrect but related
            "nonsense"   # Score 0: Completely irrelevant
        ], p=[0.2, 0.2, 0.2, 0.2, 0.1, 0.1])
        
        # Generate prediction based on accuracy level
        if accuracy_level == "exact":
            # Score 5: Exact match
            pred = gold_standard
            score = 5
            
        elif accuracy_level == "subtle":
            # Score 4: Subtle differences in activation state or minor variations
            if "T cell" in gold_standard:
                if "regulatory" in gold_standard.lower():
                    pred = gold_standard.replace("regulatory", "suppressor")
                elif "cytotoxic" in gold_standard.lower():
                    pred = gold_standard.replace("cytotoxic", "effector")
                elif "helper" in gold_standard.lower():
                    pred = gold_standard.replace("helper", "CD4+")
                else:
                    pred = gold_standard + " (resting)"
            elif "macrophage" in gold_standard.lower():
                if "M1" in gold_standard:
                    pred = "Inflammatory macrophage"
                elif "M2" in gold_standard:
                    pred = "Tissue-resident macrophage"
                else:
                    pred = gold_standard + " (activated)"
            elif "endothelial" in gold_standard.lower():
                if "vascular" in gold_standard.lower():
                    pred = "Blood vessel endothelial cell"
                elif "lymphatic" in gold_standard.lower():
                    pred = "Lymph vessel endothelial cell"
                else:
                    pred = gold_standard + " (fenestrated)"
            else:
                # Add a minor state variation
                states = ["activated", "mature", "resting"]
                pred = np.random.choice(states) + " " + gold_standard
            score = 4
            
        elif accuracy_level == "close":
            # Score 3: Close on ontology tree (same lineage, different subtype)
            if general_type in lineage_map:
                similar_types = lineage_map[general_type]
                if similar_types:
                    pred = np.random.choice(similar_types)
                else:
                    pred = general_type + " (subtype unclear)"
            else:
                # Fallback if no specific mapping
                pred = general_type + " (variant)"
            score = 3
            
        elif accuracy_level == "general":
            # Score 2: Only general type correct, subtype far off
            if "T cell" in gold_standard:
                pred = "Myeloid cell"  # Far from T cell but still immune
            elif "macrophage" in gold_standard.lower() or "dendritic" in gold_standard.lower():
                pred = "Lymphocyte"  # Far from myeloid but still immune
            elif "endothelial" in gold_standard.lower():
                pred = "Mesenchymal cell"  # Different but still structural
            elif "epithelial" in gold_standard.lower():
                pred = "Connective tissue cell"  # Different but still tissue
            else:
                # Generic distant type in same system
                pred = general_type + "-like cell"
            score = 2
            
        elif accuracy_level == "distant":
            # Score 1: Incorrect general type but somewhat related
            if "T cell" in gold_standard or "B cell" in gold_standard or "NK" in gold_standard:
                pred = "Epithelial cell"  # Completely different lineage
            elif "macrophage" in gold_standard.lower() or "dendritic" in gold_standard.lower():
                pred = "Fibroblast"  # Different system
            elif "endothelial" in gold_standard.lower() or "epithelial" in gold_standard.lower():
                pred = "Hematopoietic cell"  # Different system
            else:
                # Generic distantly related cell
                distant_types = ["Adipocyte", "Melanocyte", "Chondrocyte", "Erythrocyte"]
                pred = np.random.choice(distant_types)
            score = 1
            
        else:  # nonsense
            # Score 0: Completely irrelevant
            nonsense_terms = [
                "Unknown structure", 
                "Cellular debris",
                "Extracellular matrix component",
                "Non-cellular feature",
                "Artifact",
                "In vitro cell line",
                "Undifferentiated precursor"
            ]
            pred = np.random.choice(nonsense_terms)
            score = 0
        
        # Create the data entry
        data.append({
            "species": np.random.choice(species_list),
            "tissue": np.random.choice(tissue_list),
            "gold_standard": gold_standard,
            "predicted_celltype": pred,
            "true_accuracy": score
        })
    
    return pd.DataFrame(data)

def generate_multiple_celltype_samples(n_samples: int = 5, n_types_per_sample: int = 3) -> List[Dict[str, Any]]:
    """
    Generate samples with multiple cell types for testing with ontology-based scoring.
    
    Args:
        n_samples (int): Number of samples to generate
        n_types_per_sample (int): Number of cell types per sample
        
    Returns:
        List[Dict]: List of samples with predicted and gold standard cell types
    """
    species_list = ["human", "mouse"]
    tissue_list = ["brain", "lung", "liver", "kidney"]
    
    # Gold standard cell types with ontology hierarchy
    # Format: [(gold_standard, general_type, lineage)]
    cell_types = [
        # Lymphoid cells
        ("CD8+ cytotoxic T cell", "T cell", "Lymphoid"),
        ("CD4+ helper T cell", "T cell", "Lymphoid"),
        ("CD4+ regulatory T cell", "T cell", "Lymphoid"),
        ("Memory B cell", "B cell", "Lymphoid"),
        ("Plasma cell", "B cell", "Lymphoid"),
        ("CD56+ NK cell", "NK cell", "Lymphoid"),
        ("CD56bright NK cell", "NK cell", "Lymphoid"),
        
        # Myeloid cells
        ("M1 macrophage", "Macrophage", "Myeloid"),
        ("M2 macrophage", "Macrophage", "Myeloid"),
        ("Alveolar macrophage", "Macrophage", "Myeloid"),
        ("Conventional dendritic cell", "Dendritic cell", "Myeloid"),
        ("Plasmacytoid dendritic cell", "Dendritic cell", "Myeloid"),
        ("Neutrophil", "Granulocyte", "Myeloid"),
        ("Eosinophil", "Granulocyte", "Myeloid"),
        
        # Epithelial cells
        ("Type I pneumocyte", "Pneumocyte", "Epithelial"),
        ("Type II pneumocyte", "Pneumocyte", "Epithelial"),
        ("Ciliated epithelial cell", "Airway epithelial cell", "Epithelial"),
        ("Goblet cell", "Secretory epithelial cell", "Epithelial"),
        
        # Endothelial cells
        ("Vascular endothelial cell", "Endothelial cell", "Endothelial"),
        ("Lymphatic endothelial cell", "Endothelial cell", "Endothelial"),
        
        # Stromal cells
        ("Fibroblast", "Stromal cell", "Mesenchymal"),
        ("Myofibroblast", "Stromal cell", "Mesenchymal"),
        ("Pericyte", "Perivascular cell", "Mesenchymal")
    ]
    
    samples = []
    for i in range(n_samples):
        # Select tissue and species for this sample
        tissue = np.random.choice(tissue_list)
        species = np.random.choice(species_list)
        
        # Select n_types_per_sample random cell types
        selected_cells = np.random.choice(len(cell_types), n_types_per_sample, replace=False)
        
        gold_standards = []
        predicted_types = []
        true_scores = []
        
        for cell_idx in selected_cells:
            # Get the gold standard info
            gold_standard, general_type, lineage = cell_types[cell_idx]
            gold_standards.append(gold_standard)
            
            # Generate prediction with varying accuracy
            accuracy_level = np.random.choice([
                "exact",     # Score 5: Exact match
                "subtle",    # Score 4: Subtle differences
                "close",     # Score 3: Close on ontology tree
                "general",   # Score 2: Only general type correct
                "distant",   # Score 1: Incorrect but related
                "nonsense"   # Score 0: Completely irrelevant
            ], p=[0.2, 0.2, 0.2, 0.2, 0.1, 0.1])
            
            if accuracy_level == "exact":
                # Score 5: Perfect match
                pred = gold_standard
                score = 5
                
            elif accuracy_level == "subtle":
                # Score 4: Subtle differences in activation state or minor variations
                if "T cell" in gold_standard:
                    if "regulatory" in gold_standard.lower():
                        pred = gold_standard.replace("regulatory", "suppressor")
                    elif "cytotoxic" in gold_standard.lower():
                        pred = gold_standard.replace("cytotoxic", "effector")
                    elif "helper" in gold_standard.lower():
                        pred = gold_standard.replace("helper", "CD4+")
                    else:
                        pred = gold_standard + " (resting)"
                elif "macrophage" in gold_standard.lower():
                    if "M1" in gold_standard:
                        pred = "Inflammatory macrophage"
                    elif "M2" in gold_standard:
                        pred = "Tissue-resident macrophage"
                    elif "alveolar" in gold_standard.lower():
                        pred = "Lung macrophage"
                    else:
                        pred = gold_standard + " (activated)"
                elif "NK cell" in gold_standard:
                    if "CD56bright" in gold_standard:
                        pred = "Cytokine-producing NK cell"
                    else:
                        pred = "Natural killer cell"
                elif "pneumocyte" in gold_standard.lower():
                    if "Type I" in gold_standard:
                        pred = "Squamous alveolar epithelial cell"
                    elif "Type II" in gold_standard:
                        pred = "Surfactant-producing epithelial cell"
                elif "dendritic" in gold_standard.lower():
                    if "conventional" in gold_standard.lower():
                        pred = "cDC"
                    elif "plasmacytoid" in gold_standard.lower():
                        pred = "pDC"
                else:
                    # Add a minor state variation
                    states = ["activated", "mature", "resting"]
                    pred = np.random.choice(states) + " " + gold_standard
                score = 4
                
            elif accuracy_level == "close":
                # Score 3: Close on ontology tree (same lineage, different subtype)
                if lineage == "Lymphoid":
                    if "T cell" in gold_standard:
                        other_options = ["CD8+ T cell", "CD4+ T cell", "Memory T cell", "Naive T cell", "γδ T cell"]
                        pred = np.random.choice([o for o in other_options if o != gold_standard])
                    elif "B cell" in gold_standard:
                        other_options = ["Naive B cell", "Memory B cell", "Germinal center B cell"]
                        pred = np.random.choice([o for o in other_options if o != gold_standard])
                    elif "NK" in gold_standard:
                        pred = "Innate lymphoid cell"
                    else:
                        pred = "Lymphocyte"
                elif lineage == "Myeloid":
                    if "macrophage" in gold_standard.lower():
                        pred = "Dendritic cell"
                    elif "dendritic" in gold_standard.lower():
                        pred = "Macrophage"
                    elif "neutrophil" in gold_standard.lower() or "eosinophil" in gold_standard.lower():
                        pred = "Granulocyte"
                    else:
                        pred = "Myeloid cell"
                elif lineage == "Epithelial":
                    if "pneumocyte" in gold_standard.lower():
                        pred = "Alveolar epithelial cell"
                    else:
                        pred = "Epithelial cell"
                elif lineage == "Endothelial":
                    if "vascular" in gold_standard.lower():
                        pred = "Endothelial cell"
                    elif "lymphatic" in gold_standard.lower():
                        pred = "Endothelial cell"
                    else:
                        pred = "Endothelial cell"
                elif lineage == "Mesenchymal":
                    if "fibroblast" in gold_standard.lower():
                        pred = "Stromal cell"
                    else:
                        pred = "Mesenchymal cell"
                else:
                    # Generic close type
                    pred = general_type
                score = 3
                
            elif accuracy_level == "general":
                # Score 2: Only general type correct, subtype far off
                if lineage == "Lymphoid":
                    if "T cell" in gold_standard or "B cell" in gold_standard:
                        pred = "Myeloid cell"  # Wrong lineage but still immune
                    else:
                        pred = "Leukocyte"
                elif lineage == "Myeloid":
                    pred = "Lymphocyte"  # Wrong lineage but still immune
                elif lineage == "Epithelial":
                    pred = "Mesenchymal cell"  # Different lineage
                elif lineage == "Endothelial":
                    pred = "Epithelial cell"  # Different but still structural
                elif lineage == "Mesenchymal":
                    pred = "Endothelial cell"  # Different but still structural
                else:
                    # Generic distant type
                    pred = "Cell of " + lineage + " origin"
                score = 2
                
            elif accuracy_level == "distant":
                # Score 1: Incorrect general type but somewhat related
                if lineage in ["Lymphoid", "Myeloid"]:
                    distant_options = ["Epithelial cell", "Fibroblast", "Endothelial cell"]
                    pred = np.random.choice(distant_options)
                elif lineage in ["Epithelial", "Endothelial", "Mesenchymal"]:
                    distant_options = ["Hematopoietic cell", "Immune cell", "Blood cell"]
                    pred = np.random.choice(distant_options)
                else:
                    # Generic distantly related cell
                    distant_types = ["Adipocyte", "Melanocyte", "Chondrocyte", "Erythrocyte"]
                    pred = np.random.choice(distant_types)
                score = 1
                
            else:  # nonsense
                # Score 0: Completely irrelevant
                nonsense_terms = [
                    "Unknown structure", 
                    "Cellular debris",
                    "Extracellular matrix component",
                    "Non-cellular feature",
                    "Artifact",
                    "In vitro cell line",
                    "Undifferentiated precursor"
                ]
                pred = np.random.choice(nonsense_terms)
                score = 0
            
            predicted_types.append(pred)
            true_scores.append(score)
        
        samples.append({
            "tissue": tissue,
            "species": species,
            "gold_standards": gold_standards,
            "predicted_celltypes": predicted_types,
            "true_scores": true_scores
        })
    
    return samples


# Command-line interface
def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Evaluate cell type annotations using LLMs on a 0-5 scale')
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Common arguments
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--api-key', type=str, help='OpenRouter API key')
    parent_parser.add_argument('--model', type=str, default='anthropic/claude-3.5-sonnet', 
                              help='Model to use for evaluation')
    
    # Single cell type evaluation
    single_parser = subparsers.add_parser('single', parents=[parent_parser], 
                                       help='Evaluate a single cell type annotation (0-5 scale)')
    single_parser.add_argument('--predicted', type=str, required=True, 
                            help='Predicted cell type')
    single_parser.add_argument('--gold', type=str, required=True, 
                             help='Gold standard annotation')
    single_parser.add_argument('--tissue', type=str, default='unknown', 
                             help='Tissue context')
    single_parser.add_argument('--species', type=str, default='human', 
                             help='Species context')
    
    # Multiple cell types evaluation
    multi_parser = subparsers.add_parser('multiple', parents=[parent_parser], 
                                      help='Evaluate multiple cell type annotations (0-5 scale)')
    multi_parser.add_argument('--predicted', type=str, required=True, nargs='+', 
                           help='Predicted cell types (space-separated)')
    multi_parser.add_argument('--gold', type=str, required=True, nargs='+', 
                           help='Gold standard annotations (space-separated)')
    multi_parser.add_argument('--tissue', type=str, default='unknown', 
                           help='Tissue context')
    multi_parser.add_argument('--species', type=str, default='human', 
                           help='Species context')
    
    # Batch evaluation from CSV
    batch_parser = subparsers.add_parser('batch', parents=[parent_parser], 
                                      help='Evaluate batch from CSV file (0-5 scale)')
    batch_parser.add_argument('--input', type=str, required=True, 
                           help='Input CSV file path')
    batch_parser.add_argument('--output', type=str, required=True, 
                           help='Output CSV file path')
    batch_parser.add_argument('--predicted-col', type=str, required=True, 
                           help='Column name for predictions')
    batch_parser.add_argument('--gold-col', type=str, required=True, 
                           help='Column name for gold standards')
    batch_parser.add_argument('--tissue-col', type=str, 
                           help='Column name for tissue context')
    batch_parser.add_argument('--species-col', type=str, 
                           help='Column name for species context')
    batch_parser.add_argument('--default-tissue', type=str, default='unknown', 
                           help='Default tissue if not in CSV')
    batch_parser.add_argument('--default-species', type=str, default='human', 
                           help='Default species if not in CSV')
    
    # Simulate data
    sim_parser = subparsers.add_parser('simulate', 
                                     help='Generate simulated data for testing (with 0-5 score scale)')
    sim_parser.add_argument('--output', type=str, required=True, 
                          help='Output CSV file path')
    sim_parser.add_argument('--samples', type=int, default=10, 
                          help='Number of samples to generate')
    
    args = parser.parse_args()
    
    if args.command == 'single':
        evaluator = LLMEvaluator(api_key=args.api_key, model=args.model)
        result = evaluator.evaluate_single_celltype(
            predicted_celltype=args.predicted,
            gold_standard=args.gold,
            tissue=args.tissue,
            species=args.species
        )
        print(json.dumps(result, indent=2))
        
    elif args.command == 'multiple':
        if len(args.predicted) != len(args.gold):
            print("Error: Number of predicted cell types must match number of gold standards")
            return
        
        evaluator = LLMEvaluator(api_key=args.api_key, model=args.model)
        result = evaluator.evaluate_multiple_celltypes(
            predicted_celltypes=args.predicted,
            gold_standards=args.gold,
            tissue=args.tissue,
            species=args.species
        )
        print(json.dumps(result, indent=2))
        
    elif args.command == 'batch':
        try:
            df = pd.read_csv(args.input)
            evaluator = LLMEvaluator(api_key=args.api_key, model=args.model)
            result_df = evaluator.batch_evaluate_from_dataframe(
                df=df,
                predicted_col=args.predicted_col,
                gold_col=args.gold_col,
                tissue_col=args.tissue_col,
                species_col=args.species_col,
                default_tissue=args.default_tissue,
                default_species=args.default_species
            )
            result_df.to_csv(args.output, index=False)
            print(f"Batch evaluation completed. Results saved to {args.output}")
            
        except Exception as e:
            print(f"Error processing batch: {str(e)}")
            
    elif args.command == 'simulate':
        simulated_data = generate_simulated_data(args.samples)
        simulated_data.to_csv(args.output, index=False)
        print(f"Generated {args.samples} simulated samples. Saved to {args.output}")
        
    else:
        parser.print_help()

if __name__ == "__main__":
    # If run as a script, use the command-line interface
    main()

def calculate_evaluation_metrics(eval_df: pd.DataFrame, score_col: str = 'evaluation_score', 
                         true_score_col: str = None) -> Dict[str, float]:
    """
    Calculate metrics from batch evaluation results.
    
    Args:
        eval_df (pd.DataFrame): DataFrame with evaluation results
        score_col (str): Column name for evaluation scores (0-5 scale)
        true_score_col (str, optional): Column name for true scores if available
        
    Returns:
        Dict[str, float]: Dictionary with evaluation metrics
    """
    metrics = {
        'mean_score': eval_df[score_col].mean(),
        'median_score': eval_df[score_col].median(),
        'min_score': eval_df[score_col].min(),
        'max_score': eval_df[score_col].max(),
        'std_score': eval_df[score_col].std(),
        'count': len(eval_df),
        'perfect_ratio': (eval_df[score_col] == 5).mean(),
        'very_good_ratio': (eval_df[score_col] == 4).mean(),
        'good_ratio': (eval_df[score_col] == 3).mean(),
        'partial_ratio': (eval_df[score_col] == 2).mean(),
        'poor_ratio': (eval_df[score_col] == 1).mean(),
        'nonsensical_ratio': (eval_df[score_col] == 0).mean(),
    }
    
    # If true scores are available, calculate correlation and error metrics
    if true_score_col and true_score_col in eval_df.columns:
        metrics.update({
            'correlation': eval_df[[score_col, true_score_col]].corr().iloc[0, 1],
            'mae': np.abs(eval_df[score_col] - eval_df[true_score_col]).mean(),
            'rmse': np.sqrt(((eval_df[score_col] - eval_df[true_score_col]) ** 2).mean()),
        })
    
    return metrics
