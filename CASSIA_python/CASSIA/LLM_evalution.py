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
    
    def __init__(self, api_key: str = None, model: str = "anthropic/claude-3.5-sonnet"):
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
Evaluate how accurate and specific the prediction is compared to the gold standard.

Score the prediction on a scale of 0-100 where:
- 0-20: Completely incorrect with no relation to the correct cell type
- 21-40: Incorrect but has some relation to the correct cell type
- 41-60: Partially correct, identifying the general cell type but missing specificity
- 61-80: Mostly correct but missing some nuance or specificity
- 81-100: Highly accurate match to the gold standard

Your response must include:
1. A numeric score (0-100)
2. A brief explanation of why you assigned this score
3. A JSON-formatted result with the format:
{"score": X, "explanation": "Your explanation here"}
"""
        
        user_prompt = f"""
Evaluate the following cell type annotation prediction against the gold standard:

Context:
- Tissue: {tissue}
- Species: {species}

Predicted cell type: {predicted_celltype}
Gold standard annotation: {gold_standard}

Please provide your score (0-100) and explanation, followed by a JSON in the format:
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
Evaluate how accurate and specific each prediction is compared to its gold standard.

For each pair, score the prediction on a scale of 0-100 where:
- 0-20: Completely incorrect with no relation to the correct cell type
- 21-40: Incorrect but has some relation to the correct cell type
- 41-60: Partially correct, identifying the general cell type but missing specificity
- 61-80: Mostly correct but missing some nuance or specificity
- 81-100: Highly accurate match to the gold standard

Also provide an overall score across all predictions.

Your response must include:
1. Individual scores for each prediction
2. An overall score
3. Brief explanations for your scoring decisions
4. A JSON-formatted result with the format:
{
  "overall_score": X,
  "individual_scores": [score1, score2, ...],
  "explanations": ["explanation1", "explanation2", ...],
  "overall_explanation": "Your overall explanation"
}
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

Please provide your evaluation with individual scores, an overall score, and explanations, followed by a JSON in the specified format.
"""
        
        return system_prompt, user_prompt
    
    def get_hierarchy_aware_prompts(self, 
                                   predicted_celltype: str, 
                                   gold_standard: str, 
                                   tissue: str, 
                                   species: str) -> Tuple[str, str]:
        """
        Generate prompts for hierarchy-aware evaluation that understands cell type ontology.
        
        Args:
            predicted_celltype (str): The cell type predicted by the annotation system
            gold_standard (str): The gold standard annotation
            tissue (str): The tissue context
            species (str): The species context
            
        Returns:
            Tuple[str, str]: System prompt and user prompt
        """
        
        system_prompt = """You are an expert cell biologist with deep understanding of cell type hierarchies and ontologies.
You will evaluate a predicted cell type against a gold standard (correct) annotation, accounting for hierarchical relationships.

When evaluating, consider these principles:
1. Parent-child relationships: If the predicted cell type is a direct parent of the gold standard (more general but correct lineage), this is better than an unrelated cell type.
2. Sibling relationships: If the predicted cell type is a sibling of the gold standard (same parent), this is better than a completely unrelated cell type.
3. Specificity level: Annotations at the correct level of specificity should be rewarded.

Score the prediction on a scale of 0-100 where:
- 0-20: Completely incorrect with no relation to the correct cell type
- 21-40: Wrong cell type but in related lineage or with similar function
- 41-60: Correct general cell type (parent) but missing specificity
- 61-80: Mostly correct, possibly a sibling cell type or missing minor specificity
- 81-100: Highly accurate match, identifying the correct specific cell type

Your response must include:
1. A numeric score (0-100)
2. A brief explanation of the hierarchical relationship (if any) between the prediction and gold standard
3. A JSON-formatted result with the format:
{"score": X, "explanation": "Your explanation", "relationship": "[exact/parent/child/sibling/unrelated]"}
"""
        
        user_prompt = f"""
Evaluate the following cell type annotation prediction against the gold standard, considering hierarchical relationships:

Context:
- Tissue: {tissue}
- Species: {species}

Predicted cell type: {predicted_celltype}
Gold standard annotation: {gold_standard}

Please provide your score (0-100) and explanation, followed by a JSON in the format:
{{
  "score": X,
  "explanation": "Your explanation",
  "relationship": "[exact/parent/child/sibling/unrelated]"
}}
"""
        
        return system_prompt, user_prompt
    
    def evaluate_single_celltype(self, 
                                predicted_celltype: str, 
                                gold_standard: str, 
                                tissue: str, 
                                species: str, 
                                hierarchy_aware: bool = False) -> Dict[str, Any]:
        """
        Evaluate a single cell type annotation against a gold standard.
        
        Args:
            predicted_celltype (str): The cell type predicted by the annotation system
            gold_standard (str): The gold standard annotation
            tissue (str): The tissue context
            species (str): The species context
            hierarchy_aware (bool): Whether to use hierarchy-aware evaluation
            
        Returns:
            Dict containing score and explanation
        """
        
        if hierarchy_aware:
            system_prompt, user_prompt = self.get_hierarchy_aware_prompts(
                predicted_celltype, gold_standard, tissue, species
            )
        else:
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
            Dict containing overall score, individual scores, and explanation
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
                                     default_species: str = "human",
                                     hierarchy_aware: bool = False) -> pd.DataFrame:
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
            hierarchy_aware (bool): Whether to use hierarchy-aware evaluation
            
        Returns:
            pd.DataFrame with original data plus evaluation results
        """
        
        # Create copies of the dataframe columns for modification
        result_df = df.copy()
        result_df['evaluation_score'] = None
        result_df['evaluation_explanation'] = None
        if hierarchy_aware:
            result_df['relationship'] = None
        
        for idx, row in df.iterrows():
            # Get tissue and species, using default if not in dataframe
            tissue = row.get(tissue_col, default_tissue) if tissue_col else default_tissue
            species = row.get(species_col, default_species) if species_col else default_species
            
            # Evaluate this row
            eval_result = self.evaluate_single_celltype(
                predicted_celltype=row[predicted_col],
                gold_standard=row[gold_col],
                tissue=tissue,
                species=species,
                hierarchy_aware=hierarchy_aware
            )
            
            # Store results
            result_df.at[idx, 'evaluation_score'] = eval_result.get('score', 0)
            result_df.at[idx, 'evaluation_explanation'] = eval_result.get('explanation', '')
            if hierarchy_aware and 'relationship' in eval_result:
                result_df.at[idx, 'relationship'] = eval_result.get('relationship', 'unknown')
        
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
            Dict with overall score, individual scores, and explanations
        """
        try:
            # First try to find JSON content
            import re
            json_match = re.search(r'\{[\s\S]*?\}', response)
            if json_match:
                json_str = json_match.group(0)
                result = json.loads(json_str)
                
                # Validate that we have the expected fields
                if 'overall_score' in result and 'individual_scores' in result:
                    return result
            
            # If JSON parsing failed, return a default response
            return {
                "overall_score": 0,
                "individual_scores": [],
                "explanations": [],
                "overall_explanation": "Failed to extract scores and explanations from LLM response"
            }
            
        except Exception as e:
            print(f"Error extracting scores from response: {str(e)}")
            return {
                "overall_score": 0,
                "individual_scores": [],
                "explanations": [],
                "overall_explanation": f"Error extracting scores: {str(e)}"
            }


# Example functions to generate simulated data for testing

def generate_simulated_data(n_samples: int = 10) -> pd.DataFrame:
    """
    Generate simulated data for testing the evaluator.
    
    Args:
        n_samples (int): Number of samples to generate
        
    Returns:
        pd.DataFrame: Dataframe with simulated data
    """
    species_list = ["human", "mouse", "rat"]
    tissue_list = ["brain", "lung", "liver", "kidney", "heart", "spleen"]
    
    # Gold standard cell types
    gold_celltypes = [
        "T cell",
        "B cell",
        "Neutrophil",
        "Macrophage",
        "Dendritic cell",
        "Natural killer cell",
        "Monocyte",
        "Fibroblast",
        "Epithelial cell",
        "Endothelial cell"
    ]
    
    # Generate data
    data = []
    for i in range(n_samples):
        # Select random gold standard
        gold = gold_celltypes[i % len(gold_celltypes)]
        
        # Generate simulated prediction with varying accuracy
        accuracy_level = np.random.choice([
            "correct",  # Exactly correct
            "partial",  # Partially correct (general type but not specific)
            "related",  # Related but incorrect
            "wrong"     # Completely wrong
        ], p=[0.4, 0.3, 0.2, 0.1])
        
        # Generate prediction based on accuracy level
        if accuracy_level == "correct":
            pred = gold
        elif accuracy_level == "partial":
            if "cell" in gold.lower():
                pred = gold.split()[0] + " cell"  # Just the first part
            else:
                pred = "Immature " + gold
        elif accuracy_level == "related":
            # Map to a related but different cell type
            related_map = {
                "T cell": "Lymphocyte",
                "B cell": "Plasma cell",
                "Neutrophil": "Granulocyte",
                "Macrophage": "Monocyte",
                "Dendritic cell": "Antigen-presenting cell",
                "Natural killer cell": "Cytotoxic lymphocyte",
                "Monocyte": "Myeloid cell",
                "Fibroblast": "Stromal cell",
                "Epithelial cell": "Basal cell",
                "Endothelial cell": "Vascular cell"
            }
            pred = related_map.get(gold, "Unknown cell")
        else:  # wrong
            # Choose a completely different cell type
            other_celltypes = [ct for ct in gold_celltypes if ct != gold]
            pred = np.random.choice(other_celltypes)
        
        # Add some random variation to predictions
        if np.random.random() < 0.3 and accuracy_level != "correct":
            modifiers = ["activated", "mature", "immature", "proliferating", "resting"]
            pred = np.random.choice(modifiers) + " " + pred
        
        # Create the data entry
        data.append({
            "species": np.random.choice(species_list),
            "tissue": np.random.choice(tissue_list),
            "gold_standard": gold,
            "predicted_celltype": pred,
            "true_accuracy": {
                "correct": 90 + np.random.randint(0, 11),
                "partial": 50 + np.random.randint(0, 31),
                "related": 20 + np.random.randint(0, 31),
                "wrong": np.random.randint(0, 21)
            }[accuracy_level]
        })
    
    return pd.DataFrame(data)

def generate_multiple_celltype_samples(n_samples: int = 5, n_types_per_sample: int = 3) -> List[Dict[str, Any]]:
    """
    Generate samples with multiple cell types for testing.
    
    Args:
        n_samples (int): Number of samples to generate
        n_types_per_sample (int): Number of cell types per sample
        
    Returns:
        List[Dict]: List of samples with predicted and gold standard cell types
    """
    species_list = ["human", "mouse"]
    tissue_list = ["brain", "lung", "liver", "kidney"]
    
    # Base cell types
    base_celltypes = [
        "T cell",
        "B cell",
        "Neutrophil",
        "Macrophage",
        "Dendritic cell",
        "Natural killer cell",
        "Monocyte",
        "Fibroblast",
        "Epithelial cell",
        "Endothelial cell"
    ]
    
    # Specific subtypes for each base cell type
    subtype_map = {
        "T cell": ["CD4+ T cell", "CD8+ T cell", "Regulatory T cell", "Memory T cell"],
        "B cell": ["Naive B cell", "Memory B cell", "Plasma cell", "Germinal center B cell"],
        "Neutrophil": ["Immature neutrophil", "Segmented neutrophil", "Activated neutrophil"],
        "Macrophage": ["M1 macrophage", "M2 macrophage", "Alveolar macrophage", "Tumor-associated macrophage"],
        "Dendritic cell": ["Conventional DC", "Plasmacytoid DC", "Follicular DC", "Langerhans cell"],
        "Natural killer cell": ["CD56bright NK cell", "CD56dim NK cell", "Cytotoxic NK cell"],
        "Monocyte": ["Classical monocyte", "Non-classical monocyte", "Intermediate monocyte"],
        "Fibroblast": ["Activated fibroblast", "Myofibroblast", "Cardiac fibroblast"],
        "Epithelial cell": ["Basal epithelial cell", "Ciliated epithelial cell", "Goblet cell", "Club cell"],
        "Endothelial cell": ["Vascular endothelial cell", "Lymphatic endothelial cell", "Sinusoidal endothelial cell"]
    }
    
    samples = []
    for i in range(n_samples):
        # Select tissue and species for this sample
        tissue = np.random.choice(tissue_list)
        species = np.random.choice(species_list)
        
        # Select n_types_per_sample random base cell types
        selected_base_types = np.random.choice(base_celltypes, n_types_per_sample, replace=False)
        
        gold_standards = []
        predicted_types = []
        
        for base_type in selected_base_types:
            # Select a subtype as the gold standard
            subtypes = subtype_map[base_type]
            gold_standard = np.random.choice(subtypes)
            gold_standards.append(gold_standard)
            
            # Generate prediction with varying accuracy
            accuracy_level = np.random.choice([
                "exact",   # Exactly correct
                "general", # General type only
                "related", # Related but incorrect
                "wrong"    # Completely wrong
            ], p=[0.4, 0.3, 0.2, 0.1])
            
            if accuracy_level == "exact":
                pred = gold_standard
            elif accuracy_level == "general":
                pred = base_type
            elif accuracy_level == "related":
                # Choose a different subtype of the same base type
                other_subtypes = [st for st in subtypes if st != gold_standard]
                pred = np.random.choice(other_subtypes) if other_subtypes else base_type
            else:  # wrong
                # Choose a completely different cell type
                other_bases = [bt for bt in base_celltypes if bt != base_type]
                wrong_base = np.random.choice(other_bases)
                wrong_subtypes = subtype_map[wrong_base]
                pred = np.random.choice(wrong_subtypes)
            
            predicted_types.append(pred)
        
        samples.append({
            "tissue": tissue,
            "species": species,
            "gold_standards": gold_standards,
            "predicted_celltypes": predicted_types
        })
    
    return samples


# Command-line interface
def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Evaluate cell type annotations using LLMs')
    subparsers = parser.add_subparsers(dest='command', help='Command to run')
    
    # Common arguments
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--api-key', type=str, help='OpenRouter API key')
    parent_parser.add_argument('--model', type=str, default='anthropic/claude-3.5-sonnet', 
                              help='Model to use for evaluation')
    
    # Single cell type evaluation
    single_parser = subparsers.add_parser('single', parents=[parent_parser], 
                                       help='Evaluate a single cell type annotation')
    single_parser.add_argument('--predicted', type=str, required=True, 
                            help='Predicted cell type')
    single_parser.add_argument('--gold', type=str, required=True, 
                             help='Gold standard annotation')
    single_parser.add_argument('--tissue', type=str, default='unknown', 
                             help='Tissue context')
    single_parser.add_argument('--species', type=str, default='human', 
                             help='Species context')
    single_parser.add_argument('--hierarchy', action='store_true', 
                             help='Use hierarchy-aware evaluation')
    
    # Multiple cell types evaluation
    multi_parser = subparsers.add_parser('multiple', parents=[parent_parser], 
                                      help='Evaluate multiple cell type annotations')
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
                                      help='Evaluate batch from CSV file')
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
    batch_parser.add_argument('--hierarchy', action='store_true', 
                           help='Use hierarchy-aware evaluation')
    
    # Simulate data
    sim_parser = subparsers.add_parser('simulate', 
                                     help='Generate simulated data for testing')
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
            species=args.species,
            hierarchy_aware=args.hierarchy
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
                default_species=args.default_species,
                hierarchy_aware=args.hierarchy
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
        score_col (str): Column name for evaluation scores
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
        'excellent_ratio': (eval_df[score_col] >= 80).mean(),
        'good_ratio': ((eval_df[score_col] >= 60) & (eval_df[score_col] < 80)).mean(),
        'fair_ratio': ((eval_df[score_col] >= 40) & (eval_df[score_col] < 60)).mean(),
        'poor_ratio': ((eval_df[score_col] >= 20) & (eval_df[score_col] < 40)).mean(),
        'bad_ratio': (eval_df[score_col] < 20).mean(),
    }
    
    # If true scores are available, calculate correlation and error metrics
    if true_score_col and true_score_col in eval_df.columns:
        metrics.update({
            'correlation': eval_df[[score_col, true_score_col]].corr().iloc[0, 1],
            'mae': np.abs(eval_df[score_col] - eval_df[true_score_col]).mean(),
            'rmse': np.sqrt(((eval_df[score_col] - eval_df[true_score_col]) ** 2).mean()),
        })
    
    return metrics
