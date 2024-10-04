from openai import OpenAI
import pandas as pd

def evaluate_annotation(prompt, system_prompt='''You are an professional biologist specialized in single cell analysis. Your task is to analyze the given annotation result and compare it to the true cell type based on your expertise. You will then categorize the result as either correct, partially correct, or incorrect. Include a short one sentence explanation for each evaluation. Note that sometimes celltypes names can be slightly different, it is okay, use your judgement to see if they are actually the same celltype.

## Input Format

You will receive information like below:

True cell type: Cell type name,
main_cell_type: Predicted main cell type,
sub_cell_types: subcelltype1, subcelltype2, subcelltype3
 

## Evaluation steps

1. If the main_cell_type matches the true cell type, consider it fully correct. If the main_cell_type is a broader category that includes the true cell type, look at the sub_cell_types.
2. If the first celltype in the sub_celltypes matche the true cell type, consider it fully correct.
3. If the first celltype in the sub_celltypes is closely related to the true cell type, consider it partially correct.
4. Otherwise, consider it incorrect.

Provide an summary at last

Correct: number of correct
Partially correct: number of partially correct
Incorrect: number of incorrect


Example1

True Cell Type: Epithelial cell (club)

Predicted Main Cell Type: secretory cells
Predicted Sub Cell Types: club cells, goblet cells, multiciliated cells
Evaluation: fully correct: first celltype in sub cell type is club cell matches the true celltype

Example2

True Cell Type:Fibroblast
Predicted Main Cell Type:mesenchymal cells
Predicted Sub Cell Types:fibroblasts, smooth muscle cells, myofibroblasts

Evaluation:
Fully Correct (First sub cell type matches true cell type)

Example3

True Cell Type:Epithelial cell (alveolar type II)
Predicted Main Cell Type:alveolar type II (AT2) cells
Predicted Sub Cell Types:dendritic cells, lung epithelial cells

Evaluation:Fully Correct (The main type is a fully match though first subtype is different, it does not matter, as long as the main type is correct this is fully correct)

''', model="gpt-4o", temperature=0):
    client = OpenAI()

    # Convert prompt to string if it's a DataFrame or Series
    if isinstance(prompt, (pd.DataFrame, pd.Series)):
        prompt = prompt.to_string(index=False)

    try:
        response = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": system_prompt},
                {"role": "user", "content": prompt}
            ],
            temperature=temperature
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"An error occurred: {str(e)}"