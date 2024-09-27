def format_annotation_output(raw_output):
    formatted = "Final Annotation Agent:\n"
    
    # Split the raw output into sections
    sections = raw_output.split("Step ")
    
    # Format each step
    for i, section in enumerate(sections[1:], 1):
        lines = section.strip().split("\n")
        formatted += f"\nStep {i}: {lines[0]}\n"
        for line in lines[1:]:
            formatted += f"  {line}\n"
    
    # Format the Coupling Validator section if present
    if "Coupling Validator:" in raw_output:
        validator_sections = raw_output.split("Coupling Validator:")
        if len(validator_sections) > 1:
            validator_section = validator_sections[1].strip()
            formatted += "\nCoupling Validator:\n"
            for line in validator_section.split("\n"):
                formatted += f"  {line}\n"
    
    return formatted

# Example usage
raw_output = """Final Annotation Agent: Step 1: List the Key Functional Markers The key functional markers from the provided list are: 1. CCL21: This is a chemokine that plays a role in the immune response by attracting lymphocytes and dendritic cells towards the sites of inflammation. 2. MMRN1: Multimerin 1 is a protein that is involved in platelet formation and coagulation. 3. PKHD1L1: Polycystic kidney and hepatic disease 1-like 1 protein is involved in the development of kidneys and the liver. 4. PROX1: Prospero homeobox 1 is a transcription factor involved in the development of the lymphatic system. 5. TBX1: T-box 1 is a transcription factor that plays a crucial role in the development of the heart and the pharyngeal system. 6. LRCOL1: Leucine rich colipase-like 1 is a protein with unknown function. 7. GPR182: G protein-coupled receptor 182 is involved in various cellular processes including cell proliferation and survival. 8. TFF3: Trefoil factor 3 is involved in mucosal healing and has been associated with various cancers. 9. FLT4: Fms related tyrosine kinase 4 is a receptor for vascular endothelial growth factor C (VEGF-C), and it's involved in lymphangiogenesis. 10. HOXD3: Homeobox D3 is a transcription factor involved in the regulation of body plan development. Step 2: List the Key Cell Type Markers The key cell type markers from the provided list are: 1. CCL21: Expressed in lymphoid tissues, endothelial cells of blood vessels, and fibroblastic reticular cells. 2. MMRN1: Expressed in endothelial cells. 3. PKHD1L1: Expressed in renal and hepatic cells. 4. PROX1: Expressed in lymphatic endothelial cells. 5. TBX1: Expressed in cells of the pharyngeal system and heart. 6. LRCOL1: Expression pattern is not well defined. 7. GPR182: Expressed in endothelial cells. 8. TFF3: Expressed in goblet cells of the intestine and respiratory tract. 9. FLT4: Expressed in lymphatic endothelial cells. 10. HOXD3: Expressed in cells involved in body plan development. Step 3: Cross-reference Known Databases Cross-referencing these markers with the Human Cell Atlas and relevant literature, it is evident that these markers are associated with endothelial cells, particularly lymphatic endothelial cells, and cells involved in organ development and function. Step 4: Determine the Most Probable General Cell Type Based on the expression of these markers, the most likely general cell type of the cluster is endothelial cells. Step 5: Identify the Top 3 Most Probable Sub Cell Types Based on the expression of these markers, the top three most probable sub cell types within the endothelial cells are: 1. Lymphatic endothelial cells (PROX1, FLT4) 2. Blood vessel endothelial cells (CCL21, MMRN1, GPR182) 3. Organ-specific endothelial cells (PKHD1L1, TBX1) The most likely subtype is lymphatic endothelial cells due to the expression of PROX1 and FLT4, which are specific markers for these cells. Step 6: Identify the Most Probable Sub-Sub Cell Type Within the lymphatic endothelial cells, the most specific cell type cannot be determined based on the provided markers. Step 7: Provide a Concise Summary of Your Analysis The provided markers indicate that the cell cluster is most likely composed of endothelial cells, specifically lymphatic endothelial cells. This is based on the expression of PROX1 and FLT4, which are specific markers for lymphatic endothelial cells. Other markers suggest the presence of blood vessel endothelial cells and organ-specific endothelial cells, but the strongest evidence points towards lymphatic endothelial cells. FINAL ANNOTATION COMPLETED.

Coupling Validator: Validation result: VALIDATION PASSED Confidence score: 8 Brief justification: The annotation is consistent with the provided marker list and the identified cell type (endothelial cells, specifically lymphatic endothelial cells) is plausible given the tissue type (lung). The lung contains a network of blood and lymphatic vessels, which would contain endothelial cells. The markers PROX1 and FLT4 are indeed associated with lymphatic endothelial cells. However, the presence of other markers like TBX1 (associated with heart and pharyngeal system) and PKHD1L1 (associated with renal and hepatic cells) suggest that there might be some mixed cell types present, which could be further investigated."""

print(format_annotation_output(raw_output))