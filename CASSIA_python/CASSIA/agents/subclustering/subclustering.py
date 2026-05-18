try:
    from CASSIA.engine.tools_function import *
    from CASSIA.core.model_settings import get_agent_default
except ImportError:
    try:
        from ...engine.tools_function import *
        from ...core.model_settings import get_agent_default
    except ImportError:
        from tools_function import *
        from model_settings import get_agent_default

try:
    from CASSIA.core.llm_utils import *
except ImportError:
    try:
        from ...core.llm_utils import *
    except ImportError:
        from llm_utils import *

import pandas as pd
import re

try:
    from CASSIA.agents.reference_agent import ReferenceAgent
except ImportError:
    try:
        from ...agents.reference_agent import ReferenceAgent
    except ImportError:
        try:
            from reference_agent import ReferenceAgent
        except ImportError:
            ReferenceAgent = None


def _get_get_top_markers():
    """Lazy import of get_top_markers to avoid circular imports."""
    try:
        from CASSIA.core.marker_utils import get_top_markers
        return get_top_markers
    except ImportError:
        try:
            from ...core.marker_utils import get_top_markers
            return get_top_markers
        except ImportError:
            from marker_utils import get_top_markers
            return get_top_markers


def _prepare_subcluster_marker_dataframe(marker, n_genes=50):
    """Return a two-column marker dataframe suitable for subcluster prompts."""
    if isinstance(marker, pd.DataFrame):
        marker_df = marker.copy()
    elif isinstance(marker, str):
        marker_df = pd.read_csv(marker)
    else:
        raise ValueError("marker must be a pandas DataFrame or a CSV file path")

    if len(marker_df.columns) > 2:
        get_top_markers = _get_get_top_markers()
        marker_df = get_top_markers(marker_df, n_genes=n_genes)

    return marker_df


def _parse_marker_values(marker_value, n_genes=50):
    """Parse a marker-list cell into ordered marker symbols."""
    if marker_value is None:
        return []

    if isinstance(marker_value, (list, tuple, set)):
        raw_markers = list(marker_value)
    else:
        marker_text = str(marker_value)
        raw_markers = re.split(r"[,;|\n]+", marker_text)
        if len(raw_markers) <= 1:
            raw_markers = re.split(r"\s+", marker_text)

    markers = []
    seen = set()
    for marker in raw_markers:
        marker = str(marker).strip().strip("'\"`")
        if not marker:
            continue
        marker_key = marker.upper()
        if marker_key in seen:
            continue
        seen.add(marker_key)
        markers.append(marker)
        if n_genes and len(markers) >= n_genes:
            break
    return markers


def _combine_context(additional_context, reference_context):
    if additional_context and reference_context:
        return f"{additional_context}\n\n{reference_context}"
    return additional_context or reference_context


def build_subcluster_reference_context(
    marker,
    major_cluster_info,
    provider="openrouter",
    n_genes=50,
    tissue=None,
    species=None,
    reference_provider=None,
    reference_model=None,
    reference_cell_type_hint=None,
    reference_depth="detailed",
    reference_max_content_length=5000,
    reference_max_context_length=12000,
    verbose=False,
):
    """
    Build an expert-reference context block for a subclustering run.

    The function runs reference retrieval per subcluster marker set, deduplicates
    selected documents, and returns a single text block that can be appended to
    the subclustering prompt.
    """
    info = {
        "reference_used": False,
        "references_used": [],
        "clusters": [],
        "reason": "",
    }

    if ReferenceAgent is None:
        info["reason"] = "Reference agent not available"
        return "", info

    marker_df = _prepare_subcluster_marker_dataframe(marker, n_genes=n_genes)
    ref_provider = reference_provider or provider
    cell_type_hint = reference_cell_type_hint or major_cluster_info
    agent = ReferenceAgent(provider=ref_provider, model=reference_model)

    marker_sets = []
    for _, row in marker_df.iterrows():
        cluster_id = str(row.iloc[0])
        markers = _parse_marker_values(row.iloc[1], n_genes=n_genes)
        if markers:
            marker_sets.append({"cluster_id": cluster_id, "markers": markers})

    if hasattr(agent, "get_reference_brief_for_subclusters"):
        try:
            brief_result = agent.get_reference_brief_for_subclusters(
                marker_sets=marker_sets,
                major_cluster_info=major_cluster_info,
                tissue=tissue,
                species=species,
                cell_type_hint=cell_type_hint,
                depth=reference_depth,
                max_reference_content_length=reference_max_content_length,
                max_brief_length=reference_max_context_length,
            )
        except Exception as exc:
            brief_result = {
                "should_use_reference": False,
                "content": "",
                "references_used": [],
                "reasoning": f"Agentic reference brief failed: {exc}",
            }

        info["references_used"] = brief_result.get("references_used", []) or []
        info["clusters"] = brief_result.get("planning", {}).get("cluster_hypotheses", [])
        info["tool_trace"] = brief_result.get("tool_trace", [])
        info["reason"] = brief_result.get("reasoning", "")

        if brief_result.get("should_use_reference") and brief_result.get("content"):
            info["reference_used"] = True
            context = (
                "<expert_reference>\n"
                "Agent-generated subtype reference brief for this subclustering run. "
                "Use it as literature-grounded guidance, but prioritize the observed "
                "marker genes and parent-cluster context when they conflict. Prefer "
                "consensus marker-program labels as the primary subtype; keep paper-specific "
                "labels only as traceability evidence when the marker support is strong.\n\n"
                f"Parent cluster context: {major_cluster_info}\n\n"
                f"{brief_result['content']}\n"
                "</expert_reference>"
            )
            if verbose:
                print(f"Reference context added for subclustering: {', '.join(info['references_used'])}")
            return context, info

    seen_references = set()
    cluster_blocks = []

    for marker_set in marker_sets:
        cluster_id = marker_set["cluster_id"]
        markers = marker_set["markers"]
        try:
            ref_result = agent.get_reference_for_markers(
                markers=markers[:20],
                tissue=tissue,
                species=species,
                cell_type_hint=cell_type_hint,
                depth=reference_depth,
                max_content_length=reference_max_content_length,
            )
        except Exception as exc:
            info["clusters"].append({
                "cluster_id": cluster_id,
                "reference_used": False,
                "references_used": [],
                "reason": str(exc),
            })
            continue

        references_used = ref_result.get("references_used", []) or []
        cluster_info = {
            "cluster_id": cluster_id,
            "reference_used": bool(ref_result.get("should_use_reference") and references_used),
            "preliminary_cell_type": ref_result.get("preliminary_cell_type"),
            "cell_type_range": ref_result.get("cell_type_range", []),
            "references_used": references_used,
            "reason": ref_result.get("reasoning", ""),
        }
        info["clusters"].append(cluster_info)

        if not ref_result.get("should_use_reference") or not ref_result.get("content"):
            continue

        new_references = [ref for ref in references_used if ref not in seen_references]
        if not new_references:
            continue

        seen_references.update(new_references)
        marker_preview = ", ".join(markers[:12])
        cluster_blocks.append(
            f"## Cluster {cluster_id} reference match\n"
            f"- Marker preview: {marker_preview}\n"
            f"- Preliminary cell type: {ref_result.get('preliminary_cell_type', 'Unknown')}\n"
            f"- References used: {', '.join(references_used)}\n\n"
            f"{ref_result['content']}"
        )

    if not cluster_blocks:
        info["reason"] = "No relevant references selected"
        return "", info

    combined = "\n\n---\n\n".join(cluster_blocks)
    if reference_max_context_length and len(combined) > reference_max_context_length:
        combined = combined[:reference_max_context_length] + "\n\n[... subcluster reference context truncated ...]"

    info["reference_used"] = True
    info["references_used"] = sorted(seen_references)
    info["reason"] = f"Selected {len(seen_references)} reference document(s)"

    context = (
        "<expert_reference>\n"
        "Expert-curated subtype references for this subclustering run. Use these "
        "as additional evidence for subtype differentiation, but prioritize the "
        "provided markers and tissue/species context when they conflict. Prefer "
        "consensus marker-program labels as the primary subtype; keep paper-specific "
        "labels only as traceability evidence when the marker support is strong.\n\n"
        f"Parent cluster context: {major_cluster_info}\n\n"
        f"{combined}\n"
        "</expert_reference>"
    )

    if verbose:
        print(f"Reference context added for subclustering: {', '.join(info['references_used'])}")

    return context, info


_REASONING_MODEL_HINTS = (
    "reasoner", "thinking", "think",
    "deepseek-v4", "deepseek-r1",
    "o1", "o3", "o4",
    "gpt-5", "gpt5",
    "claude-opus-4-5", "claude-opus-4-6", "claude-opus-4-7",
    "claude-sonnet-4-5", "claude-sonnet-4-6", "claude-sonnet-4-7",
    "gemini-3", "gemini-2.5-pro",
)


_THINK_TAG_PATTERN = re.compile(
    r'<\s*(think|thinking|reasoning)\b[^>]*>.*?</\s*\1\s*>',
    re.DOTALL | re.IGNORECASE,
)


def _strip_reasoning_blocks(text):
    """Remove <think>/<thinking>/<reasoning> wrapper blocks from LLM output.

    Some reasoning models (including newer DeepSeek and OpenRouter providers
    that forward reasoning to content) inline their hidden chain-of-thought
    inside the response. Stripping it before XML parsing prevents the
    reasoning text from being mistakenly matched or pushing the real answer
    out of the regex window.
    """
    if not text:
        return text
    return _THINK_TAG_PATTERN.sub("", str(text))


_CLUSTER_TAG_PROBE = re.compile(r'<\s*cluster\b', re.IGNORECASE)


def _xml_extract_with_retry(prompt, analysis_text, provider, model, temperature):
    """Call the LLM extractor; retry once with a stricter prompt if no <cluster> tag appears.

    Returns the best (most cluster-tag-rich) response so the downstream parser
    can still inspect raw output if both attempts fail.
    """
    first = subcluster_agent_annotate_subcluster(
        prompt, provider=provider, model=model, temperature=temperature
    )
    first_clean = _strip_reasoning_blocks(first or "")
    if _CLUSTER_TAG_PROBE.search(first_clean):
        return first

    retry_prompt = (
        "Your previous response did not contain the required XML structure.\n"
        "Output ONLY the XML blocks below. No preamble, no markdown fences, no commentary, "
        "no <think> tags. Wrap every cluster strictly as:\n\n"
        "<cluster id=\"CLUSTER_ID\">\n"
        "<celltype1>first cell type</celltype1>\n"
        "<celltype2>second cell type</celltype2>\n"
        "<reason>concise reason grounded in the markers</reason>\n"
        "</cluster>\n\n"
        "Use the exact Cluster ID from the analysis (quote it if it contains spaces).\n"
        "Include every cluster mentioned. Begin your reply with the first <cluster> tag.\n\n"
        "Analysis to convert:\n"
        f"{analysis_text}\n"
    )
    second = subcluster_agent_annotate_subcluster(
        retry_prompt, provider=provider, model=model, temperature=temperature
    )
    second_clean = _strip_reasoning_blocks(second or "")
    if _CLUSTER_TAG_PROBE.search(second_clean):
        print("  Stage-2 retry succeeded: recovered <cluster> XML on second attempt.")
        return second
    # Both attempts failed; return whichever has more raw content so downstream
    # debug output is most informative.
    return second if len(second or "") > len(first or "") else first


def _is_reasoning_model(model):
    """Best-effort detection for thinking/reasoning models.

    Reasoning models consume the token budget for hidden chain-of-thought, so
    a 4096 cap routinely truncates the final answer (empty content or
    half-written XML). Bumping to 12288 for these models avoids that without
    inflating cost for plain chat models.
    """
    if not model:
        return False
    name = str(model).lower()
    return any(hint in name for hint in _REASONING_MODEL_HINTS)


def subcluster_agent_annotate_subcluster(user_message, model=None, temperature=None, provider="openrouter"):
    """
    Unified function to call LLM for subcluster annotation.

    Args:
        user_message: The prompt message for subcluster annotation
        model: Model to use (defaults to provider's default if None)
        temperature: Temperature for generation (0-1)
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)

    Returns:
        The generated annotation as a string
    """
    # Set default model and temperature based on provider if not specified
    if model is None or temperature is None:
        defaults = get_agent_default("subclustering", provider)
        if model is None:
            model = defaults["model"]
        if temperature is None:
            temperature = defaults["temperature"]

    max_tokens = 12288 if _is_reasoning_model(model) else 4096

    # Use the unified call_llm function
    result = call_llm(
        prompt=user_message,
        provider=provider,
        model=model,
        temperature=temperature,
        max_tokens=max_tokens
    )

    return result if result else ''



def construct_prompt_from_csv_subcluster(marker, major_cluster_info, n_genes=50, additional_context=None):
    # Process DataFrame if it has more than 2 columns
    if len(marker.columns) > 2:
        print(f"Processing input dataframe to get top {n_genes} markers")
        get_top_markers = _get_get_top_markers()
        marker = get_top_markers(marker, n_genes=n_genes)
    else:
        print("Using input dataframe directly as it appears to be pre-processed (2 columns)")
        marker = marker.copy()

    # Initialize the prompt with the major cluster information
    prompt = f"""

You are an expert biologist specializing in cell type annotation, with deep expertise in immunology, cancer biology, and developmental biology. You will be given sets of highly expressed markers ranked by significance for some subclusters from the {major_cluster_info} cluster, identify what is the most likely top2 cell type each marker set implies.

Work step by step and ground every subtype call in the provided marker genes and parent-cluster context.
If additional context provides consensus macrophage/TAM programs, use those consensus program names as primary subtype labels when supported. If it provides paper-specific aliases, mention them only as evidence in the explanation.

For each output, provide:
1. Key marker:
2. Explanation:
3. Most likely top2 cell types:

Remember these subclusters are from a {major_cluster_info} big cluster. You must include all clusters mentioned in the analysis.
Return exactly one result for every Cluster ID listed below. Do not omit a
cluster, even if the markers look ambiguous, contaminating, or technically
stressed; instead, include that Cluster ID and explain the uncertainty.

The clusters are identified by their Cluster ID below:
"""

    # Iterate over each row in the DataFrame using actual cluster IDs
    for index, row in marker.iterrows():
        cluster_id = row.iloc[0]  # Use iloc for positional indexing
        markers = row.iloc[1]     # Use iloc for positional indexing
        prompt += f"Cluster {cluster_id}: {markers}\n"

    if additional_context:
        prompt += f"\n\nAdditional context that may help with the analysis:\n{additional_context}"

    return prompt



def annotate_subclusters(
    marker,
    major_cluster_info,
    model=None,
    temperature=None,
    provider="openrouter",
    n_genes=50,
    additional_context=None,
    tissue=None,
    species=None,
    use_reference=False,
    reference_provider=None,
    reference_model=None,
    reference_cell_type_hint=None,
    reference_depth="detailed",
    reference_max_content_length=5000,
    reference_max_context_length=12000,
):
    """
    Annotate subclusters using an LLM.

    Args:
        marker: DataFrame containing marker data
        major_cluster_info: Description of the major cluster type
        model: Model to use (defaults to provider's subclustering default)
        temperature: Temperature for generation (0-1)
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
        n_genes: Number of top genes to use
        additional_context: Optional context appended to the prompt
        tissue: Tissue type being analyzed. Optional.
        species: Species being analyzed. Optional.
        use_reference: Whether to retrieve expert subtype references for the
            subcluster marker sets before annotation.
        reference_provider: Provider for reference selection (default: provider)
        reference_model: Model for reference selection
        reference_cell_type_hint: Optional parent-lineage hint, e.g. "macrophage"

    Returns:
        The generated annotation as a string
    """
    context_parts = []
    if tissue:
        context_parts.append(f"Tissue: {tissue}")
    if species:
        context_parts.append(f"Species: {species}")
    if context_parts:
        tissue_species_context = ". ".join(context_parts) + "."
        additional_context = _combine_context(tissue_species_context, additional_context)

    if use_reference:
        reference_context, _ = build_subcluster_reference_context(
            marker=marker,
            major_cluster_info=major_cluster_info,
            provider=provider,
            n_genes=n_genes,
            tissue=tissue,
            species=species,
            reference_provider=reference_provider,
            reference_model=reference_model,
            reference_cell_type_hint=reference_cell_type_hint,
            reference_depth=reference_depth,
            reference_max_content_length=reference_max_content_length,
            reference_max_context_length=reference_max_context_length,
        )
        additional_context = _combine_context(additional_context, reference_context)

    prompt = construct_prompt_from_csv_subcluster(marker, major_cluster_info, n_genes=n_genes, additional_context=additional_context)
    output_text = subcluster_agent_annotate_subcluster(prompt, model=model, temperature=temperature, provider=provider)
    return output_text



def extract_subcluster_results_with_llm_multiple_output(analysis_text, provider="openrouter", model=None, temperature=None):
    """
    Extract multiple output results from subcluster analysis text using XML tags.
    Uses brief reasoning to save tokens in batch runs.

    Args:
        analysis_text: Text containing the analysis results
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
        model: Model to use
        temperature: Temperature for generation (0-1)

    Returns:
        Extracted results as string with XML-tagged format
    """
    # Define the prompt to instruct the LLM with XML output format (brief reasoning for batch)
    prompt = f"""You are an expert in analyzing celltype annotation for subclusters.

Extract the cell type annotations from the following analysis. For each cluster, output in this exact XML format:

<cluster id="CLUSTER_ID">
<celltype1>first cell type</celltype1>
<celltype2>second cell type</celltype2>
<reason>brief reason</reason>
</cluster>

IMPORTANT: Use the exact Cluster ID from the analysis (e.g., if the analysis mentions "Cluster 0", use id="0"; if it mentions "Cluster ABC", use id="ABC"). Do not renumber the clusters.

You should include all clusters mentioned in the analysis.

{analysis_text}
"""

    return _xml_extract_with_retry(prompt, analysis_text, provider, model, temperature)




def extract_subcluster_results_with_llm(analysis_text, provider="openrouter", model=None, temperature=None):
    """
    Extract results with reasons from subcluster analysis text using XML tags.

    Args:
        analysis_text: Text containing the analysis results
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
        model: Model to use
        temperature: Temperature for generation (0-1)

    Returns:
        Extracted results as string with XML-tagged format
    """
    # Define the prompt to instruct the LLM with XML output format
    prompt = f"""You are an expert in analyzing celltype annotation for subclusters.

Extract the cell type annotations from the following analysis. For each cluster, output in this exact XML format:

<cluster id="CLUSTER_ID">
<celltype1>first cell type</celltype1>
<celltype2>second cell type</celltype2>
<reason>the complete explanation from the analysis</reason>
</cluster>

IMPORTANT: Use the exact Cluster ID from the analysis (e.g., if the analysis mentions "Cluster 0", use id="0"; if it mentions "Cluster ABC", use id="ABC"). Do not renumber the clusters.

You should include all clusters mentioned in the analysis.

{analysis_text}
"""

    return _xml_extract_with_retry(prompt, analysis_text, provider, model, temperature)



def write_results_to_csv(results, output_name='subcluster_results', marker_map=None, expected_cluster_ids=None):
    """
    Extract cell type results from LLM output (XML format) and write to CSV file

    Args:
        results (str): LLM analysis results with XML-tagged clusters
        output_name (str): Base name for output file (will add .csv if not present)
        marker_map (dict): Optional mapping of cluster_id (str) -> marker genes string
        expected_cluster_ids (list): Optional list of expected cluster IDs from the marker data.
            When provided and LLM returns different IDs, remaps by position.

    Returns:
        pandas.DataFrame: DataFrame containing the extracted results
    """
    # Add .csv suffix if not present
    if not output_name.lower().endswith('.csv'):
        output_name = output_name + '.csv'

    results_str = _strip_reasoning_blocks(str(results))

    # Parse XML-tagged clusters: <cluster id="1">...</cluster>
    # Quoted ids may contain spaces ("cd8-positive, alpha-beta t cell"), so capture
    # quoted and unquoted variants separately and merge afterwards.
    cluster_pattern = r'<cluster\b[^>]*\bid=(?:"([^"]*)"|\'([^\']*)\'|([^\s>]+))[^>]*>(.*?)</cluster>'
    cluster_matches = re.findall(cluster_pattern, results_str, re.DOTALL | re.IGNORECASE)

    rows = []
    if cluster_matches:
        for q1, q2, unq, content in cluster_matches:
            cluster_id = (q1 or q2 or unq).strip()
            # Extract celltype1
            ct1_match = re.search(r'<celltype1>(.*?)</celltype1>', content, re.DOTALL | re.IGNORECASE)
            celltype1 = ct1_match.group(1).strip() if ct1_match else 'Unknown'

            # Extract celltype2
            ct2_match = re.search(r'<celltype2>(.*?)</celltype2>', content, re.DOTALL | re.IGNORECASE)
            celltype2 = ct2_match.group(1).strip() if ct2_match else 'Unknown'

            # Extract reason
            reason_match = re.search(r'<reason>(.*?)</reason>', content, re.DOTALL | re.IGNORECASE)
            reason = reason_match.group(1).strip() if reason_match else ''

            markers = marker_map.get(str(cluster_id), '') if marker_map else ''
            rows.append([cluster_id, celltype1, celltype2, markers, reason])

    if rows:
        df = pd.DataFrame(rows, columns=['Result ID', 'main_cell_type', 'sub_cell_type', 'key_markers', 'reason'])

        # Remap cluster IDs if LLM returned different ones than expected
        if expected_cluster_ids is not None:
            df['Result ID'] = df['Result ID'].astype(str)
            expected_list = [str(x) for x in expected_cluster_ids]
            expected_set = set(str(x) for x in expected_cluster_ids)
            result_set = set(df['Result ID'].astype(str))
            if len(df) == len(expected_cluster_ids) and expected_set != result_set:
                print(f"  Remapping cluster IDs: {df['Result ID'].tolist()} -> {expected_list}")
                df['Result ID'] = expected_list
                # Also re-populate key_markers with correct mapping
                if marker_map:
                    df['key_markers'] = [marker_map.get(cid, '') for cid in expected_list]
                result_set = set(df['Result ID'].astype(str))

            missing_clusters = [cid for cid in expected_list if cid not in result_set]
            if missing_clusters:
                print(f"  Warning: LLM output omitted cluster IDs: {missing_clusters}")
                missing_rows = pd.DataFrame([
                    {
                        'Result ID': cid,
                        'main_cell_type': 'Unknown',
                        'sub_cell_type': 'Missing result',
                        'key_markers': marker_map.get(cid, '') if marker_map else '',
                        'reason': 'MISSING_RESULT: LLM did not return an annotation for this cluster.',
                    }
                    for cid in missing_clusters
                ])
                df = pd.concat([df, missing_rows], ignore_index=True)

            order = {cid: i for i, cid in enumerate(expected_list)}
            if order:
                df['_expected_order'] = df['Result ID'].map(lambda cid: order.get(str(cid), len(order)))
                df = df.sort_values('_expected_order').drop(columns=['_expected_order']).reset_index(drop=True)

        df.to_csv(output_name, index=False)
        print(f"Results have been written to {output_name}")
        return df

    # If XML parsing failed, save raw results for debugging
    raw_output_file = f"{output_name}.txt"
    with open(raw_output_file, "w") as f:
        f.write(results_str)

    raise RuntimeError(
        f"\n{'='*60}\n"
        f"SUBCLUSTERING FAILED - Could not parse LLM results\n"
        f"{'='*60}\n"
        f"Raw output saved to: {raw_output_file}\n"
        f"Expected XML format: <cluster id=\"1\"><celltype1>...</celltype1>...</cluster>\n"
        f"{'='*60}"
    )



def runCASSIA_subclusters(marker, major_cluster_info, output_name,
                       model=None, temperature=None, provider="openrouter", n_genes=50,
                       tissue=None, species=None, additional_context=None,
                       use_reference=False, reference_provider=None, reference_model=None,
                       reference_cell_type_hint=None, reference_depth="detailed",
                       reference_max_content_length=5000,
                       reference_max_context_length=12000):
    """
    Process subclusters from marker data and generate annotated results.

    Args:
        marker: DataFrame containing marker data
        major_cluster_info: Description of the major cluster type
        output_name: Base name for output file (will add .csv if not present)
        model: Model name to use (defaults to provider's subclustering default)
        temperature: Temperature parameter for API calls (0-1)
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
        n_genes: Number of top genes to use for analysis
        tissue: Tissue type being analyzed (e.g., "lung", "brain"). Optional.
        species: Species being analyzed (e.g., "human", "mouse"). Optional.
        additional_context: Optional context appended to the prompt.
        use_reference: Whether to retrieve expert subtype references for
            subcluster marker sets before annotation.
        reference_provider: Provider for reference selection (default: provider).
        reference_model: Model for reference selection.
        reference_cell_type_hint: Optional parent-lineage hint, e.g. "macrophage".

    Returns:
        None: Results are saved to a CSV file
    """
    # Build tissue/species context and prepend to additional_context if provided
    context_parts = []
    if tissue:
        context_parts.append(f"Tissue: {tissue}")
    if species:
        context_parts.append(f"Species: {species}")
    if context_parts:
        tissue_species_context = ". ".join(context_parts) + "."
        if additional_context:
            additional_context = tissue_species_context + " " + additional_context
        else:
            additional_context = tissue_species_context

    if use_reference:
        reference_context, reference_info = build_subcluster_reference_context(
            marker=marker,
            major_cluster_info=major_cluster_info,
            provider=provider,
            n_genes=n_genes,
            tissue=tissue,
            species=species,
            reference_provider=reference_provider,
            reference_model=reference_model,
            reference_cell_type_hint=reference_cell_type_hint,
            reference_depth=reference_depth,
            reference_max_content_length=reference_max_content_length,
            reference_max_context_length=reference_max_context_length,
            verbose=True,
        )
        additional_context = _combine_context(additional_context, reference_context)
        if reference_info.get("reference_used"):
            print(f"Using reference documents: {', '.join(reference_info['references_used'])}")
        else:
            print(f"Reference retrieval did not add context: {reference_info.get('reason', 'No match')}")

    # Apply agent defaults if model or temperature not specified
    if model is None or temperature is None:
        defaults = get_agent_default("subclustering", provider)
        if model is None:
            model = defaults["model"]
        if temperature is None:
            temperature = defaults["temperature"]

    # Construct prompt and get analysis from LLM
    prompt = construct_prompt_from_csv_subcluster(marker, major_cluster_info, n_genes=n_genes, additional_context=additional_context)
    output_text = subcluster_agent_annotate_subcluster(prompt, model=model, temperature=temperature, provider=provider)

    # Extract structured results from the analysis text
    results = extract_subcluster_results_with_llm(output_text, provider=provider, model=model, temperature=temperature)

    # Build marker_map from input data for populating key_markers column
    try:
        get_top_markers = _get_get_top_markers()
        marker_df = get_top_markers(marker, n_genes=n_genes) if len(marker.columns) > 2 else marker.copy()
    except Exception:
        marker_df = marker.copy()
    marker_map = {str(row.iloc[0]): str(row.iloc[1]) for _, row in marker_df.iterrows()}
    expected_cluster_ids = [str(row.iloc[0]) for _, row in marker_df.iterrows()]

    # Save results to CSV
    write_results_to_csv(results, output_name, marker_map=marker_map, expected_cluster_ids=expected_cluster_ids)

    # --- Generate HTML report for the single run CSV ---
    try:
        from CASSIA.reports.generate_reports import process_evaluation_csv
    except ImportError:
        try:
            from ...reports.generate_reports import process_evaluation_csv
        except ImportError:
            from reports.generate_reports import process_evaluation_csv
    import os
    csv_file = output_name if output_name.lower().endswith('.csv') else output_name + '.csv'
    if os.path.exists(csv_file):
        process_evaluation_csv(csv_file, overwrite=True, model_name=model)
    
    return None



def runCASSIA_n_subcluster(n, marker, major_cluster_info, base_output_name,
                          model=None, temperature=None,
                          provider="openrouter", max_workers=5, n_genes=50,
                          tissue=None, species=None, additional_context=None,
                          use_reference=False, reference_provider=None, reference_model=None,
                          reference_cell_type_hint=None, reference_depth="detailed",
                          reference_max_content_length=5000,
                          reference_max_context_length=12000):
    """
    Run multiple subcluster analyses in parallel and save results.

    Args:
        n: Number of analyses to run
        marker: DataFrame containing marker data
        major_cluster_info: Description of the major cluster type
        base_output_name: Base name for output files
        model: Model name to use (defaults to provider's subclustering_n default)
        temperature: Temperature parameter for API calls (defaults to subclustering_n default: 0.3)
        provider: LLM provider ("openai", "anthropic", "openrouter", or a custom API URL)
        max_workers: Maximum number of parallel workers
        n_genes: Number of top genes to use for analysis
        tissue: Tissue type being analyzed (e.g., "lung", "brain"). Optional.
        species: Species being analyzed (e.g., "human", "mouse"). Optional.
        additional_context: Optional context appended to the prompt.
        use_reference: Whether to retrieve expert subtype references for the
            subcluster marker sets before repeated annotation.
        reference_provider: Provider for reference selection (default: provider).
        reference_model: Model for reference selection.
        reference_cell_type_hint: Optional parent-lineage hint, e.g. "macrophage".

    Returns:
        None: Results are saved to CSV files
    """
    # Build tissue/species context and prepend to additional_context if provided
    context_parts = []
    if tissue:
        context_parts.append(f"Tissue: {tissue}")
    if species:
        context_parts.append(f"Species: {species}")
    if context_parts:
        tissue_species_context = ". ".join(context_parts) + "."
        if additional_context:
            additional_context = tissue_species_context + " " + additional_context
        else:
            additional_context = tissue_species_context

    if use_reference:
        reference_context, reference_info = build_subcluster_reference_context(
            marker=marker,
            major_cluster_info=major_cluster_info,
            provider=provider,
            n_genes=n_genes,
            tissue=tissue,
            species=species,
            reference_provider=reference_provider,
            reference_model=reference_model,
            reference_cell_type_hint=reference_cell_type_hint,
            reference_depth=reference_depth,
            reference_max_content_length=reference_max_content_length,
            reference_max_context_length=reference_max_context_length,
            verbose=True,
        )
        additional_context = _combine_context(additional_context, reference_context)
        if reference_info.get("reference_used"):
            print(f"Using reference documents: {', '.join(reference_info['references_used'])}")
        else:
            print(f"Reference retrieval did not add context: {reference_info.get('reason', 'No match')}")

    # Apply agent defaults for n-times subclustering (uses subclustering_n for variability)
    if model is None or temperature is None:
        defaults = get_agent_default("subclustering_n", provider)
        if model is None:
            model = defaults["model"]
        if temperature is None:
            temperature = defaults["temperature"]
    
    # Build marker_map from input data for populating key_markers column
    try:
        get_top_markers = _get_get_top_markers()
        marker_df_for_map = get_top_markers(marker, n_genes=n_genes) if len(marker.columns) > 2 else marker.copy()
    except Exception:
        marker_df_for_map = marker.copy()
    marker_map = {str(row.iloc[0]): str(row.iloc[1]) for _, row in marker_df_for_map.iterrows()}

    def run_single_analysis(i):
        # Run the annotation process
        output_text = annotate_subclusters(marker, major_cluster_info,
                                         model=model, temperature=temperature, provider=provider, n_genes=n_genes,
                                         additional_context=additional_context,
                                         tissue=None, species=None, use_reference=False)

        # Extract results using XML format
        results = extract_subcluster_results_with_llm_multiple_output(output_text, provider=provider, model=model, temperature=temperature)
        results_str = _strip_reasoning_blocks(str(results))

        # Parse XML-tagged clusters: <cluster id="1">...</cluster>
        # Quoted ids may contain spaces, so capture quoted and unquoted variants separately.
        cluster_pattern = r'<cluster\b[^>]*\bid=(?:"([^"]*)"|\'([^\']*)\'|([^\s>]+))[^>]*>(.*?)</cluster>'
        cluster_matches = re.findall(cluster_pattern, results_str, re.DOTALL | re.IGNORECASE)

        rows = []
        if cluster_matches:
            for q1, q2, unq, content in cluster_matches:
                cluster_id = (q1 or q2 or unq).strip()
                # Extract celltype1
                ct1_match = re.search(r'<celltype1>(.*?)</celltype1>', content, re.DOTALL | re.IGNORECASE)
                celltype1 = ct1_match.group(1).strip() if ct1_match else 'Unknown'

                # Extract celltype2
                ct2_match = re.search(r'<celltype2>(.*?)</celltype2>', content, re.DOTALL | re.IGNORECASE)
                celltype2 = ct2_match.group(1).strip() if ct2_match else 'Unknown'

                # Extract reason
                reason_match = re.search(r'<reason>(.*?)</reason>', content, re.DOTALL | re.IGNORECASE)
                reason = reason_match.group(1).strip() if reason_match else ''

                markers = marker_map.get(str(cluster_id), '')
                rows.append([cluster_id, celltype1, celltype2, markers, reason])

        if rows:
            df = pd.DataFrame(rows, columns=['Result ID', 'main_cell_type', 'sub_cell_type', 'key_markers', 'reason'])
        else:
            raise RuntimeError(
                f"Iteration {i+1}: Could not parse LLM response. "
                f"Expected XML format: <cluster id=\"1\"><celltype1>...</celltype1>...</cluster>. "
                f"Response preview: {results_str[:200]}"
            )

        try:
            # Get the marker dataframe with cluster IDs for validation
            try:
                get_top_markers = _get_get_top_markers()
                marker_df = get_top_markers(marker, n_genes=n_genes)
            except KeyError as e:
                # If get_top_markers fails due to missing columns, use the original marker dataframe
                print(f"Warning: {str(e)}. Using original marker dataframe.")
                marker_df = marker.copy()

            # Convert types to ensure compatibility
            df['Result ID'] = df['Result ID'].astype(str)

            # Extract expected cluster IDs from marker_df
            expected_cluster_ids_list = [str(x) for x in marker_df.iloc[:, 0].tolist()]
            expected_cluster_ids = set(expected_cluster_ids_list)
            result_cluster_ids = set(str(x) for x in df['Result ID'].tolist())

            # Check for mismatched cluster IDs
            missing_clusters = expected_cluster_ids - result_cluster_ids
            unexpected_clusters = result_cluster_ids - expected_cluster_ids

            if missing_clusters:
                print(f"Warning: Iteration {i+1} - LLM output missing cluster IDs: {missing_clusters}")
            if unexpected_clusters:
                print(f"Warning: Iteration {i+1} - LLM output has unexpected cluster IDs: {unexpected_clusters}")

            # If LLM returned different IDs (e.g., numeric "1","2" instead of
            # "monocyte","plasma cell"), remap by position to the original IDs
            if missing_clusters and unexpected_clusters and len(df) == len(expected_cluster_ids_list):
                print(f"  Remapping cluster IDs by position: {df['Result ID'].tolist()} -> {expected_cluster_ids_list}")
                df['Result ID'] = expected_cluster_ids_list[:len(df)]
                # Re-populate key_markers using the correct (remapped) cluster IDs
                df['key_markers'] = [marker_map.get(cid, '') for cid in expected_cluster_ids_list[:len(df)]]

        except Exception as e:
            print(f"Warning: Error during cluster ID validation: {str(e)}. Continuing anyway.")

        # Write the DataFrame to a CSV file with an index
        indexed_csv_file_path = f'{base_output_name}_{i+1}.csv'
        df.to_csv(indexed_csv_file_path, index=False)
        
        return indexed_csv_file_path

    failed_iterations = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(run_single_analysis, i): i for i in range(n)}
        result_files = []
        for future in as_completed(futures):
            i = futures[future]
            try:
                result_file = future.result()
                print(f"Results for iteration {i+1} have been written to {result_file}")
                result_files.append(result_file)
            except Exception as exc:
                print(f"Iteration {i+1} generated an exception: {exc}")
                failed_iterations.append((i+1, str(exc)))

    # Check if ALL iterations failed - this is a total failure
    if len(result_files) == 0 and len(failed_iterations) > 0:
        error_sample = failed_iterations[0][1][:200] if failed_iterations else "Unknown"
        raise RuntimeError(
            f"\n{'='*60}\n"
            f"SUBCLUSTERING BATCH FAILED - All {n} iterations failed\n"
            f"{'='*60}\n"
            f"Sample error: {error_sample}\n"
            f"{'='*60}"
        )

    # Warn about partial failures
    if failed_iterations:
        print(f"\nWarning: {len(failed_iterations)} of {n} iterations failed:")
        for iter_num, err in failed_iterations[:5]:
            print(f"  - Iteration {iter_num}: {err[:80]}...")
        if len(failed_iterations) > 5:
            print(f"  ... and {len(failed_iterations) - 5} more")

    # --- Generate HTML reports for all batch CSVs ---
    try:
        from CASSIA.reports.generate_reports import process_evaluation_csv, create_index_html
    except ImportError:
        try:
            from ...reports.generate_reports import process_evaluation_csv, create_index_html
        except ImportError:
            from reports.generate_reports import process_evaluation_csv, create_index_html
    import os

    # Generate HTML report for each CSV
    for csv_file in result_files:
        if os.path.exists(csv_file):
            process_evaluation_csv(csv_file, overwrite=True)

    # Create an index.html summary in the same directory as the first result file
    if result_files:
        output_dir = os.path.dirname(result_files[0]) or '.'
        create_index_html(result_files, output_dir)
        print(f"Batch HTML reports and index generated in {output_dir}")


def test_custom_api_parsing():
    """
    Test function to simulate a response from a custom API provider and test the parsing functionality.
    This is useful for debugging the parsing logic without making actual API calls.
    """
    import pandas as pd
    import os
    try:
        from CASSIA.engine.tools_function import get_top_markers
    except ImportError:
        try:
            from .tools_function import get_top_markers
        except ImportError:
            from tools_function import get_top_markers

    # Sample structured response that mimics what DeepSeek or other custom APIs might return
    sample_response = [
        {
            'cluster': 1,
            'key_markers': 'IL7R, CD8A, CD8B, CCL4, KLRB1, ITK',
            'explanation': 'The presence of IL7R, CD8A, and CD8B suggests a CD8+ T cell identity. CCL4 is associated with effector functions, while KLRB1 (CD161) and ITK indicate a memory-like or tissue-resident phenotype.',
            'most_likely_top2_cell_types': ['CD8+ memory T cells', 'Tissue-resident memory CD8+ T cells (TRM)']
        },
        {
            'cluster': 2,
            'key_markers': 'LAYN, HAVCR2 (TIM-3), TIGIT, IKZF2, KLRC2, KLRC3',
            'explanation': 'LAYN (Lag-3) and HAVCR2 (TIM-3) are markers of exhausted or chronically stimulated CD8+ T cells.',
            'most_likely_top2_cell_types': ['Exhausted CD8+ T cells', 'NK-like CD8+ T cells']
        },
        {
            'cluster': 3,
            'key_markers': 'GZMK, GZMH, PRF1, NKG7, CCR7, CD27',
            'explanation': 'GZMK, GZMH, PRF1, and NKG7 are markers of cytotoxic activity, typical of effector CD8+ T cells.',
            'most_likely_top2_cell_types': ['Effector CD8+ T cells', 'Central memory CD8+ T cells']
        },
        {
            'cluster': 4,
            'key_markers': 'WFDC2, CEACAM7, CLDN8, PPARG, HOXD13, HOXB13',
            'explanation': 'WFDC2, CEACAM7, and CLDN8 are markers associated with epithelial or secretory cells.',
            'most_likely_top2_cell_types': ['Regulatory CD8+ T cells', 'Epithelial-like CD8+ T cells (rare subset)']
        }
    ]
    
    # Test the write_results_to_csv function
    print("Testing write_results_to_csv with structured data...")
    df = write_results_to_csv(sample_response, output_name='test_parsing_result')
    print(f"Generated DataFrame:\n{df}")
    
    # Test the extract_subcluster_results_with_llm function
    print("\nTesting extract_subcluster_results_with_llm with structured data...")
    result = extract_subcluster_results_with_llm(sample_response)
    print(f"Result type: {type(result)}")
    
    # Test the extract_subcluster_results_with_llm_multiple_output function
    print("\nTesting extract_subcluster_results_with_llm_multiple_output with structured data...")
    result = extract_subcluster_results_with_llm_multiple_output(sample_response)
    print(f"Result type: {type(result)}")
    
    # Create a simple marker dataframe for testing
    print("\nCreating test marker dataframe...")
    test_marker_df = pd.DataFrame({
        'cluster': ['cluster1', 'cluster2', 'cluster3', 'cluster4', 'cluster5', 'cluster6'],
        'markers': [
            'IL7R, CD8A, CD8B, CCL4, KLRB1, ITK',
            'LAYN, HAVCR2, TIGIT, IKZF2, KLRC2, KLRC3',
            'GZMK, GZMH, PRF1, NKG7, CCR7, CD27',
            'WFDC2, CEACAM7, CLDN8, PPARG, HOXD13, HOXB13',
            'GNLY, KLRF1, FCER1G, TYROBP, CD38, KIR2DL4',
            'LPL, SNAI2, HAND2, SOX2, NES, PDGFRA'
        ]
    })
    print(f"Created test marker dataframe with shape: {test_marker_df.shape}")
    
    # Save the test dataframe to a temporary file
    temp_file = 'test_markers_temp.csv'
    test_marker_df.to_csv(temp_file, index=False)
    print(f"Saved test marker dataframe to {temp_file}")
    
    try:
        # Create a mock function to simulate annotate_subclusters
        def mock_annotate(*args, **kwargs):
            return sample_response
            
        # Save the original function
        original_annotate = globals()['annotate_subclusters']
        
        # Replace with mock function
        globals()['annotate_subclusters'] = mock_annotate
        
        print("\nTesting runCASSIA_n_subcluster with mock data...")
        # Run with n=1 to test a single iteration
        runCASSIA_n_subcluster(
            n=1,
            marker=test_marker_df,
            major_cluster_info="cd8 t cell",
            base_output_name="test_n_subcluster",
            model="dummy-model",
            provider="dummy-provider"
        )
        
        # Check if the output file was created
        output_file = "test_n_subcluster_1.csv"
        if os.path.exists(output_file):
            print(f"Successfully created output file: {output_file}")
            result_df = pd.read_csv(output_file)
            print(f"Output file contents:\n{result_df}")
        else:
            print(f"Error: Output file {output_file} was not created")
        
        # Restore the original function
        globals()['annotate_subclusters'] = original_annotate
        
        print("\nAll tests completed.")
        
    except Exception as e:
        print(f"Error during testing: {str(e)}")
        import traceback
        traceback.print_exc()
    finally:
        # Clean up temporary files
        if os.path.exists(temp_file):
            os.remove(temp_file)
            print(f"Removed temporary file: {temp_file}")

# Uncomment to run the test function when this file is executed directly
# if __name__ == "__main__":
#     test_custom_api_parsing()
