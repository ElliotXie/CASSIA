"""
Reference Agent - Main Orchestrator for CASSIA.

Intelligent reference-document retrieval and context injection for cell-type
annotation. For subclustering, the agent reads the lineage overview/router,
plans which detailed documents to read for all subclusters together, then
synthesizes a compact literature-grounded brief for prompt injection.
Single-cluster annotation still uses the compatible one-call selector.

The caller opts into reference mode explicitly (use_reference=True in the
main pipeline), so there is no separate "do we need a reference?" gating
step — if the library has nothing relevant the selector returns an empty
list and the agent reports should_use_reference=False.
"""

import json
import re
from pathlib import Path
from typing import Dict, List, Optional
try:
    from CASSIA.core.llm_utils import call_llm
except ImportError:
    try:
        from ...core.llm_utils import call_llm
    except ImportError:
        from llm_utils import call_llm

try:
    from .complexity_scorer import select_references_llm
    from .section_extractor import extract_sections
    from .utils import (
        get_references_dir,
        load_reference_index,
        load_markdown_file,
        format_reference_content,
    )
except ImportError:
    from complexity_scorer import select_references_llm
    from section_extractor import extract_sections
    from utils import (
        get_references_dir,
        load_reference_index,
        load_markdown_file,
        format_reference_content,
    )


HUMAN_CANCER_MACROPHAGE_CONSENSUS_STATES = [
    {
        "id": "spp1_areg_inflammatory_angiogenic_tam",
        "consensus_name": "SPP1/AREG inflammatory angiogenic TAM",
        "markers": ["SPP1", "AREG", "EREG", "CXCL3", "CXCL2", "CXCL8", "CXCL1", "CCL20", "IL1B", "TIMP1"],
        "aliases": ["Cheng SPP1+ angiogenesis-associated TAM", "Coulton 6_SPP1AREGMac"],
        "confidence": "high",
    },
    {
        "id": "apoe_trem2_lipid_phagolysosomal_tam",
        "consensus_name": "APOE/TREM2 lipid-phagolysosomal TAM",
        "markers": ["TREM2", "APOE", "APOC1", "GPNMB", "LIPA", "CTSD", "LGMN", "PLA2G7", "ACP5", "PSAP"],
        "aliases": ["TREM2+ lipid-associated TAM", "Coulton 4_ICIMac2"],
        "confidence": "high",
    },
    {
        "id": "areg_thbs1_angiogenic_remodeling_tam",
        "consensus_name": "AREG/THBS1 angiogenic remodeling TAM",
        "markers": ["AREG", "THBS1", "EREG", "NAMPT", "ZNF331", "IL1B", "GPR183", "NR4A3", "G0S2", "BTG1"],
        "aliases": ["EGF-ligand angiogenic TAM", "Coulton 9_AngioMac"],
        "confidence": "medium-high",
    },
    {
        "id": "folr2_selenop_resident_like_tam",
        "consensus_name": "FOLR2/SELENOP resident-like iron-handling TAM",
        "markers": ["FOLR2", "SELENOP", "SLC40A1", "F13A1", "STAB1", "RNASE1", "CD163", "MRC1", "LYVE1", "DAB2"],
        "aliases": ["resident-like TAM", "Coulton 1_MetM2Mac"],
        "confidence": "high",
    },
    {
        "id": "c1qc_c3_antigen_presenting_tam",
        "consensus_name": "C1QC/C3 complement antigen-presenting TAM",
        "markers": ["C1QA", "C1QB", "C1QC", "C3", "CD74", "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CX3CR1", "PLD4"],
        "aliases": ["Cheng C1QC+ phagocytic TAM", "Coulton 2_C3Mac"],
        "confidence": "high",
    },
    {
        "id": "cxcl9_cxcl10_ifng_response_tam",
        "consensus_name": "CXCL9/CXCL10 IFNG-response TAM",
        "markers": ["CXCL9", "CXCL10", "GBP1", "GBP5", "STAT1", "WARS1", "SLAMF7", "LGALS2", "VAMP5", "MMP9"],
        "aliases": ["IFNGMac", "Coulton 8_IFNGMac"],
        "confidence": "high",
    },
    {
        "id": "isg15_ifit_type_i_ifn_tam",
        "consensus_name": "ISG15/IFIT type-I interferon TAM",
        "markers": ["ISG15", "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM3", "MX1", "RSAD2", "CXCL10", "TNFSF10"],
        "aliases": ["Cheng ISG15+ TAM", "Coulton 17_IFNMac3"],
        "confidence": "high",
    },
    {
        "id": "ifi27_apoe_c1q_interferon_lipid_tam",
        "consensus_name": "IFI27/APOE/C1Q interferon-lipid TAM",
        "markers": ["IFI27", "APOE", "APOC1", "C1QA", "C1QB", "C1QC", "GPNMB", "TREM2", "A2M", "FTL"],
        "aliases": ["Li 2024 Macro_IFI27", "interferon-exposed lipid/complement TAM"],
        "confidence": "medium",
    },
    {
        "id": "il1b_tnf_inflammatory_tam",
        "consensus_name": "IL1B/TNF inflammatory TAM",
        "markers": ["IL1B", "TNF", "CXCL8", "CXCL1", "CXCL2", "CXCL3", "CCL3", "CCL4", "CCL3L3", "CCL4L2"],
        "aliases": ["acute inflammatory TAM", "Coulton 10_InflamMac"],
        "confidence": "high",
    },
    {
        "id": "hmox1_slc40a1_heme_iron_tam",
        "consensus_name": "HMOX1/SLC40A1 heme-iron handling TAM",
        "markers": ["HMOX1", "SLC40A1", "CD163", "CCL18", "LGMN", "CTSB", "CTSL", "HAMP", "STAB1", "SELENOP"],
        "aliases": ["erythrophagocytic macrophage", "Coulton 21_HemeMac"],
        "confidence": "medium-high",
    },
    {
        "id": "mt1_mt2_metallothionein_stress_tam",
        "consensus_name": "MT1/MT2 metallothionein stress TAM",
        "markers": ["MT1G", "MT1X", "MT2A", "MT1E", "MT1H", "MT1F", "MT1M", "MIF", "SPP1", "LDHA"],
        "aliases": ["metal-ion stress TAM", "Coulton 11_MetalloMac"],
        "confidence": "medium-high",
    },
    {
        "id": "hspa_heat_shock_stress_macrophage",
        "consensus_name": "HSPA heat-shock stress macrophage",
        "markers": ["HSPA6", "HSPA1B", "HSPA1A", "DNAJB1", "HSPB1", "BAG3", "HSPH1", "HSP90AA1", "ZFAND2A", "HSPD1"],
        "aliases": ["dissociation/stress-response macrophage", "Coulton 5_StressMac"],
        "confidence": "medium",
    },
    {
        "id": "mki67_proliferating_macrophage",
        "consensus_name": "MKI67/TOP2A proliferating macrophage",
        "markers": ["MKI67", "TOP2A", "STMN1", "PCLAF", "UBE2C", "TYMS", "CENPF", "HMGB2", "TUBB", "H2AZ1"],
        "aliases": ["cycling macrophage", "Coulton 14_ProliMac"],
        "confidence": "medium",
    },
    {
        "id": "fcn1_s100a8_monocyte_like_macrophage",
        "consensus_name": "FCN1/S100A8 monocyte-like inflammatory macrophage",
        "markers": ["FCN1", "S100A8", "S100A9", "S100A12", "VCAN", "LYZ", "LST1", "CCR2", "IL1B", "EREG"],
        "aliases": ["recruited monocyte-like macrophage", "Coulton 19_ClassMono"],
        "confidence": "medium",
    },
    {
        "id": "col1a1_sparc_ecm_macrophage_like_state",
        "consensus_name": "COL1A1/SPARC ECM-remodeling macrophage-like state",
        "markers": ["COL1A1", "COL1A2", "COL3A1", "SPARC", "LUM", "DCN", "POSTN", "BGN", "MGP", "SFRP2"],
        "aliases": ["ECMMac/METAM-like", "Coulton 18_ECMMac"],
        "confidence": "ambiguous",
    },
]


SUBCLUSTER_REFERENCE_PLAN_PROMPT = """You are the CASSIA Reference Agent. You do not annotate the final cell types yet. Your job is to decide which expert reference documents should be read before subtype annotation.

You have already been given all subcluster marker sets from one parent cluster. First use the lineage overview as your router, then choose only the detailed reference files that are needed. Do not select every file by default.

## Parent Cluster Context
{major_cluster_info}

## Tissue / Species
- Tissue: {tissue}
- Species: {species}
- Cell type hint: {cell_type_hint}

## All Subcluster Marker Sets
{cluster_marker_text}

## Lineage Overview / Router
{overview_content}

## Available Detailed Reference Files
{available_reference_text}

## Task
1. Form a short global hypothesis for the subtype landscape across these clusters.
2. Pick the smallest set of detailed reference files to read next.
3. For each cluster, write a preliminary hypothesis and mention which reference file(s) would help.

## Response Format
Return JSON only:
{{
  "global_hypothesis": "...",
  "selected_references": ["myeloid/macrophage/tam_pan_cancer.md"],
  "cluster_hypotheses": [
    {{
      "cluster_id": "...",
      "preliminary_pattern": "...",
      "supporting_markers": ["GENE1", "GENE2"],
      "reference_targets": ["myeloid/macrophage/tam_pan_cancer.md"]
    }}
  ],
  "reasoning": "..."
}}"""


SUBCLUSTER_REFERENCE_BRIEF_PROMPT = """You are the CASSIA Reference Agent writing a case-specific expert reference brief for a downstream cell-type annotation agent.

Do not dump the reference documents. Use them to synthesize a compact, objective, literature-grounded brief. The downstream agent will still see the original cluster marker sets, so your job is to explain what the literature says and how it maps to this specific subclustering problem.

## Parent Cluster Context
{major_cluster_info}

## Tissue / Species
- Tissue: {tissue}
- Species: {species}
- Cell type hint: {cell_type_hint}

## All Subcluster Marker Sets
{cluster_marker_text}

## Reference Planning Result
{planning_json}

## Reference Documents Read
{reference_documents}

## Reference Marker Match Summary
{reference_marker_matches}

## Required Output
Write a concise markdown block with these sections:

<reference_brief>
## Reference Agent Summary
One paragraph summarizing the likely subtype landscape.

## Objective Literature Facts
Bullets with factual statements from the cited papers. Include the paper name or DOI in each bullet when available. Prefer facts about marker programs, named subtype states, tissue/disease context, and known distinctions.

## Cluster-Specific Guidance
For each cluster, give:
- likely consensus lineage/state program as the primary label
- supporting markers
- important alternatives or conflicts to avoid
- paper-specific aliases only as evidence or traceability, not as universal ontology

## Cross-Cluster Distinctions
Short bullets explaining how to distinguish confusing states in this run.
</reference_brief>

Keep the brief specific to these marker sets. Avoid unsupported claims. If evidence is weak for a cluster, say so. Prefer consensus program labels such as "TOX/CXCL13 exhausted CD8 T cell", "FOXP3/CTLA4 regulatory T cell", or "CXCL9/CXCL10 IFNG-response TAM"; write paper-specific labels as parenthetical evidence, not as the primary universal subtype."""


class ReferenceAgent:
    """Agent for reference document retrieval and context injection."""

    def __init__(
        self,
        reference_dir: Optional[str] = None,
        provider: str = "openrouter",
        model: Optional[str] = None,
        api_key: Optional[str] = None,
    ):
        if reference_dir:
            self.reference_dir = Path(reference_dir)
        else:
            self.reference_dir = get_references_dir()

        self.provider = provider
        self.model = model
        self.api_key = api_key
        self._index = None

    @property
    def index(self) -> Dict:
        if self._index is None:
            self._index = load_reference_index()
        return self._index

    def get_reference_for_markers(
        self,
        markers: List[str],
        tissue: Optional[str] = None,
        species: Optional[str] = None,
        cell_type_hint: Optional[str] = None,
        depth: str = "detailed",
        max_content_length: Optional[int] = 8000,
    ) -> Dict:
        """Select references for a marker set and return formatted content.

        Returns a dict with:
            - should_use_reference: bool
            - content: str (markdown to inject, empty if no reference used)
            - references_used: List[str]
            - preliminary_cell_type: str
            - cell_type_range: List[str]
            - selected_references: List[str]
            - reasoning: str
        """
        if not markers:
            return self._empty_result("No markers provided")

        markers = markers[:20]

        selection = select_references_llm(
            markers=markers,
            tissue=tissue,
            species=species,
            cell_type_hint=cell_type_hint,
            provider=self.provider,
            model=self.model,
            api_key=self.api_key,
        )

        result = {
            "should_use_reference": False,
            "content": "",
            "references_used": [],
            "preliminary_cell_type": selection.get("preliminary_cell_type", "Unknown"),
            "cell_type_range": selection.get("cell_type_range", []),
            "selected_references": selection.get("selected_references", []),
            "reasoning": selection.get("reasoning", ""),
        }

        selected_paths = selection.get("selected_references", [])
        if not selected_paths:
            result["reasoning"] += " No relevant references in library."
            return result

        all_content = []
        for ref_path_str in selected_paths:
            ref_path = self.reference_dir / ref_path_str
            if not ref_path.exists():
                ref_path_str = ref_path_str.lstrip("/").lstrip("\\")
                ref_path = self.reference_dir / ref_path_str
            if not ref_path.exists():
                continue

            extraction = extract_sections(
                file_path=ref_path,
                cell_type_guess=result["preliminary_cell_type"],
                markers=markers,
                depth=depth,
            )

            if extraction["full_content"]:
                ref_name = (
                    ref_path_str.replace("/", "_").replace("\\", "_").replace(".md", "")
                )
                all_content.append(
                    f"### Reference: {ref_name}\n\n{extraction['full_content']}"
                )
                result["references_used"].append(ref_path_str)

        if all_content:
            combined = "\n\n---\n\n".join(all_content)
            result["content"] = format_reference_content(combined, max_content_length)
            result["should_use_reference"] = True
        else:
            result["reasoning"] += " No content extracted from selected references."

        return result

    def get_reference_brief_for_subclusters(
        self,
        marker_sets: List[Dict],
        major_cluster_info: str,
        tissue: Optional[str] = None,
        species: Optional[str] = None,
        cell_type_hint: Optional[str] = None,
        depth: str = "detailed",
        max_reference_content_length: int = 12000,
        max_brief_length: Optional[int] = 8000,
    ) -> Dict:
        """Generate an agentic, case-specific reference brief for subclusters.

        The agent first reads a lineage overview/router and all cluster marker
        sets, then selects the detailed reference files to read. A second LLM
        call synthesizes a compact reference brief with objective literature
        facts and cluster-specific guidance.

        Args:
            marker_sets: List of dicts with ``cluster_id`` and ``markers``.
            major_cluster_info: Parent cluster context.
            tissue: Optional tissue context.
            species: Optional species context.
            cell_type_hint: Optional parent-lineage hint, e.g. "macrophage".
            depth: Reserved for compatibility with the single-cluster API.
            max_reference_content_length: Maximum combined reference text sent
                to the synthesis call.
            max_brief_length: Maximum returned brief length.

        Returns:
            Dict with ``should_use_reference``, ``content``,
            ``references_used``, ``planning``, and ``tool_trace``.
        """
        del depth  # The agentic workflow reads selected documents directly.

        clean_marker_sets = _normalize_marker_sets(marker_sets)
        if not clean_marker_sets:
            return self._empty_result("No subcluster marker sets provided")

        overview_paths = self._select_overview_paths(cell_type_hint, major_cluster_info)
        overview_content = self._load_overviews(overview_paths)
        available_paths = self._list_available_reference_paths(overview_paths)

        tool_trace = [
            {"tool": "read_overview", "paths": overview_paths},
            {"tool": "list_reference_files", "paths": available_paths},
        ]

        planning = self._plan_subcluster_reference_reads(
            marker_sets=clean_marker_sets,
            major_cluster_info=major_cluster_info,
            tissue=tissue,
            species=species,
            cell_type_hint=cell_type_hint,
            overview_content=overview_content,
            available_paths=available_paths,
        )

        selected_paths = self._valid_selected_paths(
            planning.get("selected_references", []),
            available_paths,
        )
        if not selected_paths:
            selected_paths = self._fallback_reference_plan(clean_marker_sets, available_paths)
            planning.setdefault("reasoning", "")
            planning["reasoning"] = (
                planning["reasoning"] + " Fallback marker-overlap planning selected references."
            ).strip()
            planning["selected_references"] = selected_paths

        if not selected_paths:
            result = self._empty_result("No relevant detailed reference files selected")
            result["planning"] = planning
            result["tool_trace"] = tool_trace
            return result

        reference_documents = self._load_reference_documents(
            selected_paths,
            max_reference_content_length=max_reference_content_length,
        )
        tool_trace.append({"tool": "read_reference", "paths": selected_paths})

        if not reference_documents:
            result = self._empty_result("Selected references could not be loaded")
            result["planning"] = planning
            result["references_used"] = selected_paths
            result["tool_trace"] = tool_trace
            return result

        reference_marker_matches = self._build_reference_marker_match_summary(
            clean_marker_sets,
            selected_paths,
        )
        tool_trace.append({"tool": "reference_marker_lookup", "references": selected_paths})

        brief = self._synthesize_subcluster_reference_brief(
            marker_sets=clean_marker_sets,
            major_cluster_info=major_cluster_info,
            tissue=tissue,
            species=species,
            cell_type_hint=cell_type_hint,
            planning=planning,
            reference_documents=reference_documents,
            reference_marker_matches=reference_marker_matches,
            max_brief_length=max_brief_length,
        )
        tool_trace.append({"tool": "synthesize_reference_brief", "references": selected_paths})

        if not brief:
            result = self._empty_result("Reference brief synthesis returned no content")
            result["planning"] = planning
            result["references_used"] = selected_paths
            result["tool_trace"] = tool_trace
            return result

        return {
            "should_use_reference": True,
            "content": (
                "## Reference Marker Match Summary\n"
                f"{reference_marker_matches}\n\n"
                f"{brief}"
            ),
            "references_used": selected_paths,
            "preliminary_cell_type": planning.get("global_hypothesis", "Unknown"),
            "cell_type_range": [
                cluster.get("preliminary_pattern", "")
                for cluster in planning.get("cluster_hypotheses", [])
                if cluster.get("preliminary_pattern")
            ],
            "selected_references": selected_paths,
            "reasoning": planning.get("reasoning", ""),
            "planning": planning,
            "tool_trace": tool_trace,
        }

    def _empty_result(self, reason: str) -> Dict:
        return {
            "should_use_reference": False,
            "content": "",
            "references_used": [],
            "preliminary_cell_type": "Unknown",
            "cell_type_range": [],
            "selected_references": [],
            "reasoning": reason,
        }

    def _select_overview_paths(
        self,
        cell_type_hint: Optional[str],
        major_cluster_info: Optional[str],
    ) -> List[str]:
        context = f"{cell_type_hint or ''} {major_cluster_info or ''}".lower()
        if any(term in context for term in ["macrophage", "tam", "monocyte", "myeloid"]):
            return ["myeloid/macrophage/_overview.md", "myeloid/_overview.md"]
        if any(term in context for term in ["cd8", "cytotoxic t"]):
            return ["t_cell/cd8/_overview.md", "t_cell/_overview.md"]
        if any(term in context for term in ["cd4", "helper t", "treg"]):
            return ["t_cell/cd4/_overview.md", "t_cell/_overview.md"]
        if "t cell" in context or "t-cell" in context:
            return ["t_cell/_overview.md"]
        if "b cell" in context or "plasma" in context:
            return ["b_cell/_overview.md"]
        return ["_router.md"]

    def _load_overviews(self, overview_paths: List[str]) -> str:
        parts = []
        for path in overview_paths:
            content = self._load_reference_text(path)
            if content:
                parts.append(f"### Overview: {path}\n\n{content}")
        if parts:
            return format_reference_content("\n\n---\n\n".join(parts), 10000)

        router_content = self._load_reference_text("_router.md")
        return format_reference_content(router_content, 10000)

    def _list_available_reference_paths(self, overview_paths: List[str]) -> List[str]:
        roots = []
        for overview_path in overview_paths:
            if overview_path == "_router.md":
                continue
            root = str(Path(overview_path).parent)
            if root == ".":
                continue
            if root not in roots:
                roots.append(root)

        candidates = []
        if roots:
            for root in roots:
                root_path = self.reference_dir / root
                if root_path.exists():
                    candidates.extend(root_path.rglob("*.md"))
        else:
            candidates.extend(self.reference_dir.rglob("*.md"))

        paths = []
        for candidate in candidates:
            rel_path = candidate.relative_to(self.reference_dir).as_posix()
            if rel_path.endswith("_router.md") or rel_path.endswith("_overview.md"):
                continue
            paths.append(rel_path)

        return sorted(dict.fromkeys(paths))

    def _plan_subcluster_reference_reads(
        self,
        marker_sets: List[Dict],
        major_cluster_info: str,
        tissue: Optional[str],
        species: Optional[str],
        cell_type_hint: Optional[str],
        overview_content: str,
        available_paths: List[str],
    ) -> Dict:
        prompt = SUBCLUSTER_REFERENCE_PLAN_PROMPT.format(
            major_cluster_info=major_cluster_info,
            tissue=tissue or "Unknown",
            species=species or "Unknown",
            cell_type_hint=cell_type_hint or "",
            cluster_marker_text=_format_marker_sets(marker_sets),
            overview_content=overview_content or "No overview available.",
            available_reference_text="\n".join(f"- {path}" for path in available_paths) or "None",
        )

        try:
            response = call_llm(
                prompt=prompt,
                provider=self.provider,
                model=self.model,
                temperature=0,
                max_tokens=1400,
                api_key=self.api_key,
            )
            planning = _parse_json_response(response)
        except Exception as exc:
            planning = {
                "global_hypothesis": "Unknown",
                "selected_references": [],
                "cluster_hypotheses": [],
                "reasoning": f"Planning LLM call failed: {exc}",
            }

        if not isinstance(planning, dict):
            planning = {}

        planning.setdefault("global_hypothesis", "Unknown")
        planning.setdefault("selected_references", [])
        planning.setdefault("cluster_hypotheses", [])
        planning.setdefault("reasoning", "")
        return planning

    def _valid_selected_paths(self, selected_paths: List[str], available_paths: List[str]) -> List[str]:
        available_set = set(available_paths)
        valid = []
        for path in selected_paths or []:
            normalized = _normalize_reference_path(path)
            if normalized in available_set and normalized not in valid:
                valid.append(normalized)
        return valid[:4]

    def _fallback_reference_plan(self, marker_sets: List[Dict], available_paths: List[str]) -> List[str]:
        available_set = set(available_paths)
        markers = {
            str(marker).upper()
            for marker_set in marker_sets
            for marker in marker_set.get("markers", [])
        }
        selected = []

        def add(path: str):
            if path in available_set and path not in selected:
                selected.append(path)

        if markers & {
            "CD8A", "CD8B", "GZMK", "GZMA", "GZMB", "GZMH", "PRF1",
            "NKG7", "GNLY", "CCL5", "IFNG", "PDCD1", "LAG3", "HAVCR2",
            "TOX", "CXCL13", "ENTPD1", "TIGIT", "TCF7", "SLAMF6", "XCL1",
            "XCL2",
        }:
            add("t_cell/cd8/tumor_exhaustion_cytotoxic.md")
        if markers & {
            "CD4", "IL7R", "CCR7", "SELL", "TCF7", "LEF1", "ANXA1",
            "FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2", "CCR8",
            "TNFRSF18", "TNFRSF4", "ENTPD1", "LAYN", "CXCR5", "BCL6",
            "ICOS", "IL21", "CD40LG", "TOX2", "TBX21", "CXCR3", "CCR5",
            "RORC", "IL17A", "IL17F", "CCR6", "IL23R", "GATA3", "IL4",
            "IL5", "IL13", "CCR4",
        }:
            add("t_cell/cd4/helper_treg_tfh.md")
        if markers & {
            "KLRB1", "SLC4A10", "TRAV1-2", "ZBTB16", "RORA", "RORC",
            "IL18RAP", "CCR6", "CXCR6", "DPP4", "NCR3", "TRDC", "TRGC1",
            "TRGC2", "TRDV1", "TRDV2", "TRGV9", "TRGV10", "KLRD1", "NCAM1",
        }:
            add("t_cell/innate_like/mait_gamma_delta_nkt.md")
        if markers & {
            "CD69", "ITGAE", "CXCR6", "ITGA1", "ZNF683", "ISG15", "IFIT1",
            "IFIT2", "IFIT3", "IFITM1", "IFITM3", "MX1", "OAS1", "OAS2",
            "RSAD2", "STAT1", "MKI67", "TOP2A", "STMN1", "TYMS", "UBE2C",
            "PCLAF", "CENPF", "HSPA1A", "HSPA1B", "HSPA6", "DNAJB1",
            "FOS", "JUN",
        }:
            add("t_cell/state_modules/resident_ifn_cellcycle_stress.md")
        if markers & {
            "CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2",
            "CD4", "CD8A", "CD8B", "IL7R", "CCR7", "SELL", "TCF7",
            "LEF1", "GZMK", "GZMB", "PRF1", "NKG7", "GNLY", "PDCD1",
            "LAG3", "HAVCR2", "TOX", "CXCL13", "FOXP3", "IL2RA",
            "CTLA4", "TIGIT", "CCR8", "CXCR5", "BCL6", "ICOS", "CD69",
            "ITGAE", "CXCR6", "TRDC", "TRGC1", "TRGC2", "KLRB1",
            "SLC4A10", "MKI67", "TOP2A", "HSPA1A", "HSPA1B",
        }:
            add("t_cell/immune_t_cell_consensus.md")

        if markers & {
            "SPP1", "AREG", "EREG", "TREM2", "GPNMB", "VEGFA", "MMP9",
            "COL1A1", "COL1A2", "SPARC", "HMOX1", "SLC40A1",
        }:
            add("myeloid/macrophage/tam_pan_cancer.md")
        if markers & {
            "FOLR2", "SELENOP", "STAB1", "LYVE1", "MRC1", "CD163",
            "C1QA", "C1QB", "C1QC", "RNASE1",
        }:
            add("myeloid/macrophage/resident_like.md")
        if markers & {
            "CXCL9", "CXCL10", "GBP1", "GBP5", "STAT1", "ISG15",
            "IFIT1", "IFITM1", "IFITM3", "MX1", "IL1B", "TNF", "NLRP3",
        }:
            add("myeloid/macrophage/inflammatory_interferon.md")
        return selected[:4]

    def _load_reference_documents(
        self,
        selected_paths: List[str],
        max_reference_content_length: int,
    ) -> str:
        parts = []
        per_doc_limit = max(2000, max_reference_content_length // max(1, len(selected_paths)))
        for path in selected_paths:
            content = self._load_reference_text(path)
            if not content:
                continue
            parts.append(
                f"### Reference Document: {path}\n\n"
                f"{format_reference_content(content, per_doc_limit)}"
            )
        return format_reference_content("\n\n---\n\n".join(parts), max_reference_content_length)

    def _build_reference_marker_match_summary(
        self,
        marker_sets: List[Dict],
        selected_paths: List[str],
    ) -> str:
        if "myeloid/macrophage/tam_pan_cancer.md" not in selected_paths:
            return "No structured reference-marker lookup was available for the selected references."

        lines = [
            "Structured lookup against the human cancer macrophage consensus layer.",
            "Use the consensus program as the primary label. Use paper-specific aliases only as traceability evidence.",
        ]
        reference_sets = {
            state["id"]: {
                "state": state,
                "markers": {marker.upper() for marker in state["markers"]},
            }
            for state in HUMAN_CANCER_MACROPHAGE_CONSENSUS_STATES
        }
        for marker_set in marker_sets:
            cluster_id = marker_set.get("cluster_id", "unknown")
            query_markers = {str(marker).upper() for marker in marker_set.get("markers", [])}
            ranked = []
            for state_id, reference in reference_sets.items():
                overlap = sorted(query_markers & reference["markers"])
                if overlap:
                    ranked.append((len(overlap), state_id, reference["state"], overlap))
            confidence_rank = {"high": 0, "medium-high": 1, "medium": 2, "ambiguous": 3}
            ranked.sort(key=lambda item: (
                -item[0],
                confidence_rank.get(item[2].get("confidence"), 9),
                item[2]["consensus_name"],
            ))
            if not ranked:
                lines.append(f"- Cluster {cluster_id}: no consensus macrophage/TAM program hit.")
                continue

            best = ranked[0]
            alternatives = ranked[1:3]
            alt_text = ""
            if alternatives:
                alt_text = "; alternatives: " + "; ".join(
                    f"{state['consensus_name']} ({count} markers: {', '.join(overlap[:6])})"
                    for count, _state_id, state, overlap in alternatives
                )
            aliases = ", ".join(best[2].get("aliases", [])[:3]) or "none"
            lines.append(
                f"- Cluster {cluster_id}: best consensus program {best[2]['consensus_name']} "
                f"(confidence: {best[2]['confidence']}; {best[0]} markers: {', '.join(best[3][:10])}; "
                f"paper aliases: {aliases}){alt_text}."
            )
        return "\n".join(lines)

    def _synthesize_subcluster_reference_brief(
        self,
        marker_sets: List[Dict],
        major_cluster_info: str,
        tissue: Optional[str],
        species: Optional[str],
        cell_type_hint: Optional[str],
        planning: Dict,
        reference_documents: str,
        reference_marker_matches: str,
        max_brief_length: Optional[int],
    ) -> str:
        prompt = SUBCLUSTER_REFERENCE_BRIEF_PROMPT.format(
            major_cluster_info=major_cluster_info,
            tissue=tissue or "Unknown",
            species=species or "Unknown",
            cell_type_hint=cell_type_hint or "",
            cluster_marker_text=_format_marker_sets(marker_sets),
            planning_json=json.dumps(planning, ensure_ascii=False, indent=2),
            reference_documents=reference_documents,
            reference_marker_matches=reference_marker_matches,
        )

        try:
            brief = call_llm(
                prompt=prompt,
                provider=self.provider,
                model=self.model,
                temperature=0,
                max_tokens=2500,
                api_key=self.api_key,
            )
        except Exception:
            return self._fallback_reference_brief(marker_sets, planning, reference_documents)

        return format_reference_content(brief, max_brief_length)

    def _fallback_reference_brief(
        self,
        marker_sets: List[Dict],
        planning: Dict,
        reference_documents: str,
    ) -> str:
        facts = _extract_source_facts(reference_documents)
        cluster_lines = []
        for marker_set in marker_sets:
            cluster_id = marker_set.get("cluster_id", "unknown")
            markers = ", ".join(marker_set.get("markers", [])[:12])
            cluster_lines.append(
                f"- Cluster {cluster_id}: review selected lineage/state programs against markers {markers}."
            )

        fact_text = "\n".join(f"- {fact}" for fact in facts[:6]) or "- Reference documents were selected, but no source facts could be extracted automatically."
        return format_reference_content(
            "<reference_brief>\n"
            "## Reference Agent Summary\n"
            f"{planning.get('global_hypothesis', 'Subtype reference documents were selected for this run.')}\n\n"
            "## Objective Literature Facts\n"
            f"{fact_text}\n\n"
            "## Cluster-Specific Guidance\n"
            f"{chr(10).join(cluster_lines)}\n\n"
            "## Cross-Cluster Distinctions\n"
            "- Prefer specific marker-program labels over broad lineage-only labels when evidence supports them.\n"
            "</reference_brief>"
        )

    def _load_reference_text(self, reference_path: str) -> str:
        normalized = _normalize_reference_path(reference_path)
        path = self.reference_dir / normalized
        if not path.exists():
            return ""
        return load_markdown_file(path)

    def list_available_references(self, category: Optional[str] = None) -> List[Dict]:
        if category:
            try:
                from .reference_selector import find_references_by_category
            except ImportError:
                from reference_selector import find_references_by_category
            return find_references_by_category(category)

        results = []
        if self.index and "references" in self.index:
            for ref_id, ref_info in self.index["references"].items():
                results.append(
                    {
                        "id": ref_id,
                        "category": ref_info.get("category", ""),
                        "cell_types": ref_info.get("cell_types", []),
                        "path": ref_info.get("path", ""),
                    }
                )
        return results

    def get_reference_content_direct(
        self,
        reference_id: str,
        section_path: Optional[str] = None,
    ) -> Optional[str]:
        if not self.index or "references" not in self.index:
            return None

        ref_info = self.index["references"].get(reference_id)
        if not ref_info:
            return None

        ref_path = self.reference_dir / ref_info.get("path", "")
        if not ref_path.exists():
            return None

        if section_path:
            try:
                from .section_extractor import get_section_by_path
            except ImportError:
                from section_extractor import get_section_by_path
            return get_section_by_path(ref_path, section_path)

        try:
            from .utils import load_markdown_file, parse_yaml_frontmatter
        except ImportError:
            from utils import load_markdown_file, parse_yaml_frontmatter
        content = load_markdown_file(ref_path)
        _, body = parse_yaml_frontmatter(content)
        return body


def get_reference_content(
    markers: List[str],
    tissue: Optional[str] = None,
    species: Optional[str] = None,
    cell_type_hint: Optional[str] = None,
    provider: str = "openrouter",
    model: Optional[str] = None,
    **kwargs,
) -> Dict:
    """Convenience wrapper around ReferenceAgent.get_reference_for_markers."""
    agent = ReferenceAgent(provider=provider, model=model)
    return agent.get_reference_for_markers(
        markers=markers,
        tissue=tissue,
        species=species,
        cell_type_hint=cell_type_hint,
        **kwargs,
    )


def get_subcluster_reference_brief(
    marker_sets: List[Dict],
    major_cluster_info: str,
    tissue: Optional[str] = None,
    species: Optional[str] = None,
    cell_type_hint: Optional[str] = None,
    provider: str = "openrouter",
    model: Optional[str] = None,
    **kwargs,
) -> Dict:
    """Convenience wrapper for agentic subcluster reference briefing."""
    agent = ReferenceAgent(provider=provider, model=model)
    return agent.get_reference_brief_for_subclusters(
        marker_sets=marker_sets,
        major_cluster_info=major_cluster_info,
        tissue=tissue,
        species=species,
        cell_type_hint=cell_type_hint,
        **kwargs,
    )


def format_reference_for_prompt(reference_result: Dict) -> str:
    """Format a reference_result for injection into the annotation prompt."""
    if not reference_result.get("should_use_reference") or not reference_result.get("content"):
        return ""

    return f"""<expert_reference>
The following expert-curated reference information is provided to assist with cell type annotation.
Use this technical guidance to inform your analysis, particularly for subtype differentiation
and marker interpretation.

Preliminary assessment: {reference_result.get('preliminary_cell_type', 'Unknown')}

{reference_result['content']}
</expert_reference>"""


def _normalize_marker_sets(marker_sets: List[Dict]) -> List[Dict]:
    normalized = []
    for marker_set in marker_sets or []:
        cluster_id = marker_set.get("cluster_id", marker_set.get("id", "unknown"))
        markers = marker_set.get("markers", [])
        if isinstance(markers, str):
            raw_markers = re.split(r"[,;|\n\s]+", markers)
        else:
            raw_markers = list(markers or [])

        clean_markers = []
        seen = set()
        for marker in raw_markers:
            marker = str(marker).strip().strip("'\"`")
            if not marker:
                continue
            key = marker.upper()
            if key in seen:
                continue
            seen.add(key)
            clean_markers.append(marker)
        if clean_markers:
            normalized.append({"cluster_id": str(cluster_id), "markers": clean_markers})
    return normalized


def _format_marker_sets(marker_sets: List[Dict]) -> str:
    lines = []
    for marker_set in marker_sets:
        markers = ", ".join(marker_set.get("markers", [])[:30])
        lines.append(f"- Cluster {marker_set.get('cluster_id')}: {markers}")
    return "\n".join(lines)


def _parse_json_response(response: str) -> Dict:
    if not response:
        return {}
    json_match = re.search(r"\{[\s\S]*\}", str(response))
    if not json_match:
        return {}
    try:
        return json.loads(json_match.group())
    except json.JSONDecodeError:
        return {}


def _normalize_reference_path(path: str) -> str:
    normalized = (path or "").strip().lstrip("/").lstrip("\\")
    if normalized.startswith("references/"):
        normalized = normalized[len("references/"):]
    if normalized.startswith("references_brain/"):
        normalized = normalized[len("references_brain/"):]
    return normalized.replace("\\", "/")


def _extract_source_facts(reference_documents: str) -> List[str]:
    facts = []
    for line in reference_documents.splitlines():
        stripped = line.strip(" -*")
        if not stripped:
            continue
        if "doi.org/" in stripped or " et al." in stripped:
            facts.append(stripped)
    return list(dict.fromkeys(facts))
