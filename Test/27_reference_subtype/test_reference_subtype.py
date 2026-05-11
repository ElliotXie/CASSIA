"""
CASSIA Test 27: Reference Subtype Module
========================================
Focused non-live tests for macrophage subtype reference retrieval and
subclustering prompt injection.
"""

import json
import importlib.util
import sys
import tempfile
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "shared" / "python"))

from test_utils import setup_cassia_imports, print_test_header, print_test_result

setup_cassia_imports()


def test_reference_agent_selects_macrophage_docs():
    import CASSIA.agents.reference_agent.complexity_scorer as complexity_scorer
    from CASSIA.agents.reference_agent import ReferenceAgent

    original_call_llm = complexity_scorer.call_llm

    def fake_call_llm(**kwargs):
        return json.dumps({
            "preliminary_cell_type": "SPP1/AREG tumor-associated macrophage",
            "cell_type_range": [
                "SPP1/AREG inflammatory angiogenic TAM",
                "APOE/TREM2 lipid-associated TAM",
            ],
            "selected_references": [
                "myeloid/macrophage/tam_pan_cancer.md",
                "myeloid/macrophage/_overview.md",
            ],
            "reasoning": "SPP1, AREG, CXCL8, and IL1B support TAM subtype reference use.",
        })

    complexity_scorer.call_llm = fake_call_llm
    try:
        agent = ReferenceAgent(provider="openrouter", model="moonshotai/kimi-k2.6")
        result = agent.get_reference_for_markers(
            markers=["SPP1", "AREG", "CXCL8", "IL1B", "TIMP1", "LYZ", "C1QA"],
            tissue="tumor",
            species="human",
            cell_type_hint="macrophage",
        )
    finally:
        complexity_scorer.call_llm = original_call_llm

    assert result["should_use_reference"] is True
    assert "myeloid/macrophage/tam_pan_cancer.md" in result["references_used"]
    assert "SPP1/AREG inflammatory angiogenic TAM" in result["content"]
    assert "M1/M2" in result["content"]


def test_reference_agent_builds_agentic_subcluster_brief():
    import CASSIA.agents.reference_agent.reference_agent as reference_agent_module
    from CASSIA.agents.reference_agent import get_subcluster_reference_brief

    original_call_llm = reference_agent_module.call_llm
    calls = []

    def fake_call_llm(**kwargs):
        prompt = kwargs["prompt"]
        calls.append(prompt)
        if "## Reference Documents Read" in prompt:
            assert "myeloid/macrophage/tam_pan_cancer.md" in prompt
            return """<reference_brief>
## Reference Agent Summary
The marker landscape fits tumor macrophage states.

## Objective Literature Facts
- Cheng et al. Cell 2021 reported SPP1+ TAMs as angiogenesis-associated.
- Coulton et al. Nat Commun 2024 describes IFNGMac marked by CXCL9/CXCL10/GBP genes.

## Cluster-Specific Guidance
- Cluster spp1_tam: likely SPP1/AREG inflammatory angiogenic TAM.
- Cluster ifng_tam: likely IFNG/CXCL9 interferon-activated TAM.

## Cross-Cluster Distinctions
- Separate SPP1/AREG inflammatory TAM from IFNG/CXCL9 TAM by chemokine and GBP markers.
</reference_brief>"""

        return json.dumps({
            "global_hypothesis": "macrophage TAM subtype landscape",
            "selected_references": [
                "myeloid/macrophage/tam_pan_cancer.md",
                "myeloid/macrophage/inflammatory_interferon.md",
            ],
            "cluster_hypotheses": [
                {
                    "cluster_id": "spp1_tam",
                    "preliminary_pattern": "SPP1/AREG TAM",
                    "supporting_markers": ["SPP1", "AREG"],
                    "reference_targets": ["myeloid/macrophage/tam_pan_cancer.md"],
                },
                {
                    "cluster_id": "ifng_tam",
                    "preliminary_pattern": "IFNG/CXCL9 TAM",
                    "supporting_markers": ["CXCL9", "CXCL10", "GBP1"],
                    "reference_targets": ["myeloid/macrophage/inflammatory_interferon.md"],
                },
            ],
            "reasoning": "Macrophage overview points to TAM and IFN subtype references.",
        })

    reference_agent_module.call_llm = fake_call_llm
    try:
        result = get_subcluster_reference_brief(
            marker_sets=[
                {"cluster_id": "spp1_tam", "markers": ["SPP1", "AREG", "CXCL8", "LYZ"]},
                {"cluster_id": "ifng_tam", "markers": ["CXCL9", "CXCL10", "GBP1", "LYZ"]},
            ],
            major_cluster_info="human tumor macrophage",
            tissue="tumor",
            species="human",
            cell_type_hint="macrophage",
            provider="openrouter",
            model="moonshotai/kimi-k2.6",
        )
    finally:
        reference_agent_module.call_llm = original_call_llm

    assert result["should_use_reference"] is True
    assert len(calls) == 2
    assert "macrophage/_overview.md" in calls[0]
    assert "myeloid/macrophage/tam_pan_cancer.md" in result["references_used"]
    assert "Objective Literature Facts" in result["content"]
    assert "Cheng et al. Cell 2021" in result["content"]
    assert result["tool_trace"][0]["tool"] == "read_overview"


def test_subclustering_reference_context_is_injected():
    import CASSIA.agents.subclustering.subclustering as subclustering

    marker_df = pd.DataFrame({
        "cluster": ["spp1_tam", "resident"],
        "markers": [
            "SPP1, AREG, CXCL8, IL1B, TIMP1, LYZ, C1QA",
            "FOLR2, SELENOP, SLC40A1, C1QA, C1QB, C1QC, CD163",
        ],
    })

    class DummyReferenceAgent:
        def __init__(self, provider, model):
            self.provider = provider
            self.model = model

        def get_reference_brief_for_subclusters(self, marker_sets, **kwargs):
            return {
                "should_use_reference": True,
                "content": """<reference_brief>
## Reference Agent Summary
This run contains SPP1/AREG TAM and FOLR2 resident-like macrophage patterns.

## Objective Literature Facts
- Cheng et al. Cell 2021 reported SPP1+ TAMs as angiogenesis-associated.

## Cluster-Specific Guidance
- Cluster spp1_tam: SPP1/AREG inflammatory angiogenic TAM.
- Cluster resident: FOLR2/SELENOP resident-like macrophage.

## Cross-Cluster Distinctions
- Separate inflammatory angiogenic TAM from resident-like macrophage by AREG/CXCL8 versus FOLR2/SELENOP.
</reference_brief>""",
                "references_used": [
                    "myeloid/macrophage/tam_pan_cancer.md",
                    "myeloid/macrophage/resident_like.md",
                ],
                "planning": {
                    "cluster_hypotheses": [
                        {"cluster_id": item["cluster_id"], "preliminary_pattern": "macrophage subtype"}
                        for item in marker_sets
                    ]
                },
                "tool_trace": [{"tool": "read_overview"}],
                "reasoning": "Selected macrophage subtype references.",
            }

        def get_reference_for_markers(self, markers, **kwargs):
            if "SPP1" in markers:
                return {
                    "should_use_reference": True,
                    "content": "### Reference: tam_pan_cancer\nSPP1/AREG inflammatory angiogenic TAM",
                    "references_used": ["myeloid/macrophage/tam_pan_cancer.md"],
                    "preliminary_cell_type": "Macrophage",
                    "cell_type_range": ["TAM"],
                    "reasoning": "SPP1 TAM markers",
                }
            return {
                "should_use_reference": True,
                "content": "### Reference: resident_like\nFOLR2/SELENOP resident-like macrophage",
                "references_used": ["myeloid/macrophage/resident_like.md"],
                "preliminary_cell_type": "Macrophage",
                "cell_type_range": ["resident-like macrophage"],
                "reasoning": "Resident-like markers",
            }

    original_reference_agent = subclustering.ReferenceAgent
    subclustering.ReferenceAgent = DummyReferenceAgent
    try:
        context, info = subclustering.build_subcluster_reference_context(
            marker=marker_df,
            major_cluster_info="human tumor macrophage",
            provider="openrouter",
            reference_model="moonshotai/kimi-k2.6",
            reference_cell_type_hint="macrophage",
            tissue="tumor",
            species="human",
        )
    finally:
        subclustering.ReferenceAgent = original_reference_agent

    assert info["reference_used"] is True
    assert info["references_used"] == [
        "myeloid/macrophage/tam_pan_cancer.md",
        "myeloid/macrophage/resident_like.md",
    ]
    assert "spp1_tam" in context
    assert "Objective Literature Facts" in context
    assert "SPP1/AREG inflammatory angiogenic TAM" in context
    assert "FOLR2/SELENOP resident-like macrophage" in context


def test_subcluster_prompt_contains_reference_block():
    from CASSIA.agents.subclustering.subclustering import construct_prompt_from_csv_subcluster

    marker_df = pd.DataFrame({
        "cluster": ["ifng_tam"],
        "markers": ["CXCL9, CXCL10, GBP1, GBP5, STAT1, LYZ, C1QA"],
    })

    prompt = construct_prompt_from_csv_subcluster(
        marker=marker_df,
        major_cluster_info="human tumor macrophage",
        additional_context="<expert_reference>IFNG/CXCL9 interferon-activated macrophage</expert_reference>",
    )

    assert "Cluster ifng_tam" in prompt
    assert "<expert_reference>" in prompt
    assert "IFNG/CXCL9 interferon-activated macrophage" in prompt


def test_llm_usage_tracking_records_cost():
    import CASSIA.core.llm_utils as llm_utils

    llm_utils.reset_llm_usage_log()
    llm_utils._record_llm_usage(
        provider="openrouter",
        model="moonshotai/kimi-k2.6",
        response_id="gen-test",
        usage={
            "prompt_tokens": 100,
            "completion_tokens": 25,
            "total_tokens": 125,
            "cost": 0.00125,
            "completion_tokens_details": {"reasoning_tokens": 5},
            "prompt_tokens_details": {"cached_tokens": 10},
        },
    )

    summary = llm_utils.get_llm_usage_summary()
    log = llm_utils.get_llm_usage_log()
    assert summary["requests"] == 1
    assert summary["prompt_tokens"] == 100
    assert summary["completion_tokens"] == 25
    assert summary["reasoning_tokens"] == 5
    assert summary["cached_tokens"] == 10
    assert summary["cost"] == 0.00125
    assert log[0]["response_id"] == "gen-test"

    llm_utils.reset_llm_usage_log()


def test_benchmark_loads_external_case_csv():
    benchmark_path = (
        Path(__file__).resolve().parents[2]
        / "Benchmark"
        / "reference_subtype"
        / "macrophage_subtype_benchmark.py"
    )
    spec = importlib.util.spec_from_file_location("macrophage_subtype_benchmark", benchmark_path)
    benchmark = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(benchmark)

    with tempfile.TemporaryDirectory() as tmpdir:
        cases_csv = Path(tmpdir) / "cases.csv"
        pd.DataFrame([{
            "case_id": "paper_cluster_1",
            "markers": "CXCL9, CXCL10, GBP1, LYZ",
            "expected_terms": "IFNGMac/IFN-gamma;CXCL9;macrophage/TAM",
            "min_score": 2,
        }]).to_csv(cases_csv, index=False)

        cases = benchmark.load_cases(str(cases_csv))
        marker_df = benchmark.marker_dataframe(cases)

    assert cases[0]["id"] == "paper_cluster_1"
    assert cases[0]["expected_terms"] == [
        ["IFNGMac", "IFN-gamma"],
        ["CXCL9"],
        ["macrophage", "TAM"],
    ]
    assert marker_df.loc[0, "cluster"] == "paper_cluster_1"
    assert "CXCL10" in marker_df.loc[0, "markers"]
    assert benchmark.paper_label_hit("This matches Coulton 8_IFNGMac", "8_IFNGMac")


def run_reference_subtype_tests():
    print_test_header("27 - Reference Subtype Module")
    tests = [
        test_reference_agent_selects_macrophage_docs,
        test_reference_agent_builds_agentic_subcluster_brief,
        test_subclustering_reference_context_is_injected,
        test_subcluster_prompt_contains_reference_block,
        test_llm_usage_tracking_records_cost,
        test_benchmark_loads_external_case_csv,
    ]
    errors = []

    for test in tests:
        try:
            test()
            print(f"  PASS: {test.__name__}")
        except Exception as exc:
            errors.append(f"{test.__name__}: {exc}")
            print(f"  FAIL: {test.__name__}: {exc}")

    success = not errors
    print_test_result(success, f"{len(tests) - len(errors)}/{len(tests)} tests passed")
    return success


if __name__ == "__main__":
    sys.exit(0 if run_reference_subtype_tests() else 1)
