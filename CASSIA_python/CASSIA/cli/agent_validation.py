"""Validator orchestration for subscription-backed CASSIA annotation.

The functions in this module deliberately do not know how an agent is
executed.  Callers provide a callback, which lets the production CLI and the
benchmark harness share the exact same validation/revision mechanism.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Sequence

from CASSIA.engine.main_function_code import (
    coupling_validator_system_v0,
    coupling_validator_system_v1,
    coupling_validator_system_v2,
)


AgentCall = Callable[[str, int, str], str]


SELF_REFLECT_VALIDATOR_SYSTEM_V2 = """
You are an independent, skeptical expert reviewer of a single-cell RNA-seq cell-type annotation.
Your job is not to check whether the proposed story sounds plausible. Your job is to try to
falsify the proposed final annotation and decide whether another interpretation is better
supported by the ranked marker list.

Perform the following review in order:

1. Independent reconstruction: Before accepting the proposed conclusion, infer the broad
   lineage and the three strongest cell-type/subtype hypotheses directly from the ranked
   markers and tissue context.
2. Evidence audit: Separate strong identity markers from state, stress, housekeeping,
   ribosomal, mitochondrial, ambient-RNA, and spatial-mixture signals. Never treat one
   isolated marker as decisive when the corresponding marker program is absent.
3. Head-to-head challenge: Compare the proposed final annotation against the strongest
   alternative. State the positive markers for each, contradictory markers, and important
   expected markers that are missing. Absence is soft evidence because of dropout, but a
   subtype claim without any positive defining markers is not adequately supported.
4. Resolution check: Validate the exact final subtype, not merely its broad lineage. Fail an
   over-specific or sibling-subtype call when a different subtype is better supported. If the
   evidence only supports a broader label, require the annotation to be broadened and its
   uncertainty stated.
5. Context challenge: Use species and tissue as priors, not as proof. Explicitly examine
   contradictions between markers and tissue. For spatial data or suspected contamination,
   distinguish the likely source cell from ambient or neighboring-cell transcripts instead of
   explaining away inconvenient markers without evidence.
6. Decision: Pass only if the proposed final annotation is the best-supported hypothesis at
   the claimed resolution and no alternative is materially better supported. Internal
   consistency alone is not enough. When evidence is genuinely tied, fail and request a
   broader or explicitly uncertain annotation.

Do not invent genes that are not in the supplied ranked marker list. Do not defer to the
annotation's confidence or rhetoric.

Output exactly one decision line first:
Validation result: VALIDATION PASSED
or
Validation result: VALIDATION FAILED

Then provide a compact review with these fields:
Independent top candidates: ...
Evidence for proposed call: ...
Strongest counter-evidence / alternative: ...
Resolution and context check: ...
Revision instruction: ...

For a failed result, the revision instruction must name the most defensible replacement or
the broader label that should be used, and say what uncertainty must be acknowledged.
"""


@dataclass
class ValidatedAnnotationRun:
    """Result of the CASSIA annotation-validator revision loop."""

    final_response: str
    validation_passed: bool
    validation_attempts: int
    history: List[Dict[str, object]]


def is_tissue_blind(tissue: Optional[str]) -> bool:
    """Return whether CASSIA should use its tissue-blind validator prompt."""
    return not tissue or tissue.strip().lower() in {"none", "tissue blind"}


def select_validator_system(tissue: Optional[str], involvement: str = "v1") -> str:
    """Select the same validator system prompt used by the original engine."""
    if involvement == "v0":
        return coupling_validator_system_v0.strip()
    return (
        coupling_validator_system_v2.strip()
        if is_tissue_blind(tissue)
        else coupling_validator_system_v1.strip()
    )


def build_validator_prompt(
    annotation_response: str,
    marker_list: Sequence[str],
    tissue: Optional[str],
    additional_info: Optional[str] = None,
    involvement: str = "v1",
    species: Optional[str] = None,
) -> str:
    """Build the original CASSIA coupling-validator system + user message.

    Agent CLIs generally accept one prompt rather than distinct system and user
    messages, so the two original messages are concatenated without changing
    their wording.
    """
    cleaned_markers = [
        str(marker).strip() for marker in marker_list if str(marker).strip()
    ]
    if involvement == "self-reflect-v2":
        system = SELF_REFLECT_VALIDATOR_SYSTEM_V2.strip()
        markers = "\n".join(
            f"{rank}. {marker}" for rank, marker in enumerate(cleaned_markers, start=1)
        )
        validation_message = f"""Review the following annotation independently and adversarially.

Proposed Annotation:
{annotation_response}

Dataset Context:
Species: {species or 'Not specified'}
Tissue: {tissue or 'Not specified'}
Additional Info: {additional_info or 'None'}

Ranked Marker List (highest to lowest):
{markers}

Apply the full falsification and head-to-head comparison protocol. Judge the exact primary
final call, not whether it appears somewhere among several alternatives.
"""
        return f"{system}\n\n{validation_message}"

    system = select_validator_system(tissue, involvement)
    markers = ", ".join(cleaned_markers)
    validation_message = f"""Please validate the following annotation result:

Annotation Result:
{annotation_response}

Context:

Marker List: {markers}
Additional Info: {additional_info or 'None'}

Validate the annotation based on this context.
"""
    return f"{system}\n\n{validation_message}"


def build_revision_prompt(
    original_prompt: str,
    previous_response: str,
    validation_feedback: str,
    involvement: str = "v1",
    marker_list: Optional[Sequence[str]] = None,
    tissue: Optional[str] = None,
    additional_info: Optional[str] = None,
    species: Optional[str] = None,
) -> str:
    """Build the same retry message used by the original CASSIA engine."""
    if involvement == "self-reflect-v2":
        markers = ", ".join(
            str(marker).strip() for marker in (marker_list or []) if str(marker).strip()
        )
        return f"""Your previous annotation was challenged by an independent reviewer.
Do a fresh annotation from the evidence; do not merely defend or cosmetically edit the old
answer. Compare the previous call with the reviewer's strongest alternative, then choose the
best-supported final call at the resolution justified by the markers. If subtype evidence is
insufficient, use a broader label and state the uncertainty. Do not mention the validation
process in the final report.

Previous response:
{previous_response}

Independent review:
{validation_feedback}

Evidence context:
Species: {species or 'Not specified'}
Tissue: {tissue or 'Not specified'}
Ranked markers: {markers or 'See original prompt'}
Additional info: {additional_info or 'None'}

Original annotation task:
{original_prompt}

Return one complete updated annotation, including a single unambiguous primary final call."""
    return f"""Previous annotation attempt failed validation. Please review your previous response and the validation feedback, then provide an updated annotation:

Previous response:
{previous_response}

Validation feedback:
{validation_feedback}

Original prompt:
{original_prompt}

Please provide an updated annotation addressing the validation feedback."""


def validation_passed(response: str) -> bool:
    """Match the pass condition used by the original CASSIA engine."""
    return "VALIDATION PASSED" in response


def run_validated_annotation(
    initial_prompt: str,
    marker_list: Sequence[str],
    tissue: Optional[str],
    call_agent: AgentCall,
    additional_info: Optional[str] = None,
    involvement: str = "v1",
    max_attempts: int = 3,
    initial_response: Optional[str] = None,
    species: Optional[str] = None,
) -> ValidatedAnnotationRun:
    """Run annotation -> validation -> revision, up to ``max_attempts``.

    ``call_agent`` receives ``(stage, attempt_number, prompt)`` where stage is
    either ``"annotation"`` or ``"validation"``.  Annotation attempt 1 uses
    ``initial_prompt`` verbatim, which is important for paired benchmarks that
    reuse prompts saved by an earlier one-shot run.  When ``initial_response``
    is supplied, that saved one-shot response is reused as annotation attempt
    1 instead of making a new model call.  This makes the validator benchmark
    strictly paired at the point where validation begins.
    """
    if max_attempts < 1:
        raise ValueError("max_attempts must be at least 1")

    history: List[Dict[str, object]] = []
    current_prompt = initial_prompt
    final_response = ""

    for attempt in range(1, max_attempts + 1):
        reused = attempt == 1 and initial_response is not None
        annotation_response = (
            initial_response
            if reused
            else call_agent("annotation", attempt, current_prompt)
        )
        final_response = annotation_response
        history.append({
            "attempt": attempt,
            "stage": "annotation",
            "prompt": current_prompt,
            "response": annotation_response,
            "reused": reused,
        })

        validator_prompt = build_validator_prompt(
            annotation_response=annotation_response,
            marker_list=marker_list,
            tissue=tissue,
            additional_info=additional_info,
            involvement=involvement,
            species=species,
        )
        validator_response = call_agent("validation", attempt, validator_prompt)
        passed = validation_passed(validator_response)
        history.append({
            "attempt": attempt,
            "stage": "validation",
            "prompt": validator_prompt,
            "response": validator_response,
            "passed": passed,
        })

        if passed:
            return ValidatedAnnotationRun(
                final_response=final_response,
                validation_passed=True,
                validation_attempts=attempt,
                history=history,
            )

        if attempt < max_attempts:
            current_prompt = build_revision_prompt(
                original_prompt=initial_prompt,
                previous_response=annotation_response,
                validation_feedback=validator_response,
                involvement=involvement,
                marker_list=marker_list,
                tissue=tissue,
                additional_info=additional_info,
                species=species,
            )

    return ValidatedAnnotationRun(
        final_response=final_response,
        validation_passed=False,
        validation_attempts=max_attempts,
        history=history,
    )
