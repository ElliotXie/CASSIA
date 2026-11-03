"""Deterministic HTML report for direct Fused Boost annotations."""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

import pandas as pd


def _escape(value: Any, fallback: str = "Not available") -> str:
    if value is None or value == "":
        return fallback
    if isinstance(value, (dict, list, tuple)):
        value = json.dumps(value, ensure_ascii=False, indent=2)
    return html.escape(str(value)).replace("\n", "<br>")


def _items(values: Any) -> list[str]:
    if values is None or values == "":
        return []
    if isinstance(values, str):
        values = values.split(",")
    if not isinstance(values, (list, tuple, set)):
        values = [values]
    return [str(value).strip() for value in values if str(value).strip()]


def _chips(values: Any, css_class: str = "") -> str:
    items = _items(values)
    if not items:
        return '<span class="muted">None</span>'
    cls = f"chip {css_class}".strip()
    return "".join(f'<span class="{cls}">{html.escape(item)}</span>' for item in items)


def _assistant_rounds(messages: Sequence[Mapping[str, str]]) -> list[str]:
    return [str(message.get("content", "")) for message in messages if message.get("role") == "assistant"]


def _remove_gene_request(text: str) -> str:
    import re

    return re.sub(
        r"<check_genes>\s*[\s\S]*?\s*</check_genes>",
        "",
        text,
        flags=re.IGNORECASE,
    ).strip()


def _query_table(frame: pd.DataFrame) -> str:
    if frame.empty:
        return '<p class="muted">No marker statistics were returned.</p>'
    columns = [
        column for column in ("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "status")
        if column in frame.columns
    ]
    if not columns:
        return '<p class="muted">No displayable marker statistics were returned.</p>'
    head = "".join(f"<th>{html.escape(column)}</th>" for column in columns)
    rows = []
    for _, row in frame[columns].iterrows():
        cells = []
        for column in columns:
            value = row[column]
            text = _format_table_value(value, column)
            cells.append(f"<td>{html.escape(text)}</td>")
        rows.append(f"<tr>{''.join(cells)}</tr>")
    return f'<div class="table-wrap"><table><thead><tr>{head}</tr></thead><tbody>{"".join(rows)}</tbody></table></div>'


def _format_table_value(value: Any, column: str) -> str:
    """Format report-only marker statistics without changing saved raw values."""
    if pd.isna(value):
        return ""
    if column in {"gene", "status"}:
        return str(value)
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if column in {"p_val", "p_val_adj"} and number != 0 and abs(number) < 0.001:
        return f"{number:.2e}"
    return f"{number:.3f}".rstrip("0").rstrip(".")


def _quality_card(assessment: Optional[Mapping[str, Any]]) -> str:
    if not assessment:
        return ""
    score = assessment.get("score")
    try:
        score_number = max(0.0, min(100.0, float(score)))
    except (TypeError, ValueError):
        score_number = 0.0
    risk = str(assessment.get("risk") or ("low" if score_number >= 80 else "medium" if score_number >= 60 else "high"))
    return f"""
    <section class="card quality-card">
      <div class="section-heading"><span class="eyebrow">OPTIONAL BLIND REVIEW</span><h2>Annotation quality assessment</h2></div>
      <div class="quality-grid">
        <div class="score-ring" style="--score:{score_number:.0f}"><strong>{score_number:.0f}</strong><span>/100</span></div>
        <div><span class="risk risk-{html.escape(risk.lower())}">{html.escape(risk.title())} risk</span>
        <p>{_escape(assessment.get('reasoning'))}</p>
        <div class="chip-row">{_chips(assessment.get('flags'), 'warning')}</div></div>
      </div>
      <p class="fine-print">This score is answer-blind and is not a ground-truth accuracy guarantee.</p>
    </section>"""


def write_fused_boost_html_report(
    output_path: Path,
    manifest: Mapping[str, Any],
    messages: Sequence[Mapping[str, str]],
    result: Optional[Mapping[str, Any]],
    errors: Sequence[str],
    query_frames: Sequence[pd.DataFrame],
    top_marker_rows: pd.DataFrame,
    quality_assessment: Optional[Mapping[str, Any]] = None,
) -> Path:
    """Render Fused Boost as a primary annotation, never as a revision report."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result = dict(result or {})
    status = str(manifest.get("status", "unknown"))
    confidence = str(result.get("confidence", "unknown"))
    ranked = _items(result.get("ranked_sub_cell_types") or result.get("sub_cell_types"))
    if result.get("final_sub_cell_type") and str(result["final_sub_cell_type"]) not in ranked:
        ranked.insert(0, str(result["final_sub_cell_type"]))
    top_genes = top_marker_rows.get("gene", pd.Series(dtype=str)).dropna().astype(str).head(30).tolist()
    rounds = _assistant_rounds(messages)
    timeline = []
    for index, frame in enumerate(query_frames, start=1):
        reasoning = _remove_gene_request(rounds[index - 1]) if index <= len(rounds) else ""
        genes = frame.get("gene", pd.Series(dtype=str)).dropna().astype(str).tolist()
        timeline.append(f"""
        <details class="round" {'open' if index == 1 else ''}>
          <summary><span>Evidence round {index}</span><span class="round-count">{len(genes)} genes checked</span></summary>
          <div class="round-body"><h3>Working hypotheses</h3><p>{_escape(reasoning, 'The agent requested targeted marker evidence.')}</p>
          <h3>Requested genes</h3><div class="chip-row">{_chips(genes)}</div>{_query_table(frame)}</div>
        </details>""")
    usage = manifest.get("usage") or {}
    reported_cost = manifest.get("reported_cost_per_cluster")
    backend_name = str(manifest.get("provider") or manifest.get("backend") or "")
    if reported_cost is not None:
        cost_text = str(reported_cost)
    elif backend_name == "codex-cli":
        cost_text = "Codex subscription — no per-run cost reported"
    elif backend_name.endswith("-cli") or backend_name == "shell":
        cost_text = "No per-run cost reported by CLI backend"
    else:
        cost_text = "Not reported by provider"
    error_html = ""
    if errors:
        error_html = f'<section class="card error"><h2>Run errors</h2><ul>{"".join(f"<li>{html.escape(str(error))}</li>" for error in errors)}</ul></section>'
    ranked_html = "".join(
        f'<li><span class="rank">{index}</span><span>{html.escape(label)}</span></li>'
        for index, label in enumerate(ranked[:5], start=1)
    ) or '<li class="muted">No ranked subtype was returned.</li>'

    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>CASSIA Fused Boost v2 · Cluster {html.escape(str(manifest.get('cluster', '')))}</title>
<style>
:root{{--ink:#17201d;--muted:#64706b;--paper:#f4f1ea;--card:#fffdf8;--line:#d9d4c8;--green:#176b55;--green2:#dcebe4;--amber:#a35c15;--red:#a63c36;--shadow:0 16px 44px rgba(38,48,43,.09)}}
*{{box-sizing:border-box}} body{{margin:0;background:linear-gradient(135deg,#eef3ed 0,#f8f4ec 45%,#eee8dd 100%);color:var(--ink);font:15px/1.65 Inter,ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}}
.shell{{max-width:1180px;margin:auto;padding:44px 24px 72px}} .hero{{background:#173f36;color:#f8f5ea;border-radius:24px;padding:34px 38px;box-shadow:var(--shadow);position:relative;overflow:hidden}}
.hero:after{{content:"";position:absolute;width:260px;height:260px;border:1px solid rgba(255,255,255,.18);border-radius:50%;right:-60px;top:-130px;box-shadow:0 0 0 42px rgba(255,255,255,.04),0 0 0 84px rgba(255,255,255,.025)}}
.eyebrow{{font-size:11px;letter-spacing:.16em;font-weight:800;color:#8fcab7}} h1{{font:700 clamp(30px,5vw,52px)/1.05 Georgia,serif;margin:10px 0 14px;max-width:800px}} .hero p{{max-width:760px;color:#dbe8e1;margin:0}}
.meta{{display:flex;flex-wrap:wrap;gap:9px;margin-top:24px}} .meta span,.badge{{border:1px solid rgba(255,255,255,.24);border-radius:999px;padding:6px 11px;font-size:12px;background:rgba(255,255,255,.08)}}
.layout{{display:grid;grid-template-columns:minmax(0,1.55fr) minmax(280px,.75fr);gap:20px;margin-top:20px}} .stack{{display:grid;gap:20px;align-content:start}}
.card{{background:rgba(255,253,248,.96);border:1px solid var(--line);border-radius:20px;padding:25px;box-shadow:var(--shadow)}} .section-heading h2,.card h2{{font:700 24px/1.2 Georgia,serif;margin:5px 0 18px}} .call{{font:700 clamp(27px,4vw,42px)/1.15 Georgia,serif;color:#134b3e;margin:5px 0 8px}}
.subtype{{font-size:18px;margin:0 0 18px}} .facts{{display:grid;grid-template-columns:repeat(3,1fr);gap:10px;margin-top:18px}} .fact{{padding:12px;background:#f1f4ef;border-radius:12px}} .fact b{{display:block;font-size:11px;text-transform:uppercase;color:var(--muted);letter-spacing:.08em}} .fact span{{font-weight:700}}
.chip-row{{display:flex;flex-wrap:wrap;gap:7px;margin:8px 0}} .chip{{background:var(--green2);color:#174c40;border-radius:999px;padding:5px 9px;font-size:12px}} .chip.negative{{background:#f5dfdc;color:#7d302c}} .chip.warning{{background:#f8ead3;color:#794818}}
.ranked{{list-style:none;padding:0;margin:0;display:grid;gap:9px}} .ranked li{{display:flex;gap:11px;align-items:center;padding:10px 12px;border:1px solid var(--line);border-radius:12px}} .rank{{display:grid;place-items:center;width:25px;height:25px;border-radius:50%;background:#173f36;color:white;font-size:12px;font-weight:800}}
.round{{border:1px solid var(--line);border-radius:14px;background:#fff;margin:10px 0;overflow:hidden}} .round summary{{cursor:pointer;display:flex;justify-content:space-between;gap:15px;padding:15px 17px;font-weight:800}} .round-count{{font-size:12px;color:var(--muted);font-weight:600}} .round-body{{border-top:1px solid var(--line);padding:17px}} .round-body h3{{font-size:13px;text-transform:uppercase;letter-spacing:.08em;color:var(--muted);margin:15px 0 6px}}
.table-wrap{{overflow:auto;margin-top:15px}} table{{border-collapse:collapse;width:100%;font-size:12px}} th,td{{text-align:left;border-bottom:1px solid #e7e3da;padding:8px;white-space:nowrap}} th{{color:var(--muted);font-weight:800;background:#f7f5ef}}
.provenance dl{{display:grid;grid-template-columns:auto 1fr;gap:8px 14px;margin:0}} dt{{color:var(--muted)}} dd{{margin:0;font-weight:650;overflow-wrap:anywhere}} .muted{{color:var(--muted)}} .fine-print{{font-size:12px;color:var(--muted)}} .error{{border-color:#e4aaa5;background:#fff5f4}}
.quality-grid{{display:grid;grid-template-columns:105px 1fr;gap:18px;align-items:center}} .score-ring{{width:96px;height:96px;border-radius:50%;display:grid;place-content:center;text-align:center;background:conic-gradient(var(--green) calc(var(--score)*1%),#e2e0d8 0);position:relative}} .score-ring:after{{content:"";position:absolute;inset:9px;background:var(--card);border-radius:50%}} .score-ring strong,.score-ring span{{z-index:1}} .score-ring strong{{font-size:28px;line-height:1}} .risk{{display:inline-block;padding:5px 9px;border-radius:999px;font-weight:800;font-size:12px;background:#eee}} .risk-high{{background:#f5dfdc;color:var(--red)}} .risk-medium{{background:#f8ead3;color:var(--amber)}} .risk-low{{background:var(--green2);color:var(--green)}}
@media(max-width:800px){{.shell{{padding:18px 12px 45px}}.hero{{padding:27px 23px;border-radius:18px}}.layout{{grid-template-columns:1fr}}.facts{{grid-template-columns:1fr 1fr}}.quality-grid{{grid-template-columns:1fr}}}}
</style></head><body><main class="shell">
<header class="hero"><span class="eyebrow">DIRECT ACTIVE-EVIDENCE ANNOTATION</span><h1>CASSIA Fused Boost v2</h1><p>A primary cell-type annotation built from ranked markers, explicit competing hypotheses, and targeted local evidence queries. No prior annotation is assumed or revised.</p>
<div class="meta"><span>Cluster {html.escape(str(manifest.get('cluster', '')))}</span><span>{html.escape(status.title())}</span><span>{len(query_frames)} evidence rounds</span><span>{len(_items(manifest.get('checked_genes')))} genes checked</span></div></header>
<div class="layout"><div class="stack">
<section class="card"><div class="section-heading"><span class="eyebrow">PRIMARY CALL</span><h2>Final annotation</h2></div><div class="call">{_escape(result.get('final_cell_type') or result.get('main_cell_type'), 'No final annotation')}</div><p class="subtype">{_escape(result.get('final_sub_cell_type'), 'No subtype/state returned')}</p>
<div class="facts"><div class="fact"><b>Confidence</b><span>{html.escape(confidence.title())}</span></div><div class="fact"><b>Mode</b><span>Fused Boost v2</span></div><div class="fact"><b>Schema</b><span>{_escape(result.get('schema_version') or manifest.get('result_schema_version'))}</span></div></div></section>
{_quality_card(quality_assessment)}
<section class="card"><div class="section-heading"><span class="eyebrow">EVIDENCE SYNTHESIS</span><h2>Why this call</h2></div><p>{_escape(result.get('evidence'))}</p><h3>Supporting markers</h3><div class="chip-row">{_chips(result.get('supporting_markers'))}</div><h3>Refuting markers</h3><div class="chip-row">{_chips(result.get('refuting_markers'), 'negative')}</div></section>
<section class="card"><div class="section-heading"><span class="eyebrow">ACTIVE EVIDENCE LOG</span><h2>Hypotheses and marker queries</h2></div>{''.join(timeline) or '<p class="muted">No completed evidence round.</p>'}</section>
<section class="card"><div class="section-heading"><span class="eyebrow">NEXT STEP</span><h2>Recommended validation</h2></div><p>{_escape(result.get('recommended_next_steps'), 'No additional validation was requested.')}</p></section>{error_html}
</div><aside class="stack">
<section class="card"><h2>Ranked alternatives</h2><ol class="ranked">{ranked_html}</ol><h3>Other alternatives</h3><div class="chip-row">{_chips(result.get('alternatives'))}</div><h3>Possible mixed identities</h3><div class="chip-row">{_chips(result.get('possible_mixed_cell_types'), 'warning')}</div></section>
<section class="card"><h2>Input marker context</h2><div class="chip-row">{_chips(top_genes)}</div><p class="fine-print">Top 30 ranked genes shown. The active evidence table records every targeted local query used before finalization.</p></section>
<section class="card provenance"><h2>Run provenance</h2><dl><dt>Provider</dt><dd>{_escape(manifest.get('provider') or manifest.get('backend'))}</dd><dt>Model</dt><dd>{_escape(manifest.get('model'))}</dd><dt>LLM calls</dt><dd>{_escape(manifest.get('llm_calls'))}</dd><dt>Total tokens</dt><dd>{_escape(usage.get('total_tokens'), 'Not captured in this run')}</dd><dt>Billing</dt><dd>{html.escape(cost_text)}</dd><dt>Execution</dt><dd>{_escape(manifest.get('execution_time'))} s</dd></dl></section>
</aside></div><p class="fine-print">Generated deterministically from saved CASSIA artifacts. Report formatting does not alter the annotation.</p></main></body></html>"""
    output_path.write_text(document, encoding="utf-8")
    return output_path


__all__ = ["write_fused_boost_html_report"]
