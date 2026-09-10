# Boost-100 findings: no-cap fused Annotation Boost

Final run: `cursor-agent / composer-2.5`, 100 fixed cases, no numerical gene cap.

## Outcome

| Arm | Exact | Correct + partial | Wrong | Missing |
|---|---:|---:|---:|---:|
| Historical one-shot output, freshly re-judged | 50 | 71 | 27 | 2 |
| One-shot then no-cap review Boost | 56 | 75 | 25 | 0 |
| No-baseline fused no-cap Boost | 65 | 76 | 24 | 0 |

Fused versus one-shot produced 22 exact improvements and 7 exact regressions
(two-sided exact McNemar p=0.00813). Fused versus post-annotation Boost produced
15 exact improvements and 6 regressions (p=0.07835). Correct + partial is nearly
tied between fused and post-Boost (76 versus 75).

The practical reading is that fused Boost improves exact subtype/state resolution,
but it is not a universally safer annotator. On the 50 basic controls it retained
43 exact answers, while it recovered 16/35 boost-tier and 6/15 frontier cases.
Buried-marker cases were the clearest win (10/15 exact); RNA/protein-discordant
cases remained unsolved (0/4 exact).

## Cost and evidence use

- Fused: 256 agent calls, 156 marker-query rounds, 27.4 checked genes/case on average.
- One-shot then Boost: 376 total agent calls, 177 Boost query rounds, 30.0 checked genes/case.
- Fused used about one-third fewer calls and input/output tokens than the two-stage route.
- At Cursor's published Composer 2.5 Standard rates, one-shot is approximately $2.10,
  one-shot then Boost is $5.18, fused Boost is $3.40, and the final three-arm judge is $0.61.
  These are token-equivalent on-demand values; subscription First-party-pool inclusion can
  cover some or all of the actual out-of-pocket charge.
- Per cluster, that is $0.0212 for one-shot, $0.0309 for the incremental post-Boost step,
  $0.0521 for the full one-shot-then-Boost route, and $0.0340 for fused Boost. The judge
  costs $0.0061 per benchmark cluster for all three arms together, or $0.0020 per annotation
  output judged. The one-shot mean uses its 99 available raw envelopes.
- Counting the discarded first judge pass and four parser-fix reruns, the new benchmark
  execution consumed $7.90 in token-equivalent usage. Including historical one-shot
  generation makes the fully accounted total $10.00. The phase ledger is saved separately.

## Why answers can still become wrong

Marker querying supplies more evidence, but the model chooses the hypotheses and the
genes. A wrong early hypothesis can therefore request a self-confirming panel and become
more confidently over-specific. The main regressions were sibling compartments or states:
DCT versus connecting tubule, choroid plexus versus ependymal, erythroid-primed progenitor
versus HSC, CD56dim versus CD56bright NK, macula densa versus parietal epithelium, and
neuroendocrine cell versus neuro-mimetic mTEC. Spatial paper labels and RNA evidence can
also genuinely disagree.

## Infrastructure audit

The fused prompt keeps the original `final_annotation_system_v1` text intact and appends
an active-evidence extension; the opening necessarily changes because there is no trusted
prior annotation. Every fused case completed at least one orchestrated query. Agent calls
ran in an isolated empty workspace in Cursor ask mode, so Composer could not inspect the
marker CSV directly.

Removing the cap exposed a parser edge case: explanatory prose containing a literal
`<check_genes>` mention could be mistaken for an opening request tag. The parser now rejects
nested angle brackets and accepts only gene-symbol-shaped tokens. All four affected cases
were rerun; final stored gene lists match the real request tags in all 200 Boost results.

The old 20-gene-capped post-Boost run remains a historical artifact. It is not used as the
primary control here because this comparison requires the same no-cap policy. Its higher
historical score should be treated as a hypothesis that a cap or other run differences may
regularize review Boost, not as a causal cap result without a dedicated randomized rerun.
