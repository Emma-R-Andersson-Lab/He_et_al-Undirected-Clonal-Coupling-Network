# Known issues

The analysis and plotting scripts are shipped **as they were run for the
paper**, so published results stay reproducible. The issues below were found by
review and confirmed by measurement or by reproducing the failure.

**None of them invalidate the published results.** Section 0 records the checks
that establish this against the real dataset.

The single intentional deviation from the scripts as run is that
`--min-shared-clones` now defaults to `10` in step 2 rather than `3`, matching
the criterion stated in the paper. This changes only the contents of the
intermediate `pairwise_coupling_significant_positive_edges.csv` (172 rows to 158
rows); the plotted edge set is unchanged, because the published run applied the
threshold of 10 at the plotting step instead.

---

## 0. Verification against the published dataset

The published run covers 102,077 barcode-positive cells, 8,315 embryo-specific
clones, 45 clusters and 9 embryos, giving 990 unique cluster pairs at 10,000
permutations.

Independently recomputing the method's equations from the raw
`lineage_cells.csv` and comparing to the shipped result matrices:

| Check | Result |
| --- | --- |
| Coupling score $S_{xy}$ recomputed from the paper formula vs `observed_coupling_score_matrix.csv` | max abs diff **6.8e-13** |
| Shared-clone counts recomputed vs `shared_clone_count_matrix.csv` | **exactly identical** |
| $z$ vs $(S-\mu)/\sigma$ from the reported columns | max abs diff **1.4e-14** |
| $q$ vs Benjamini-Hochberg recomputed from reported $p$ | max abs diff **2.2e-16** |
| Score matrix symmetric, and non-negative everywhere | **yes** (min 0.0) |
| Non-finite $z$-scores anywhere | **0** |
| Zero-variance nulls anywhere | **0** |

Separately, on synthetic data: a null dataset with no planted coupling gives
uniform $p$-values, $z \sim N(0,1)$ and no surviving edge; a dataset with three
planted fate modules recovers all nine within-module edges with zero false
positives and correctly isolates the unplanted clusters. `bh_fdr` is numerically
identical to `statsmodels`' `fdr_bh`. The permutation preserves exactly the three
quantities it claims, and the "backbone" is a genuine maximum spanning **forest**,
spanning each component separately rather than forcing a single tree.

The method description in `docs/methods.md` is therefore an accurate account of
what the code does, on this dataset.

---

## 1. About 88% of the permutation runtime is wasted (performance)

**Severity:** high impact on runtime, **no impact on results**.

`coupling_score_from_counts()` also builds the `shared_clones` matrix, but the
permutation loop discards it:

```python
permuted_score, _ = coupling_score_from_counts(permuted_counts, inverse_clone_sizes)
```

That discarded line is an `int64` matrix multiply, which gets no BLAS path and
is far slower than the `float64` score it accompanies.

Measured on 100,000 cells / 20,000 clones / 30 clusters, per permutation:

| Variant | Time |
| --- | --- |
| As shipped | 77.7 ms |
| Without the discarded `shared_clones` | 9.6 ms |
| Also exploiting symmetry (one matmul) | 6.1 ms |

On the real published dataset the shipped code costs 81.8 ms per permutation per
core, i.e. ~14 min for 10,000 permutations on one core and ~1 min on 16. The fix
would bring the single-core figure to roughly 1 minute:

```python
def coupling_score_only(counts, inverse_clone_sizes):
    presence = (counts > 0).astype(np.float64)
    weighted = counts.astype(np.float64) * inverse_clone_sizes[:, None]
    m = weighted.T @ presence
    return m + m.T          # identical to weighted.T @ presence + presence.T @ weighted
```

`shared_clones` is only ever needed once, for the observed data. A distant second
bottleneck is `np.add.at` in `counts_from_assignments()`; a flattened
`np.bincount` is roughly 4x faster there.

---

## 2. Permutation count can silently censor edges when there are many clusters

**Severity:** can suppress real findings in general. **Did not affect the
published result** — see the verdict below.

The smallest attainable $p$-value is $1/(R+1)$; with 10,000 permutations that
floor is 1.0e-4. Benjamini-Hochberg then multiplies by the number of pairs,
$m = K(K-1)/2$. For a single maximally-significant edge:

| Clusters | Pairs | $q$ for one lone maximal edge | Edges needed at the floor to reach $q \le 0.05$ |
| --- | --- | --- | --- |
| 10 | 45 | 0.0045 | 1 |
| 20 | 190 | 0.0190 | 1 |
| 30 | 435 | 0.0435 | 1 |
| 40 | 780 | **0.0780** | 2 |
| **45** | **990** | **0.0990** | **2** |
| 60 | 1770 | **0.1770** | 4 |
| 100 | 4950 | **0.4950** | 10 |

Rule of thumb: use `--permutations` of at least $K(K-1)/2 \div 0.05$ — about
19,800 for the 45 clusters here, 15,600 for 40, 35,400 for 60.

**Verdict for the published run.** The condition is *live* at $K = 45$: a single
isolated maximal edge would have received $q = 0.099$ and been dropped. It did
not bite, because **122 edges sat at the $p$-floor simultaneously**, and BH
divides by rank: those edges received $q = 9.999\text{e-}5 \times 990 / 122 =
8.1\text{e-}4$, three orders of magnitude inside the threshold. Observed range
across the 172 significant edges was $q \in [8.1\text{e-}4,\ 0.0351]$, all
comfortably below 0.05. No edge was censored by permutation resolution.

The caveat matters for anyone re-running on a **sparser** dataset, where few
edges reach the floor. There, raise `--permutations`.

---

## 3. `--min-z` below 0 crashes at the plotting stage

**Severity:** low; wastes a run. Not used for the paper.

`proportional_edge_width_map()` raises `ValueError` on any non-positive edge,
*after* the full permutation sweep and after every CSV has been written, so the
tables survive but no figures are produced:

```
ValueError: Graph contains a non-positive plotted edge: 'cluster_0'–'cluster_9', z=-0.39
```

Only reachable if `--fdr` is loosened at the same time: a negative $z$ implies
$p_\text{enrichment} > 0.5$, so such an edge cannot pass $q \le 0.05$ on its own.
Reproduced with `--min-z -5 --fdr 1.0`. To explore sub-threshold edges, read
`pairwise_coupling_all.csv` instead of loosening both flags.

---

## 4. Infinite z-scores are handled in one place out of three

**Severity:** low; needs a degenerate null. **Did not occur in the published
run** (0 non-finite $z$, 0 zero-variance nulls).

`safe_zscore()` deliberately emits `±inf` when the permutation null has zero
variance. `proportional_edge_width_map()` clamps that for edge widths, but two
other consumers of the same value do not:

- `nx.spring_layout(..., weight="weight")`
- `nx.maximum_spanning_tree(..., weight="weight")`

A single infinite edge weight turns **every** layout coordinate into `NaN`,
producing a garbage figure rather than an error. If you ever see the
"Non-finite z-score for edge" warning, clamp the weight before layout. The
plotting script has the same gap: `component_layout()` compresses weights with
`math.log1p(z)`, and `log1p(inf)` is still `inf`.

---

## 5. Docstring statements that contradict the code

**Severity:** documentation only. `docs/methods.md` is correct; the module
docstring in `02_clonal_coupling_network.py` is not, in two places.

- *"Negative scores are retained in the complete output matrices"* — the
  observed score cannot be negative. It is a sum of $(C_{cx}+C_{cy})/n_c$ over
  clones present in both clusters, so every term is strictly positive; confirmed
  on the real data (minimum 0.0). Only **z-scores** go negative.
- The informal FDR formula $q = p \times (m / \text{rank})$ omits the monotonic
  step-up enforcement and the clip to 1. The **code** performs full standard
  Benjamini-Hochberg, verified identical to `statsmodels` and to a direct
  recomputation from the published $p$-values. The raw ratio as written in the
  docstring is not BH. Cite `docs/methods.md`, which states it correctly.

The docstring also still says "shared clones between 2 clusters >= 10" where the
original said 3; that line was updated alongside the default.

---

## 6. `env_activate.txt` does not describe the published figure

**Severity:** documentation only, but actively misleading.

The scratch notes file kept alongside the original scripts records invocations of
`plot_clonal_coupling_components.py` and
`plot_clonal_coupling_components_clustered.py` with `--layout-mode community`.
Neither produced the published figure:

- `_clustered.py` draws **straight** edges and has a single combined legend; the
  published figure has curved edges and two separate legends.
- The published figure was made with the `_curved_with_legends` variant using
  **`--layout-mode spring`**, `--node-size-column n_unique_clones`,
  `--min-shared-clones 10` and `--node-labels-only`.

This was established by re-rendering both layouts from the shipped results:
`spring` reproduced the published figure's topology, node sizes and both legends
exactly; `community` produced a visibly different three-blob layout roughly twice
as wide. The correct command is in README section 5.

---

## 7. The node-size legend is hard-coded to say "cells"

**Severity:** affects the published figure caption.

In `plot_component()` the node-size legend entries are built as:

```python
label=f"{count:,} cells"
```

regardless of what `--node-size-column` points at. The published figure sized
nodes by `n_unique_clones`, so its legend reads "88 cells / 424 cells /
2,377 cells" when the values are in fact **clone counts**, not cell counts.

The footer normally names the column in use ("Node area uses
`'n_unique_clones'` on one global scale"), but `--node-labels-only` suppresses
the footer, so nothing on the figure discloses it.

Fix the legend text during figure assembly, or state the unit explicitly in the
figure caption. A code fix is a one-line change to that f-string.

---

## 8. Minor

- `proportional_edge_width_map()` is recomputed for each of the three output
  formats, so the non-finite-z warning can fire three times for one run.
- `build_graphs()` is called twice in `main()`; the second call exists only to
  attach degree metrics to the exported node attributes. Harmless, redundant.
- The null variance uses the sum-of-squares form $\sum x^2 - n\bar{x}^2$, prone
  to catastrophic cancellation in general. At the magnitudes these scores take it
  is fine, and the result is clamped at zero before the square root; the
  published $z$-values reproduce to 1.4e-14.

---

## Not a defect, but worth stating

FDR correction is applied across all $K(K-1)/2$ cluster pairs **before** the
`--min-shared-clones` filter. This is conservative — the filter does not reduce
the number of tests — and deliberate, but reviewers reasonably ask. Documented in
section 5 of `docs/methods.md`.
