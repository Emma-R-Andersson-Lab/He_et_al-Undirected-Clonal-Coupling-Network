# Known issues

The analysis and plotting scripts are shipped **as they were run for the
paper**, so that published results stay byte-reproducible. The issues below were
found by review and confirmed by measurement or by reproducing the failure. None
of them affect the correctness of the published numbers; they are recorded here
so that anyone reusing or modifying the code knows what to expect.

The single intentional deviation from the scripts as run is that
`--min-shared-clones` now defaults to `10` in step 2 rather than `3`, matching
the criterion stated in the paper. This changes only the contents of the
intermediate `pairwise_coupling_significant_positive_edges.csv`; the final edge
set is unchanged, because the published run applied the threshold of 10 at the
plotting step instead. See the README.

---

## What was verified as correct

Before the issues, the parts that were checked and hold up:

- The coupling score matches a brute-force transcription of
  $S_{xy}=\sum_c I_{c,xy}(C_{cx}+C_{cy})/n_c$ to within 8.9e-16.
- The permutation preserves exactly the three quantities it claims: clone sizes,
  per-embryo cluster counts, and embryo membership. Clone sizes being preserved
  is what keeps the precomputed $1/n_c$ weights valid across permutations.
- `bh_fdr` is numerically identical to `statsmodels.stats.multitest.multipletests(method="fdr_bh")`.
- The empirical $p$-value convention $(b+1)/(R+1)$ matches the paper.
- On a null dataset with no planted coupling, $p$-values are uniform, $z \sim N(0,1)$,
  and no edge survives FDR. On a dataset with three planted fate modules, the
  pipeline recovers all nine within-module edges with zero false positives and
  correctly isolates the unplanted clusters.
- The "backbone" is a genuine maximum spanning **forest**: on a disconnected
  graph it spans each component separately rather than forcing a single tree.

---

## 1. About 88% of the permutation runtime is wasted (performance)

**Severity:** high impact on runtime, no impact on results.

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

Over 10,000 permutations on one core that is roughly **13 minutes versus 1
minute**. The fix is to split the function so the permutation loop computes only
the score, and to use the fact that the score matrix is symmetric:

```python
def coupling_score_only(counts, inverse_clone_sizes):
    presence = (counts > 0).astype(np.float64)
    weighted = counts.astype(np.float64) * inverse_clone_sizes[:, None]
    m = weighted.T @ presence
    return m + m.T          # identical to weighted.T @ presence + presence.T @ weighted
```

`shared_clones` is only ever needed once, for the observed data.

A distant second bottleneck is `np.add.at` in `counts_from_assignments()`; a
flattened `np.bincount` is roughly 4x faster there.

---

## 2. Permutation count silently censors edges when there are many clusters

**Severity:** can suppress real findings. Read this before running on a
dataset with many clusters.

The smallest attainable $p$-value is $1/(R+1)$; with the default 10,000
permutations that floor is 1.0e-4. Benjamini-Hochberg then multiplies by the
number of pairs, $m = K(K-1)/2$. For a single maximally-significant edge:

| Clusters | Pairs | $q$ for one lone maximal edge | Edges needed at the floor to reach $q \le 0.05$ |
| --- | --- | --- | --- |
| 10 | 45 | 0.0045 | 1 |
| 20 | 190 | 0.0190 | 1 |
| 30 | 435 | 0.0435 | 1 |
| 40 | 780 | **0.0780** | 2 |
| 60 | 1770 | **0.1770** | 4 |
| 80 | 3160 | **0.3160** | 7 |
| 100 | 4950 | **0.4950** | 10 |

From 40 clusters upward, an isolated strong edge can fail FDR purely because the
permutation count is too low for the number of tests — not because the coupling
is weak. Nothing warns you; the edge simply does not appear.

Rule of thumb: use `--permutations` of at least $K(K-1)/2 \div 0.05$. That is
about 15,600 for 40 clusters and 35,400 for 60.

---

## 3. `--min-z` below 0 crashes at the plotting stage

**Severity:** low; wastes a run.

`proportional_edge_width_map()` raises `ValueError` on any non-positive edge.
This happens *after* the full permutation sweep and after every CSV has been
written, so the tables survive but no figures are produced:

```
ValueError: Graph contains a non-positive plotted edge: 'cluster_0'–'cluster_9', z=-0.39
```

It is only reachable if `--fdr` is loosened at the same time: a negative $z$
implies $p_\text{enrichment} > 0.5$, so such an edge cannot pass $q \le 0.05$ on
its own. Reproduced with `--min-z -5 --fdr 1.0`. If you want to explore
sub-threshold edges, read `pairwise_coupling_all.csv` rather than loosening both
flags.

---

## 4. Infinite z-scores are handled in one place out of three

**Severity:** low; needs a degenerate null, but fails silently.

`safe_zscore()` deliberately emits `±inf` when the permutation null has zero
variance. `proportional_edge_width_map()` clamps that for edge widths, but two
other consumers of the same value do not:

- `nx.spring_layout(..., weight="weight")`
- `nx.maximum_spanning_tree(..., weight="weight")`

A single infinite edge weight turns **every** layout coordinate into `NaN`,
producing a garbage figure rather than an error. If you hit the
"Non-finite z-score for edge" warning, clamp the weight before layout.

The plotting script has the same gap: `component_layout()` compresses weights
with `math.log1p(z)`, and `log1p(inf)` is still `inf`.

---

## 5. Docstring statements that contradict the code

**Severity:** documentation only. `docs/methods.md` is correct; the module
docstring in `02_clonal_coupling_network.py` is not, in two places.

- *"Negative scores are retained in the complete output matrices"* — the
  observed score cannot be negative. It is a sum of $(C_{cx}+C_{cy})/n_c$ over
  clones present in both clusters, so every term is strictly positive. Only
  **z-scores** go negative. (This is also why `log2_enrichment` never produces a
  `NaN` from a negative argument.)
- The informal FDR formula $q = p \times (m / \text{rank})$ omits the monotonic
  step-up enforcement and the clip to 1. The **code** performs full standard
  Benjamini-Hochberg, verified identical to `statsmodels`. The raw ratio as
  written in the docstring is not BH and would give non-monotonic $q$-values.
  Cite the description in `docs/methods.md`, not the docstring.

---

## 6. Minor

- `proportional_edge_width_map()` is recomputed for each of the three output
  formats, so the non-finite-z warning can fire three times for one run.
- `build_graphs()` is called twice in `main()`; the second call exists only to
  attach the degree metrics to the exported node attributes. Harmless, just
  redundant.
- The null variance uses the sum-of-squares form
  $\sum x^2 - n\bar{x}^2$, which is prone to catastrophic cancellation in
  general. At the magnitudes these scores actually take it is fine, and the
  result is clamped at zero before the square root.

---

## Not a defect, but worth stating

FDR correction is applied across all $K(K-1)/2$ cluster pairs **before** the
`--min-shared-clones` filter. This is conservative — the filter does not reduce
the number of tests — and is a deliberate choice, but reviewers reasonably ask
about it. It is documented in section 5 of `docs/methods.md`.
