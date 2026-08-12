# Changes from the original working scripts

A complete record of every difference between the scripts as they were run for
the paper and the copies in this repository. The intent is that the analysis
code stays exactly as it worked.

Verify any claim below yourself with `diff` and `md5sum`.

---

## `scripts/03_plot_clonal_coupling_components.py` — UNCHANGED

**Byte-identical** to the original `plot_clonal_coupling_components_curved_with_legends.py`.

```
d2eaaa3a771d8fa3a3c06e95e7bd20d0   original
d2eaaa3a771d8fa3a3c06e95e7bd20d0   scripts/03_plot_clonal_coupling_components.py
```

Renamed only. No content change of any kind.

---

## `scripts/02_clonal_coupling_network.py` — 3 LINES CHANGED

Original `clonal_coupling_network.py` md5 `488ad91c124940a7c2294d05ffbffb7f`;
repository copy md5 `7ba101da881e48453f0027f1e23571d5`.

The complete diff is:

```diff
@@ line 78 (module docstring) @@
-shared clones between 2 clusters >= 3
+shared clones between 2 clusters >= 10

@@ lines 212-213 (argument definition) @@
-        default=3,
-        help="Minimum observed shared clones required for an edge. Default: 3.",
+        default=10,
+        help="Minimum observed shared clones required for an edge. Default: 10.",
```

Nothing else differs. No function body, no numerical code, no output format.

**Why.** The paper states the edge criterion as ≥ 10 shared clones. The original
default of 3 meant `pairwise_coupling_significant_positive_edges.csv` contained
172 edges, of which 14 fell below the stated threshold; the criterion was
enforced later, by passing `--min-shared-clones 10` to the plotting step.

**Effect on results: none.** The plotted network is identical either way. Both
routes yield the same 158 edges across the same 4 components:

| Route | significant-edges CSV | plotted edges |
| --- | --- | --- |
| Original: analysis `--min-shared-clones 3`, plot `--min-shared-clones 10` | 172 rows | 158 |
| This repo: analysis `--min-shared-clones 10`, plot default 0 | 158 rows | 158 |

Confirmed directly against the published results: of the 172 edges in the
shipped CSV, exactly 158 have `shared_clones >= 10`, and the published
`connected_component_manifest.csv` totals 132 + 24 + 1 + 1 = 158 edges.

`pairwise_coupling_all.csv` — the supplementary table containing every pair — is
**not affected at all**, since it is written before any threshold is applied.

**To restore the original behaviour exactly**, either pass
`--min-shared-clones 3` on the command line, or revert those three lines. No
other change is needed.

---

## `scripts/01_export_lineage_input.R` — REWRITTEN (new file)

This is the one file that is not a copy. The original was an RMarkdown chunk
that assumed a Seurat object was already loaded in the session and that several
libraries were already attached. It was restructured into a script so it can be
`source()`d interactively *or* run headless with `Rscript`.

**The data logic is preserved.** Cell table construction, the invalid-clone
value list, the duplicate-cell-ID check, `clone_uid = paste(embryo, clone_id,
sep = "::")`, the `all_clustered_cells` / `lineage_cells` split, factor-level
cluster ordering, and all five summary tables are unchanged in behaviour.

Structural differences:

| Change | Rationale |
| --- | --- |
| Wrapped in a function `export_lineage_input(obj, cluster_col, clone_col, embryo_col, outdir, ...)` | Callable and testable; no reliance on globals. |
| Added `library()` calls | The chunk relied on packages attached elsewhere in the notebook. |
| Added a `--rds` command-line entry point, guarded by `sys.nframe() == 0 && !interactive()` | Allows headless runs; verified that `source()` does not trigger it. |
| `readr::write_csv` replaced by a `write_csv_compat()` helper | Uses `readr` when installed, falls back to `utils::write.csv(row.names = FALSE, na = "")` otherwise. Removes a hard dependency. Both produce CSVs the Python step parses identically. |
| Cluster-ordering branches consolidated; `setdiff()` result now sorted | The original appended unmatched clusters in encounter order. Sorting makes output deterministic. Affects only column/row ordering, never values. |
| Output filenames derived from a table list rather than five separate `write_csv` calls | Same five filenames, same contents. |
| Hard-coded `outdir` and metadata column names replaced by parameters | The originals were user-edited constants at the top of the chunk. |

**Testing.** Verified against a synthetic Seurat 5.5.1 object: factor level order
is carried into `cluster_summary.csv`, `unassigned`/`NA` clones are dropped,
`clone_uid` never spans embryos, and the resulting `lineage_cells.csv` is
accepted by step 2 unmodified.

If you would rather ship the original chunk verbatim, it can be added as
`scripts/01_export_lineage_input.Rmd` alongside this script.

---

## Files added (no effect on analysis)

`README.md`, `KNOWN_ISSUES.md`, `CHANGES_FROM_ORIGINAL.md`, `docs/methods.md`,
`LICENSE`, `CITATION.cff`, `requirements.txt`, `environment.yml`, `.gitignore`,
`.gitattributes`, and `examples/make_synthetic_input.py` (a synthetic
positive-control dataset generator used only for the demo).

---

## Files deliberately NOT included

The original folder contained six plotting-script variants. Only the one that
produced the published figure is shipped, to avoid ambiguity about which was
used:

| Original file | Status |
| --- | --- |
| `plot_clonal_coupling_components_curved_with_legends.py` | **shipped** as `scripts/03_...` |
| `plot_clonal_coupling_components.py` | not shipped — no `--min-shared-clones`, straight edges |
| `plot_clonal_coupling_components_clustered.py` | not shipped — straight edges, single legend |
| `plot_clonal_coupling_components_curved.py` | not shipped — curved, but no edge-width legend or `--node-labels-only` |
| `plot_clonal_coupling_components_spring_repulsion.py` | not shipped — superset adding a `spring_repulsion` layout mode |
| `plot_clonal_coupling_components_spring_repulsion_backup.py` | not shipped — as above; also calls `np.cross` on 2-D vectors, which errors on numpy ≥ 2.0 |

The `spring_repulsion` variant is otherwise identical to the shipped script and
was used to produce the alternative `plotting_results_spring_repulsion/` figures,
which are not the published panel. Say so if you want it added back as an
optional extra.

`env_activate.txt` is not shipped: its recorded commands do not correspond to the
published figure. See KNOWN_ISSUES #6, and README section 5 for the correct
command.
