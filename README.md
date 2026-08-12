# Undirected Clonal Coupling Network

Replicate-aware, abundance-corrected quantification of clonal coupling between
transcriptomic states, and construction of the statistically filtered undirected
lineage network reported in He *et al.*

Lineage barcodes recovered from single cells tell you which transcriptomic
states share a common progenitor. This pipeline turns that into a network: nodes
are transcriptomic clusters, and an edge means two clusters are populated by the
same clones more often than a within-replicate permutation null allows.

The method builds on the lineage-coupling framework of
[Bandler *et al.*, *Nature* 2021](https://github.com/mayer-lab/Bandler-et-al_lineage)
(itself based on Wagner *et al.*), extended here to account explicitly for
biological replicates and to emit an undirected, FDR-filtered network.

---

## Contents

```
scripts/
  01_export_lineage_input.R              Seurat  ->  cell-level CSV tables
  02_clonal_coupling_network.py          coupling scores, permutation null, FDR
  03_plot_clonal_coupling_components.py  publication figures, one per component
examples/
  make_synthetic_input.py                positive-control dataset generator
docs/
  methods.md                             formal method description
KNOWN_ISSUES.md                          documented defects and gotchas
```

---

## Installation

Python 3.10 or newer.

```bash
python -m pip install -r requirements.txt
```

or, with conda:

```bash
conda env create -f environment.yml
conda activate clonal-coupling
```

`plotly` is optional and only used for the interactive HTML output; the rest of
the pipeline runs without it.

The R step needs `Seurat`, `dplyr`, `tibble` and `tidyr`. `readr` is used for
writing CSVs when present, and the script falls back to base R otherwise.

```r
install.packages(c("Seurat", "dplyr", "tibble", "tidyr"))
```

Verified with Python 3.12.10 (numpy 2.4.6, pandas 2.3.3, networkx 3.6.1,
matplotlib 3.11.1, joblib 1.5.3) and R 4.6.1 (Seurat 5.5.1, dplyr 1.2.1,
tibble 3.3.1, tidyr 1.3.2).

---

## Quick start

Confirm the installation on synthetic data with known structure. Clones are
drawn from three overlapping "fate modules", so the correct answer is known in
advance:

```bash
python examples/make_synthetic_input.py --outdir example_input

python scripts/02_clonal_coupling_network.py \
  --input example_input/lineage_cells.csv \
  --cluster-summary example_input/cluster_summary.csv \
  --outdir example_output \
  --permutations 2000

python scripts/03_plot_clonal_coupling_components.py \
  --analysis-dir example_output \
  --outdir example_figures \
  --layout-mode community \
  --node-labels-only
```

This recovers exactly the nine within-module cluster pairs, leaves the three
unplanted clusters isolated, and produces three 3-node components. Anything else
means the installation is wrong.

---

## The pipeline

### Step 1 — export from Seurat

`scripts/01_export_lineage_input.R` reads a Seurat object carrying, per cell, a
transcriptomic cluster label, a lineage barcode, and a replicate (embryo) ID.

Interactively, with the object already in the session:

```r
source("scripts/01_export_lineage_input.R")

export_lineage_input(
  obj         = my_seurat_object,
  cluster_col = "res0_6_valve_subset_c14_refined",
  clone_col   = "cloneid_prefixed",
  embryo_col  = "embryo",
  outdir      = "input"
)
```

Or from the command line against a saved object:

```bash
Rscript scripts/01_export_lineage_input.R \
  --rds seurat_object.rds \
  --cluster-col res0_6_valve_subset_c14_refined \
  --clone-col cloneid_prefixed \
  --embryo-col embryo \
  --outdir input
```

Pass `--cluster-col IDENT` to use the object's active identities instead of a
metadata column.

**Outputs** (in `--outdir`):

| File | Purpose |
| --- | --- |
| `lineage_cells.csv` | One row per barcode-positive cell. The only file required downstream. |
| `cluster_summary.csv` | Node annotation. Supplies `n_total_cells` (all clustered cells, barcoded or not), which sizes the network nodes, and `cluster_order`, which preserves your factor level ordering. |
| `clone_summary.csv` | Per-clone size and cluster occupancy. |
| `clone_cluster_counts.csv` | The clone-by-cluster count matrix in long form. |
| `embryo_summary.csv` | Per-replicate totals. |

Two conventions are enforced here and re-validated in step 2:

- The same barcode seen in two different embryos is **two independent clones**.
  This is encoded as `clone_uid = "<embryo>::<clone_id>"`.
- Node size comes from *all* clustered cells, while coupling statistics use only
  barcode-positive cells. Keeping these separate stops sparsely-barcoded
  clusters from being drawn as small.

### Step 2 — coupling scores and permutation null

```bash
python scripts/02_clonal_coupling_network.py \
  --input input/lineage_cells.csv \
  --cluster-summary input/cluster_summary.csv \
  --outdir results \
  --permutations 10000 \
  --n-jobs -1
```

For each embryo-specific clone *c* and cluster pair (*x*, *y*), the clone
contributes only if it occupies both clusters, and its contribution is
normalised by its own total size, so a 2-cell clone split across two states
counts as much as a 100-cell clone split the same way:

$$S_{xy}=\sum_c \mathbb{1}(C_{cx}>0 \land C_{cy}>0)\,\frac{C_{cx}+C_{cy}}{n_c}$$

Cluster labels are then permuted **within each embryo**, holding fixed the clone
assignments, the clone sizes, the per-replicate cluster abundances and the
per-replicate number of barcoded cells. The observed score is standardised
against that null ($z$), and a one-sided empirical enrichment $p$-value is
computed and corrected across all unique cluster pairs with Benjamini-Hochberg.

See [docs/methods.md](docs/methods.md) for the full derivation.

**Key options**

| Flag | Default | Meaning |
| --- | --- | --- |
| `--permutations` | 10000 | Within-embryo label permutations. |
| `--min-shared-clones` | 10 | Minimum clones contributing to an edge. |
| `--fdr` | 0.05 | Maximum BH-adjusted enrichment *q*-value. |
| `--min-z` | 0 | Minimum coupling *z*-score. |
| `--n-jobs` | -1 | Parallel workers; -1 uses all cores. |
| `--seed` | 1234 | Makes the permutation run reproducible. |

An edge enters the network only when all three of $z>0$, shared clones $\ge 10$
and $q \le 0.05$ hold.

**Outputs** (in `--outdir`): `pairwise_coupling_all.csv` (every pair, whether
significant or not), `pairwise_coupling_significant_positive_edges.csv` (the
network edges), `clonal_coupling_nodes.csv`, the full score / null / *z* /
shared-clone / FDR matrices, GraphML files for Cytoscape or Gephi, a *z*-score
heatmap, and an interactive HTML network.

### Step 3 — publication figures

Step 2 already writes an overview figure. Step 3 redraws the network **one
connected component per figure**, with curved edge routing, a shared global
scale across components, and the colour scheme used in the paper.

```bash
python scripts/03_plot_clonal_coupling_components.py \
  --analysis-dir results \
  --outdir figures \
  --layout-mode community \
  --node-labels-only
```

This step never recomputes statistics. It reads
`pairwise_coupling_significant_positive_edges.csv` and `clonal_coupling_nodes.csv`
and only ever *restricts* the edge set further.

**Key options**

| Flag | Default | Meaning |
| --- | --- | --- |
| `--layout-mode` | `community` | `community` detects coupling modules and separates them; `spring` uses one global force-directed layout. |
| `--edge-routing` | `curved` | Deterministic arcs with white haloes at crossings; `straight` for plain segments. |
| `--node-labels-only` | off | Strip title, subtitle, footer, per-node counts and edge statistics, keeping only cluster names and the legends. Used for the paper figures. |
| `--min-shared-clones` | 0 | Additional plotting-time filter. 0 keeps everything step 2 accepted. |
| `--edge-width-factor` | 0.35 | Edge width = factor x *z*. |
| `--transparent` | off | Transparent PNG/SVG backgrounds. |

Node colours come from the dictionaries at the top of the script
(`C14_CLUSTER_COLORS` for `c14_`-prefixed labels, `INTEGRATED_CLUSTER_COLORS`
otherwise). Unmapped clusters are drawn grey and listed in
`unmatched_cluster_colors.csv`. Override per-cluster colours without editing the
script using `--colors-csv`, a CSV with `cluster` and `color` columns.

**How to read the figure.** Node = transcriptomic cluster; node area = total
cells in that cluster. Edge width is directly proportional to the coupling
*z*-score, so an edge with *z* = 10 is exactly twice as thick as *z* = 5. Dark
edges are the **maximum-spanning-forest backbone**: the strongest subset of
edges needed to connect the network without loops. The backbone is a readability
aid, not a measure of coupling strength and *not* a lineage trajectory, though
it is suggestive of one. Grey edges are the remaining significant couplings.

SVGs are written with editable text (`svg.fonttype = "none"`) for downstream
figure assembly.

---

## Reproducing the published figures

The paper's network figures were produced with:

```bash
python scripts/02_clonal_coupling_network.py \
  --input input/lineage_cells.csv \
  --cluster-summary input/cluster_summary.csv \
  --outdir results \
  --permutations 10000 \
  --seed 1234

python scripts/03_plot_clonal_coupling_components.py \
  --analysis-dir results \
  --outdir figures \
  --layout-mode community \
  --edge-routing curved \
  --node-labels-only
```

Two notes on exact reproduction:

1. **Clone threshold.** The published run used `--min-shared-clones 3` in step 2
   and then re-filtered at 10 during plotting. This repository defaults to 10 in
   step 2 instead, so that
   `pairwise_coupling_significant_positive_edges.csv` matches the criterion
   stated in the paper directly. The final edge set is identical either way;
   only the intermediate CSV differs.
2. **Layout.** `--layout-mode community` gives the modular layout. The
   alternative panel in the same figure series used `--layout-mode spring`. Both
   are seeded (`--layout-seed`, default 1234) and reproducible.

---

## Interpretation and caveats

- **The network is undirected.** An edge says two states share clonal origin. It
  does not say which state came first, nor that one gives rise to the other.
- **Coupling is relative to a within-replicate null.** $z > 0$ means more shared
  clonal history than expected once replicate identity and cluster abundance are
  held fixed. It is not an absolute measure.
- **Absence of an edge is weak evidence.** Sparse barcoding, small clusters or
  too few permutations can all suppress a real edge. Check
  `pairwise_coupling_all.csv` before concluding two states are uncoupled.
- **The permutation count bounds achievable significance.** With *B*
  permutations the smallest possible *p*-value is 1/(*B*+1). After BH across
  *K*(*K*-1)/2 pairs this can prevent an isolated strong edge from reaching
  *q* <= 0.05 at all. See [KNOWN_ISSUES.md](KNOWN_ISSUES.md) for the sizing rule;
  with many clusters, raise `--permutations`.

Please read [KNOWN_ISSUES.md](KNOWN_ISSUES.md) before modifying the analysis
script. It documents verified defects that are deliberately left in place so the
published results remain byte-reproducible.

---

## Citation

If you use this code, please cite the accompanying paper and the framework it
builds on. See [CITATION.cff](CITATION.cff).

## License

MIT — see [LICENSE](LICENSE).
