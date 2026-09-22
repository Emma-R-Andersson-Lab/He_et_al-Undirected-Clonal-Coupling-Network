# Undirected Clonal Coupling Network

Replicate-aware, abundance-corrected quantification of clonal coupling between
transcriptomic states, and construction of the statistically filtered undirected
clonal coupling network reported in He *et al.*

Lineage barcodes recovered from single cells tell you which transcriptomic
states share a common progenitor. This pipeline turns that into a network: nodes
are transcriptomic clusters, and an edge means two clusters are populated by the
same clones more often than a within-replicate permutation null allows.

The method builds on the lineage-coupling framework of Bandler *et al.*
([Nature 601, 404-409, 2022](https://doi.org/10.1038/s41586-021-04237-0);
[code](https://github.com/mayer-lab/Bandler-et-al_lineage)), itself based on the
approach of Wagner *et al.*, extended here to account explicitly for biological
replicates and to emit an undirected, FDR-filtered network.

Full method: [docs/methods.md](docs/methods.md).
Verified defects and gotchas: [LIMITATIONS.md](LIMITATIONS.md).

```
scripts/
  01_export_lineage_input.R              Seurat  ->  cell-level CSV tables
  02_clonal_coupling_network.py          coupling scores, permutation null, FDR
  03_plot_clonal_coupling_components.py  publication figures, one per component
examples/
  demo_data/                             small synthetic dataset for the demo
  make_synthetic_input.py                regenerates demo_data from scratch
docs/
  methods.md                             formal method description
```

**Data availability.** The real input data behind the published analysis
(102,077 barcode-positive cells, 8,315 embryo-specific clones, 45 clusters, 9
embryos) is not distributed in this repository. `[Add accession / repository
link / contact statement here once available.]` The pipeline itself is fully
exercised end to end using the synthetic dataset in `examples/demo_data/` (see
Section 3).

---

## 1. System requirements

**Operating systems.** Any OS supporting Python 3.10+ and R 4.x. No
platform-specific code is used.

- Developed and run for the paper on **Linux** (x86-64).
- Independently tested on **Windows 11** (10.0.26200).

**Software dependencies.** Exact versions used to produce the published results,
and the versions the pipeline has additionally been tested against:

| Package | Used for publication | Also tested with |
| --- | --- | --- |
| Python | 3.12.3 | 3.12.10 |
| numpy | 2.5.1 | 2.4.6 |
| pandas | 3.0.5 | 2.3.3 |
| networkx | 3.6.1 | 3.6.1 |
| matplotlib | 3.11.1 | 3.11.1 |
| joblib | 1.5.3 | 1.5.3 |
| plotly *(optional)* | 6.9.0 | 6.9.0 |
| pillow | 12.3.0 | 12.3.0 |

| R package | Used for publication | Also tested with |
| --- | --- | --- |
| R | 4.6 | 4.6.1 |
| Seurat | 5.x | 5.5.1 |
| dplyr | — | 1.2.1 |
| tibble | — | 3.3.1 |
| tidyr | — | 1.3.2 |
| readr *(optional)* | — | falls back to base R |

`plotly` is needed only for the interactive HTML network; everything else runs
without it.

**Hardware.** No non-standard hardware. A normal desktop or laptop is
sufficient. The permutation step is embarrassingly parallel and scales close to
linearly with core count, so more cores directly reduce runtime. Peak memory for
the published dataset (102,077 cells, 8,315 clones, 45 clusters) stayed under
roughly 2 GB.

---

## 2. Installation guide

```bash
git clone https://github.com/Emma-R-Andersson-Lab/He_et_al-Undirected-Clonal-Coupling-Network.git
cd He_et_al-Undirected-Clonal-Coupling-Network
python -m pip install -r requirements.txt
```

or with conda:

```bash
conda env create -f environment.yml
conda activate clonal-coupling
```

For the R export step:

```r
install.packages(c("Seurat", "dplyr", "tibble", "tidyr"))
```

**Typical install time on a normal desktop computer:** about 2-3 minutes for the
Python dependencies over a normal broadband connection (all are pre-built
wheels; no compilation). Installing Seurat from source takes considerably
longer, typically 10-20 minutes, and is only needed for step 1.

No build or compilation step is required for this repository itself; the scripts
run directly.

---

## 3. Demo

The demo generates a synthetic dataset with **known** structure, so you can
confirm the installation reproduces a correct answer. Clones are drawn from
three overlapping "fate modules" (clusters 0-2, 3-5 and 6-8), while clusters
9-11 receive clones independently of module.

A ready-made copy of this dataset is committed at `examples/demo_data/` (3,590
barcode-positive cells, 600 clones, 12 clusters, 4 embryos; 141 KB), so the demo
runs without generating anything first. To regenerate it from scratch, run
`python examples/make_synthetic_input.py --outdir examples/demo_data`.

**Instructions**

```bash
python scripts/02_clonal_coupling_network.py \
  --input examples/demo_data/lineage_cells.csv \
  --cluster-summary examples/demo_data/cluster_summary.csv \
  --outdir example_output \
  --permutations 2000

python scripts/03_plot_clonal_coupling_components.py \
  --analysis-dir example_output \
  --outdir example_figures \
  --layout-mode spring \
  --node-labels-only
```

**Expected output**

Step 2 prints:

```
Significant positive edges: 9
Backbone edges:             6
Isolated nodes:             3
```

The nine edges are exactly the within-module pairs `0-1, 0-2, 1-2, 3-4, 3-5,
4-5, 6-7, 6-8, 7-8`, with **zero** false positives, and clusters 9, 10 and 11
are left unconnected. Step 3 then writes three PNG/SVG pairs, one per module,
plus `connected_component_manifest.csv` and `isolated_nodes.csv`.

Anything other than 9 edges / 3 isolates indicates an installation problem.

**Expected run time on a normal desktop computer:** about **20 seconds** total
(measured on 16 cores: 14 s for 2,000 permutations, 5 s to plot). On 4 cores
expect roughly 50-60 seconds.

---

## 4. Instructions for use

### Step 1 — export from Seurat

Reads a Seurat object carrying, per cell, a transcriptomic cluster label, a
lineage barcode, and a replicate (embryo) ID.

Interactively, with the object already in your session:

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

Or headless against a saved object:

```bash
Rscript scripts/01_export_lineage_input.R \
  --rds seurat_object.rds \
  --cluster-col res0_6_valve_subset_c14_refined \
  --clone-col cloneid_prefixed \
  --embryo-col embryo \
  --outdir input
```

Pass `--cluster-col IDENT` to use the object's active identities.

| Output | Purpose |
| --- | --- |
| `lineage_cells.csv` | One row per barcode-positive cell. The only file required downstream. |
| `cluster_summary.csv` | Node annotation: `n_total_cells` over all clustered cells, and `cluster_order` preserving your factor level ordering. |
| `clone_summary.csv` | Per-clone size and cluster occupancy. |
| `clone_cluster_counts.csv` | Clone-by-cluster count matrix, long form. |
| `embryo_summary.csv` | Per-replicate totals. |

Two conventions are enforced here and re-validated in step 2: the same barcode
in two different embryos is **two independent clones** (encoded as
`clone_uid = "<embryo>::<clone_id>"`), and node size derives from *all* clustered
cells while coupling statistics use only barcode-positive cells.

### Step 2 — coupling scores and permutation null

```bash
python scripts/02_clonal_coupling_network.py \
  --input input/lineage_cells.csv \
  --cluster-summary input/cluster_summary.csv \
  --outdir results \
  --permutations 10000 \
  --n-jobs -1
```

| Flag | Default | Meaning |
| --- | --- | --- |
| `--permutations` | 10000 | Within-embryo label permutations. |
| `--min-shared-clones` | 10 | Minimum clones contributing to an edge. |
| `--fdr` | 0.05 | Maximum BH-adjusted enrichment *q*-value. |
| `--min-z` | 0 | Minimum coupling *z*-score. Do not set below 0 (see LIMITATIONS #3). |
| `--n-jobs` | -1 | Parallel workers; -1 uses all cores. |
| `--seed` | 1234 | Makes the run reproducible. |

An edge enters the network only when $z>0$, shared clones $\ge 10$ and
$q \le 0.05$.

**Run time on real data.** The published dataset (102,077 cells, 8,315 clones,
45 clusters, 9 embryos) at 10,000 permutations takes **about 6.5 minutes**,
measured end to end on a 16-core Windows laptop.

The permutation kernel costs 82 ms per permutation per core, so the theoretical
floor on 16 cores is ~1 minute; the rest is worker startup and per-task data
transfer, which is worse on Windows (spawn) than on Linux (fork). If the run
feels slow, raise `--batch-size`: the default of 50 creates 200 tasks, whereas
`--batch-size 500` creates 20 and amortises that overhead away. Single-core is
roughly 14 minutes.

Outputs: `pairwise_coupling_all.csv` (every pair, significant or not),
`pairwise_coupling_significant_positive_edges.csv` (the network edges),
`clonal_coupling_nodes.csv`, the full score / null / *z* / shared-clone / FDR
matrices, GraphML for Cytoscape or Gephi, a *z*-score heatmap, and an
interactive HTML network.

> **Sizing the permutation count.** The smallest attainable *p*-value is
> 1/(*B*+1). After BH across *K*(*K*-1)/2 pairs this can stop a lone strong edge
> from reaching *q* ≤ 0.05 at all. Use `--permutations` ≥ *K*(*K*-1)/2 ÷ 0.05.
> See [LIMITATIONS.md](LIMITATIONS.md) #2 — including why this did **not**
> affect the published result.

### Step 3 — publication figures

Redraws the network one connected component per figure. Never recomputes
statistics; only ever restricts the edge set further.

```bash
python scripts/03_plot_clonal_coupling_components.py \
  --analysis-dir results \
  --outdir figures \
  --layout-mode spring \
  --node-labels-only
```

| Flag | Default | Meaning |
| --- | --- | --- |
| `--layout-mode` | `community` | `community` detects coupling modules and separates them into distinct blobs; `spring` is one global force-directed layout. |
| `--edge-routing` | `curved` | Deterministic arcs with white haloes at crossings; `straight` for plain segments. |
| `--node-size-column` | `n_total_cells` | Node-table column driving node area. |
| `--node-labels-only` | off | Strip title, subtitle, footer, per-node counts and edge statistics, keeping only cluster names and the legends. |
| `--min-shared-clones` | 0 | Additional plotting-time filter. |
| `--edge-width-factor` | 0.35 | Edge width = factor x *z*. |
| `--colors-csv` | none | CSV with `cluster`,`color` to override node colours without editing the script. |
| `--transparent` | off | Transparent PNG/SVG backgrounds. |

Node colours come from the dictionaries at the top of the script
(`C14_CLUSTER_COLORS` for `c14_`-prefixed labels, `INTEGRATED_CLUSTER_COLORS`
otherwise). Unmapped clusters are drawn grey and listed in
`unmatched_cluster_colors.csv`.

**How to read the figure.** Node = transcriptomic cluster; node area is
proportional to `--node-size-column`. Edge width is directly proportional to the
coupling *z*-score, so *z* = 10 is exactly twice as thick as *z* = 5. Dark edges
are the **maximum-spanning-forest backbone**: the strongest subset of edges
needed to connect each component without loops. The backbone is a readability
aid, not a measure of coupling strength and *not* a lineage trajectory, though it
is suggestive of one. Grey edges are the remaining significant couplings.

> **Legend caveat.** The node-size legend is hard-coded to read `"<n> cells"`
> regardless of which column `--node-size-column` points at. If you size nodes by
> anything other than a cell count, the legend text is wrong and must be
> corrected during figure assembly. This affects the published figure — see
> LIMITATIONS #7.

SVGs are written with editable text (`svg.fonttype = "none"`) for downstream
figure assembly.

---

## 5. Reproducing the published results

The real input data is **not included in this repository**
(see Data availability, above). Once obtained, it takes the same two files
required by any run of the pipeline:

| File | Contents |
| --- | --- |
| `lineage_cells.csv` | 102,077 barcode-positive cells — the analysis input |
| `cluster_summary.csv` | 45 clusters with cell counts, clone counts and ordering |

The full published dataset comprises **102,077 barcode-positive cells, 8,315
embryo-specific clones, 45 transcriptomic clusters and 9 embryos**, giving 990
unique cluster pairs; `clone_summary.csv`, `clone_cluster_counts.csv` and
`embryo_summary.csv` are supporting tables, not required to run the pipeline.

```bash
# Step 2 - analysis
python scripts/02_clonal_coupling_network.py \
  --input lineage_cells.csv \
  --cluster-summary cluster_summary.csv \
  --outdir results_c14_refined \
  --permutations 10000 \
  --seed 1234

# Step 3 - the network figure in the paper
python scripts/03_plot_clonal_coupling_components.py \
  --analysis-dir results_c14_refined \
  --outdir plotting_results \
  --layout-mode spring \
  --node-size-column n_unique_clones \
  --min-shared-clones 10 \
  --node-labels-only
```

This yields 158 edges across 4 connected components plus 3 isolated clusters.
The main panel is component 1: **30 nodes, 132 edges**, *z* ranging 2.87 to
50.17, written as `component_01_nodes30_edges132.png/.svg`.

Verified properties of that run (reproduced internally against the real
dataset prior to its removal from this repository):

| Quantity | Value |
| --- | --- |
| Permutations | 10,000 (empirical *p* floor 9.999e-5) |
| Significant edges (*z*>0, *q*≤0.05, ≥10 shared clones) | 158 — the plotted network |
| Minimum *z* among plotted edges | 2.868 |
| Maximum *q* among plotted edges | 0.0351 |
| Non-finite *z*-scores | 0 |
| Zero-variance nulls | 0 |

**This reproduces exactly.** Re-running the command above on the real dataset
returns a *z*-score matrix bit-identical to the published one (maximum absolute
difference 0.0000 across all 45x45 entries) and the same 158 edges. The
permutation is fully deterministic given `--seed`, and it held across a numpy
version change (2.5.1 for the paper, 2.4.6 for the check).

Two notes on exactness:

1. **Clone threshold.** The published step-2 run used the old default of 3 and
   re-filtered at 10 during plotting. This repository defaults to **10** in step
   2 so that `pairwise_coupling_significant_positive_edges.csv` matches the
   criterion stated in the paper directly. The plotted edge set is identical
   either way (confirmed: 158 edges in both routes); only that intermediate CSV
   differs, dropping from 172 to 158 rows.
2. **Layout reproducibility.** Layouts are seeded (`--layout-seed`, default
   1234) and reproduce exactly on a fixed numpy version. Across numpy versions
   the force-directed initialisation can differ marginally: re-rendering the
   published figure under numpy 2.4.6 rather than 2.5.1 reproduced the identical
   topology, node sizes and legends, with the cropped figure width differing by
   2.2%. Node positions are also written to
   `results/network_layout_coordinates.csv` if you need to pin them exactly.

---

## Interpretation and caveats

- **The network is undirected.** An edge says two states share clonal origin. It
  does not say which came first, nor that one gives rise to the other.
- **Coupling is relative to a within-replicate null.** *z* > 0 means more shared
  clonal history than expected once replicate identity and cluster abundance are
  held fixed. It is not an absolute measure.
- **Absence of an edge is weak evidence.** Sparse barcoding, small clusters or
  too few permutations can all suppress a real edge. Check
  `pairwise_coupling_all.csv` before concluding two states are uncoupled.

---

## Citation

The manuscript describing this work is in preparation. Until it appears, please
cite this repository directly — machine-readable metadata is in
[CITATION.cff](CITATION.cff).

Please also cite the framework this method builds on:

> Bandler, R. C., Vitali, I., Delgado, R. N., Ho, M. C., Dvoretskova, E.,
> Ibarra Molinas, J. S., Frazel, P. W., Mohammadkhani, M., Machold, R.,
> Maedler, S., Liddelow, S. A., Nowakowski, T. J., Fishell, G. & Mayer, C.
> Single-cell delineation of lineage and genetic identity in the mouse brain.
> *Nature* **601**, 404-409 (2022). https://doi.org/10.1038/s41586-021-04237-0

## License

MIT — see [LICENSE](LICENSE). MIT is approved by the Open Source Initiative.
