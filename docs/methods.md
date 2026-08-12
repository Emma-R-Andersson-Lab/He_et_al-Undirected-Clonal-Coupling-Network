# Methods — Clonal coupling undirected lineage network

This is the formal description of what `scripts/02_clonal_coupling_network.py`
computes. It corresponds to the methods section of the accompanying paper.

Clonal coupling between transcriptomic states was quantified by building on the
lineage-coupling framework described by Bandler *et al.*, "Single-cell
delineation of lineage and genetic identity in the mouse brain"
([mayer-lab/Bandler-et-al_lineage](https://github.com/mayer-lab/Bandler-et-al_lineage)),
itself based on the approach of Wagner and colleagues. In that framework the
contribution of each clone to a pair of cell states is normalised by the total
size of that clone, and the observed coupling is standardised against a
permutation-derived null distribution. We adapted it to account explicitly for
biological replicates and to construct an undirected, statistically filtered
clonal-coupling lineage network.

---

## 1. Coupling score

For each embryo-specific clone $c$ and transcriptomic cluster $x$, let

$$C_{cx} = \text{number of cells from clone } c \text{ assigned to cluster } x$$

The total size of clone $c$ is

$$n_c = \sum_x C_{cx}$$

For a pair of clusters $x$ and $y$, clone $c$ counts as **shared** only when at
least one of its cells is present in each cluster:

$$
I_{c,xy} =
\begin{cases}
1 & C_{cx} > 0 \ \text{and} \ C_{cy} > 0 \\
0 & \text{otherwise}
\end{cases}
$$

For a shared clone, the number of clone-derived cells occupying either member of
the pair is

$$k_{c,xy} = C_{cx} + C_{cy}$$

and its contribution to the pairwise coupling is that count as a fraction of the
whole clone:

$$p_{c,xy} = \frac{k_{c,xy}}{n_c} = \frac{C_{cx} + C_{cy}}{n_c}$$

This normalisation is what makes the score **abundance-corrected**: a two-cell
clone split across states $x$ and $y$ contributes exactly as much as a
hundred-cell clone split the same way, so large clones cannot dominate.

The clonal coupling score between clusters $x$ and $y$ is the sum over all
clones:

$$S_{xy} = \sum_c I_{c,xy} \, \frac{C_{cx} + C_{cy}}{n_c}$$

Because every term is non-negative, $S_{xy} \ge 0$ always. Only the standardised
score $Z_{xy}$ defined below can be negative.

---

## 2. Replicate-aware permutation null

To assess whether observed coupling exceeds chance, we build an embryo-aware
permutation null following Bandler *et al.*

Cluster labels are permuted independently **within each biological replicate**.
The permutation holds the following fixed:

- the clone assignment of every cell,
- all clone sizes $n_c$,
- the number of barcode-positive cells in each cluster, within each replicate,
- the total number of barcode-positive cells in each replicate.

Only the pairing between cells and cluster labels is randomised. Permuting
within rather than across replicates is what prevents between-embryo differences
in cluster abundance from being read as coupling.

For each cluster pair, the coupling score is recomputed across $R = 10{,}000$
independent permutations:

$$S_{xy}^{(1)}, S_{xy}^{(2)}, \ldots, S_{xy}^{(R)}$$

The null mean is

$$\mu_{xy} = \frac{1}{R} \sum_{r=1}^{R} S_{xy}^{(r)}$$

and the null standard deviation is

$$\sigma_{xy} = \sqrt{\frac{1}{R-1} \sum_{r=1}^{R} \left(S_{xy}^{(r)} - \mu_{xy}\right)^2}$$

---

## 3. Standardised coupling

$$Z_{xy} = \frac{S_{xy} - \mu_{xy}}{\sigma_{xy}}$$

$Z_{xy} > 0$ means the clonal coupling between clusters $x$ and $y$ exceeds the
mean expected under the replicate-specific null; $Z_{xy} < 0$ means it falls
below it; values near zero indicate coupling indistinguishable from the null.

$Z_{xy}$ is the quantity encoded as edge width in the network figures.

---

## 4. Empirical significance

Alongside the $Z$ framework we test whether coupling is significantly greater
than the null allows. The one-sided empirical enrichment $p$-value is

$$P_{xy} = \frac{1 + \sum_{r=1}^{R} \mathbb{1}\!\left(S_{xy}^{(r)} \ge S_{xy}\right)}{R + 1}$$

The $+1$ in numerator and denominator prevents an empirical $p$-value of exactly
zero, and bounds the smallest attainable value at $1/(R+1)$.

---

## 5. Multiple-testing correction

$P$-values for all unique cluster pairs are corrected with the
Benjamini-Hochberg false discovery rate. Order them ascending, with $m$ the
number of pairs:

$$P_{(1)} \le P_{(2)} \le \ldots \le P_{(m)}$$

For rank index $i$, the initial adjusted value is

$$A_{(i)} = P_{(i)} \, \frac{m}{i}$$

The final $q$-value follows from the monotonic (step-up) correction, in which
each $A_{(i)}$ is replaced by the smallest adjusted value at its own rank or any
higher rank:

$$q_{(i)} = \min\!\left(A_{(i)}, A_{(i+1)}, \ldots, A_{(m)}\right)$$

Values are then clipped to $[0, 1]$. The monotonic step is what makes the
$q$-values non-decreasing in $p$; omitting it does not give valid BH values.

Correction is applied across all $K(K-1)/2$ unique cluster pairs, **before** the
shared-clone filter of section 6 is applied. This is deliberate and
conservative: the number of tests is not reduced by the subsequent filter.

---

## 6. Edge inclusion criteria

An edge is drawn between two cluster nodes only when all three hold:

1. $Z_{xy} > 0$
2. number of clones contributing to that edge $\ge 10$
3. $q_{xy} \le 0.05$

Criterion 2 is `--min-shared-clones`; criteria 1 and 3 are `--min-z` and
`--fdr`. Every pair, including those failing these tests, is retained in
`pairwise_coupling_all.csv`.

---

## 7. Network rendering

Nodes are transcriptomic clusters, with area proportional to the total number of
cells in the cluster (barcoded and unbarcoded alike).

Edge width is directly proportional to $Z_{xy}$, scaled by a constant factor of
0.35. The factor exists purely because unscaled widths were hard to read; it
preserves exact proportionality, so an edge with $Z = 10$ is drawn exactly twice
as thick as one with $Z = 5$.

Dark edges form the **maximum-spanning forest** of the network: the
highest-weight subset of edges that connects each component without loops,
computed per connected component. The backbone is a readability device that
highlights the strongest connective skeleton. It does not encode coupling
strength beyond the widths already shown, and it is **not** a lineage trajectory,
though it is informative about potential lineage relationships. Grey edges are
the remaining significant positive couplings.

---

## 8. Reporting checklist

When reporting results from this pipeline, state:

- the number of barcode-positive cells, embryo-specific clones, clusters and
  biological replicates;
- the number of permutations, and the random seed;
- the three edge-inclusion thresholds;
- that the same barcode in different replicates was treated as distinct clones;
- that FDR correction was applied across all unique cluster pairs.

The analysis script prints the first two groups at the top of every run.
