#!/usr/bin/env python3
"""
Abundance-corrected clonal coupling network.

Input:
    lineage_cells.csv
        One row per barcode-positive cell with columns:
        cell_id, embryo, clone_id, clone_uid, cluster

    cluster_summary.csv
        Optional node summary exported by the accompanying R script.

Method:
    1. Construct a clone-by-cluster cell-count matrix C.
    2. For clusters i and j, calculate the observed weighted coupling score:

           S_ij = sum_c I(C_ci > 0, C_cj > 0) * (C_ci + C_cj) / n_c

       where n_c is the total number of cells in clone c.

    3. Shuffle cluster labels among barcode-positive cells within each embryo.
       This preserves exactly:
           - clone sizes,
           - the number of barcode-positive cells in each cluster per embryo,
           - the number of barcode-positive cells in each embryo.

    4. Convert the observed score to a null-standardized z-score and calculate
       one-sided empirical enrichment p-values followed by BH FDR correction.

    5. Plot only positively enriched, FDR-significant edges with at least the
       requested number of independently shared clones.

Negative scores are retained in the complete output matrices but are not drawn
in the positive coupling network.


For each clone (c) and transcriptomic cluster pair x + y

k = C_cx + C_cy

Where C_cx and C_cy are the number of cells from clone c in clusters x and y.

Now we have n which is the total amount of cells that one clone has (so if a cloneid is 2 then all the cells that have cloneid 2)

The contribution this clone has to the overall lineage coupling score becomes a fraction of how many cells are part of one clone pair k :

p = k / n 

Now the total lineage coupling score between 2 transcriptomic clusters become the sum of all of these contributions p. 

s = sum( p )

Then we use the permutation shuffling that they use in the Bandler paper. Imagine we have 10 clusters, but we look at the coupling between x and y (between these we get a real sum s). We then shuffle, within each embryo, the transcriptomic labels of each cell so that they are randomly assigned to any of the other barcoded cells in these 10 clusters while keeping the cluster sizes still and then we calculate the new score between clusters x and y (which should be lower than the sum s if it is a real connection, otherwise negative value if lower than chance or very close to 0 if by chance). We do this 10 000 times so we get s 1, s 2, s 3 .... s 10 000, then we calculate the average score of these 10 000 shuffles u , we also calculate how much the random scores vary from each other o  (standard deviation).  

Then we get our final coupling z score between clusters x and y by:

z = ( s - u ) / o

Then we do it for all cluster pairs in the entire dataset.

Then we plot the node as the transcriptomic cluster, the edge as the final z-score.

What is decided to plot is based on statistical significance using the randomly shuffled 10 000 scores and the real score. We look at how many times the random shuffled score is higher than the actual real score and then calculate the p-value based on that. So imagine we have a score of 5 and in these 10 000 times we get 20 shuffles that were at least as strong or stronger than the real score of 5, then we calculate the p-value by:

p = (20 + 1) / (10 000 + 1) = 0.002 p value

So now we have a p-value for every edge, so if we have 100 edges we have 100 p-values. We then look at where the p-value is in a ranking of smallest to largest. Say edge number 37 is has the 50th smallest p-value, we take that p-value rank (so rank 50) and then we can calculate the q-value (FDR) by:

q = p x (100/50) = FDR

so it is:

q = p x (number of edges / place of p-value in list smallest to largest) = FDR

And then finally we only plot an edge if:

Clonal coupling score z > 0
shared clones between 2 clusters >= 10
q-value (FDR) =< 0.05

Edge thickness is the z-score (note, thickness is multiplied with 0.35 to make smaller, the lines became difficult to see otherwise, but they are the same ratio scale-wise), the node is the transcriptomic cluster. The dark line is the maximum-spanning-forest backbone, NOT the coupling strength, they are the edges that are REQUIRED to connect the entire network in a loop-free manner. Imagine you remove all connections, but only keep the strongest ones required to make a fully connected network, then you end up with the backbone. So it is not a lineage trajectory per definition, but does inform us of a potential lineage.

"""

from __future__ import annotations

import argparse
import math
import os
import re
import warnings
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from joblib import Parallel, delayed, effective_n_jobs

try:
    import plotly.graph_objects as go
except ImportError:
    go = None


# =============================================================================
# GENERAL UTILITIES
# =============================================================================


def natural_key(value: str) -> list[object]:
    """Natural sort key: c2v_2 comes before c2v_10."""
    return [int(token) if token.isdigit() else token.lower()
            for token in re.split(r"(\d+)", str(value))]


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR adjustment with monotonic correction."""
    p = np.asarray(p_values, dtype=float)
    q = np.full(p.shape, np.nan, dtype=float)

    valid = np.isfinite(p)
    if not np.any(valid):
        return q

    pv = p[valid]
    order = np.argsort(pv)
    ranked = pv[order]
    m = len(ranked)

    adjusted = ranked * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)

    restored = np.empty_like(adjusted)
    restored[order] = adjusted
    q[valid] = restored
    return q


def safe_zscore(
    observed: np.ndarray,
    null_mean: np.ndarray,
    null_sd: np.ndarray,
) -> np.ndarray:
    """Calculate z-scores while handling zero-variance null distributions."""
    z = np.full_like(observed, np.nan, dtype=float)
    regular = null_sd > 0
    z[regular] = (observed[regular] - null_mean[regular]) / null_sd[regular]

    zero_sd = ~regular
    equal = zero_sd & np.isclose(observed, null_mean)
    above = zero_sd & (observed > null_mean)
    below = zero_sd & (observed < null_mean)

    z[equal] = 0.0
    z[above] = np.inf
    z[below] = -np.inf
    return z


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build an abundance-corrected clonal coupling network."
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Path to lineage_cells.csv exported by the R script.",
    )
    parser.add_argument(
        "--cluster-summary",
        type=Path,
        default=None,
        help="Optional path to cluster_summary.csv exported by the R script.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Directory for all result files.",
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=10_000,
        help="Number of within-embryo label permutations. Default: 10000.",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=-1,
        help="Parallel workers. Use -1 for all available cores. Default: -1.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=50,
        help="Permutations per parallel task. Default: 50.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1234,
        help="Random seed. Default: 1234.",
    )
    parser.add_argument(
        "--min-shared-clones",
        type=int,
        default=10,
        help="Minimum observed shared clones required for an edge. Default: 10.",
    )
    parser.add_argument(
        "--fdr",
        type=float,
        default=0.05,
        help="Maximum BH-adjusted enrichment q-value for an edge. Default: 0.05.",
    )
    parser.add_argument(
        "--min-z",
        type=float,
        default=0.0,
        help="Minimum positive coupling z-score for an edge. Default: 0.",
    )
    parser.add_argument(
        "--layout-k",
        type=float,
        default=None,
        help="Optional NetworkX spring-layout k. Default: 1/sqrt(number of nodes).",
    )
    parser.add_argument(
        "--layout-iterations",
        type=int,
        default=1000,
        help="Spring-layout iterations. Default: 1000.",
    )
    return parser.parse_args()


# =============================================================================
# COUPLING SCORE
# =============================================================================


def coupling_score_from_counts(
    counts: np.ndarray,
    inverse_clone_sizes: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return weighted coupling score and shared-clone count matrices.

    counts has shape: number_of_clones x number_of_clusters.
    """
    presence = counts > 0
    weighted_fraction = counts.astype(np.float64) * inverse_clone_sizes[:, None]

    presence_float = presence.astype(np.float64, copy=False)
    score = (
        weighted_fraction.T @ presence_float
        + presence_float.T @ weighted_fraction
    )

    shared_clones = presence.astype(np.int64).T @ presence.astype(np.int64)
    return score, shared_clones


def counts_from_assignments(
    clone_codes: np.ndarray,
    cluster_codes: np.ndarray,
    n_clones: int,
    n_clusters: int,
) -> np.ndarray:
    """Construct clone-by-cluster cell counts."""
    counts = np.zeros((n_clones, n_clusters), dtype=np.int32)
    np.add.at(counts, (clone_codes, cluster_codes), 1)
    return counts


# =============================================================================
# PERMUTATION ENGINE
# =============================================================================


def run_permutation_batch(
    batch_seed: int,
    n_permutations: int,
    clone_codes: np.ndarray,
    cluster_codes: np.ndarray,
    embryo_indices: tuple[np.ndarray, ...],
    n_clones: int,
    n_clusters: int,
    inverse_clone_sizes: np.ndarray,
    observed_score: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """Run one independent permutation batch and return aggregate statistics."""
    rng = np.random.default_rng(batch_seed)

    score_sum = np.zeros((n_clusters, n_clusters), dtype=np.float64)
    score_sum_sq = np.zeros((n_clusters, n_clusters), dtype=np.float64)
    greater_equal_count = np.zeros((n_clusters, n_clusters), dtype=np.int64)
    less_equal_count = np.zeros((n_clusters, n_clusters), dtype=np.int64)

    permuted_clusters = np.empty_like(cluster_codes)

    for _ in range(n_permutations):
        for indices in embryo_indices:
            permuted_clusters[indices] = rng.permutation(cluster_codes[indices])

        permuted_counts = counts_from_assignments(
            clone_codes=clone_codes,
            cluster_codes=permuted_clusters,
            n_clones=n_clones,
            n_clusters=n_clusters,
        )

        permuted_score, _ = coupling_score_from_counts(
            permuted_counts,
            inverse_clone_sizes,
        )

        score_sum += permuted_score
        score_sum_sq += permuted_score * permuted_score
        greater_equal_count += permuted_score >= observed_score
        less_equal_count += permuted_score <= observed_score

    return (
        score_sum,
        score_sum_sq,
        greater_equal_count,
        less_equal_count,
        n_permutations,
    )


def generate_batch_plan(
    n_permutations: int,
    batch_size: int,
    seed: int,
) -> list[tuple[int, int]]:
    """Create reproducible independent seeds and permutation counts per batch."""
    if n_permutations <= 0:
        raise ValueError("--permutations must be greater than zero.")
    if batch_size <= 0:
        raise ValueError("--batch-size must be greater than zero.")

    batch_counts: list[int] = []
    remaining = n_permutations
    while remaining > 0:
        current = min(batch_size, remaining)
        batch_counts.append(current)
        remaining -= current

    seed_sequence = np.random.SeedSequence(seed)
    child_sequences = seed_sequence.spawn(len(batch_counts))
    child_seeds = [
        int(sequence.generate_state(1, dtype=np.uint32)[0])
        for sequence in child_sequences
    ]

    return list(zip(child_seeds, batch_counts))


# =============================================================================
# TABLE CONSTRUCTION
# =============================================================================


def load_inputs(
    input_path: Path,
    cluster_summary_path: Path | None,
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    cells = pd.read_csv(input_path, dtype=str)

    required = {"cell_id", "embryo", "clone_id", "cluster"}
    missing = sorted(required - set(cells.columns))
    if missing:
        raise ValueError(
            "Input is missing required column(s): " + ", ".join(missing)
        )

    if "clone_uid" not in cells.columns:
        cells["clone_uid"] = (
            cells["embryo"].astype(str) + "::" + cells["clone_id"].astype(str)
        )

    keep = ["cell_id", "embryo", "clone_id", "clone_uid", "cluster"]
    cells = cells[keep].copy()

    for column in keep:
        cells[column] = cells[column].astype(str).str.strip()

    invalid = {"", "nan", "NaN", "NA", "None", "none"}
    valid_mask = np.ones(len(cells), dtype=bool)
    for column in keep:
        valid_mask &= ~cells[column].isin(invalid)

    dropped = int((~valid_mask).sum())
    if dropped > 0:
        warnings.warn(f"Dropping {dropped:,} rows with missing/invalid values.")
        cells = cells.loc[valid_mask].copy()

    if cells.empty:
        raise ValueError("No valid rows remain in the input file.")

    if cells["cell_id"].duplicated().any():
        duplicates = cells.loc[cells["cell_id"].duplicated(), "cell_id"].head()
        raise ValueError(
            "Duplicate cell IDs found, including: " + ", ".join(duplicates)
        )

    clone_embryo_counts = cells.groupby("clone_uid")["embryo"].nunique()
    invalid_clone_uids = clone_embryo_counts[clone_embryo_counts > 1]
    if not invalid_clone_uids.empty:
        raise ValueError(
            "Some clone_uid values span multiple embryos. clone_uid must be "
            "embryo-specific. Example: " + str(invalid_clone_uids.index[0])
        )

    cluster_summary = None
    if cluster_summary_path is not None:
        cluster_summary = pd.read_csv(cluster_summary_path)
        if "cluster" not in cluster_summary.columns:
            raise ValueError("cluster_summary.csv must contain a 'cluster' column.")
        cluster_summary["cluster"] = cluster_summary["cluster"].astype(str)

    return cells, cluster_summary


def determine_cluster_order(
    cells: pd.DataFrame,
    cluster_summary: pd.DataFrame | None,
) -> list[str]:
    observed = set(cells["cluster"].astype(str))

    if cluster_summary is not None:
        summary = cluster_summary.copy()
        if "cluster_order" in summary.columns:
            summary["cluster_order"] = pd.to_numeric(
                summary["cluster_order"], errors="coerce"
            )
            summary = summary.sort_values("cluster_order", kind="stable")

        ordered = [
            cluster for cluster in summary["cluster"].astype(str)
            if cluster in observed
        ]
        extras = sorted(observed - set(ordered), key=natural_key)
        return ordered + extras

    return sorted(observed, key=natural_key)


def build_node_table(
    cells: pd.DataFrame,
    cluster_order: list[str],
    counts: np.ndarray,
    cluster_summary: pd.DataFrame | None,
) -> pd.DataFrame:
    barcode_cells = counts.sum(axis=0).astype(int)
    unique_clones = (counts > 0).sum(axis=0).astype(int)

    node_table = pd.DataFrame({
        "cluster": cluster_order,
        "cluster_order": np.arange(1, len(cluster_order) + 1),
        "n_barcoded_cells": barcode_cells,
        "n_unique_clones": unique_clones,
    })

    embryo_counts = (
        cells.groupby("cluster")["embryo"]
        .nunique()
        .reindex(cluster_order, fill_value=0)
        .to_numpy(dtype=int)
    )
    node_table["n_embryos"] = embryo_counts

    if cluster_summary is not None:
        summary = cluster_summary.drop_duplicates("cluster").copy()
        columns_to_merge = [
            column for column in summary.columns
            if column not in {
                "cluster_order",
                "n_barcoded_cells",
                "n_unique_clones",
                "n_embryos",
            }
        ]
        node_table = node_table.merge(
            summary[columns_to_merge],
            on="cluster",
            how="left",
            validate="one_to_one",
        )

    if "n_total_cells" not in node_table.columns:
        node_table["n_total_cells"] = node_table["n_barcoded_cells"]

    node_table["n_total_cells"] = pd.to_numeric(
        node_table["n_total_cells"], errors="coerce"
    ).fillna(node_table["n_barcoded_cells"]).astype(int)

    node_table["barcode_fraction"] = np.where(
        node_table["n_total_cells"] > 0,
        node_table["n_barcoded_cells"] / node_table["n_total_cells"],
        np.nan,
    )

    return node_table.sort_values("cluster_order").reset_index(drop=True)


def make_pairwise_table(
    cluster_order: list[str],
    observed_score: np.ndarray,
    shared_clones: np.ndarray,
    null_mean: np.ndarray,
    null_sd: np.ndarray,
    z_score: np.ndarray,
    p_enrichment: np.ndarray,
    p_depletion: np.ndarray,
    q_enrichment: np.ndarray,
    q_depletion: np.ndarray,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    epsilon = 1e-9

    for i in range(len(cluster_order)):
        for j in range(i + 1, len(cluster_order)):
            observed = float(observed_score[i, j])
            expected = float(null_mean[i, j])

            rows.append({
                "cluster_1": cluster_order[i],
                "cluster_2": cluster_order[j],
                "shared_clones": int(shared_clones[i, j]),
                "observed_score": observed,
                "null_mean": expected,
                "null_sd": float(null_sd[i, j]),
                "z_score": float(z_score[i, j]),
                "log2_enrichment": float(
                    np.log2((observed + epsilon) / (expected + epsilon))
                ),
                "p_enrichment": float(p_enrichment[i, j]),
                "q_enrichment": float(q_enrichment[i, j]),
                "p_depletion": float(p_depletion[i, j]),
                "q_depletion": float(q_depletion[i, j]),
            })

    return pd.DataFrame(rows)


def matrix_to_dataframe(
    matrix: np.ndarray,
    cluster_order: list[str],
) -> pd.DataFrame:
    return pd.DataFrame(matrix, index=cluster_order, columns=cluster_order)


# =============================================================================
# NETWORK AND PLOTTING
# =============================================================================


def build_graphs(
    node_table: pd.DataFrame,
    significant_edges: pd.DataFrame,
) -> tuple[nx.Graph, nx.Graph]:
    graph = nx.Graph()

    for row in node_table.itertuples(index=False):
        attrs = row._asdict()
        cluster = str(attrs.pop("cluster"))
        clean_attrs = {
            key: (value.item() if isinstance(value, np.generic) else value)
            for key, value in attrs.items()
            if pd.notna(value)
        }
        graph.add_node(cluster, **clean_attrs)

    for row in significant_edges.itertuples(index=False):
        graph.add_edge(
            str(row.cluster_1),
            str(row.cluster_2),
            weight=float(row.z_score),
            z_score=float(row.z_score),
            log2_enrichment=float(row.log2_enrichment),
            shared_clones=int(row.shared_clones),
            observed_score=float(row.observed_score),
            null_mean=float(row.null_mean),
            q_enrichment=float(row.q_enrichment),
        )

    backbone = nx.maximum_spanning_tree(graph, weight="weight")
    return graph, backbone


def get_layout(
    graph: nx.Graph,
    seed: int,
    layout_k: float | None,
    iterations: int,
) -> dict[str, np.ndarray]:
    n_nodes = max(graph.number_of_nodes(), 1)
    k = layout_k if layout_k is not None else 1.0 / math.sqrt(n_nodes)

    if graph.number_of_nodes() == 1:
        only_node = next(iter(graph.nodes))
        return {only_node: np.array([0.0, 0.0])}

    return nx.spring_layout(
        graph,
        seed=seed,
        weight="weight",
        k=k,
        iterations=iterations,
        scale=1.0,
    )


def scale_values(
    values: Iterable[float],
    low: float,
    high: float,
) -> list[float]:
    values_array = np.asarray(list(values), dtype=float)
    if values_array.size == 0:
        return []

    finite = np.isfinite(values_array)
    if not np.any(finite):
        return [low] * len(values_array)

    finite_values = values_array[finite]
    minimum = finite_values.min()
    maximum = finite_values.max()

    if np.isclose(minimum, maximum):
        scaled = np.full(values_array.shape, (low + high) / 2.0)
    else:
        scaled = low + (values_array - minimum) * (high - low) / (maximum - minimum)

    scaled[~finite] = high
    return scaled.tolist()


def node_colors_for_clusters(clusters: list[str]) -> dict[str, str]:
    cmap = plt.get_cmap("tab20")
    colors: dict[str, str] = {}
    for index, cluster in enumerate(clusters):
        rgba = cmap(index % cmap.N)
        colors[cluster] = (
            f"#{int(rgba[0] * 255):02x}"
            f"{int(rgba[1] * 255):02x}"
            f"{int(rgba[2] * 255):02x}"
        )
    return colors


EDGE_WIDTH_FACTOR = 0.35


def proportional_edge_width_map(
    graph: nx.Graph,
    factor: float = EDGE_WIDTH_FACTOR,
) -> dict[frozenset[str], float]:
    """
    Convert positive coupling z-scores directly into plotted edge widths.

    width = z_score * factor

    This preserves true proportionality: an edge with z=10 is exactly twice
    as thick as an edge with z=5. Dark backbone edges and grey secondary edges
    use the same width calculation; only colour and opacity differ.

    If an infinite z-score occurs because the permutation null has zero
    variance, it is plotted using the largest finite positive z-score so that
    the figure remains renderable.
    """
    finite_positive_z = [
        float(attrs["z_score"])
        for _, _, attrs in graph.edges(data=True)
        if np.isfinite(float(attrs["z_score"]))
        and float(attrs["z_score"]) > 0
    ]
    largest_finite_z = max(finite_positive_z, default=1.0)

    width_map: dict[frozenset[str], float] = {}

    for u, v, attrs in graph.edges(data=True):
        z = float(attrs["z_score"])

        if not np.isfinite(z):
            warnings.warn(
                f"Non-finite z-score for edge {u!r}–{v!r}; "
                f"using largest finite positive z-score ({largest_finite_z:.3f}) "
                "for plotting only."
            )
            z = largest_finite_z

        if z <= 0:
            raise ValueError(
                f"Graph contains a non-positive plotted edge: {u!r}–{v!r}, z={z}. "
                "Only positive significant edges should be passed to plotting."
            )

        width_map[frozenset((u, v))] = z * factor

    return width_map


def plot_static_network(
    graph: nx.Graph,
    backbone: nx.Graph,
    node_table: pd.DataFrame,
    positions: dict[str, np.ndarray],
    outdir: Path,
) -> None:
    cluster_order = node_table["cluster"].astype(str).tolist()
    colors = node_colors_for_clusters(cluster_order)

    node_size_map = dict(zip(
        node_table["cluster"].astype(str),
        scale_values(node_table["n_total_cells"], 700, 3200),
    ))

    backbone_edges = {frozenset(edge) for edge in backbone.edges()}
    secondary_edges = [
        edge
        for edge in graph.edges()
        if frozenset(edge) not in backbone_edges
    ]

    # One common, directly proportional width scale for every edge.
    # z=10 is exactly twice as thick as z=5.
    edge_width_map = proportional_edge_width_map(graph)

    fig, ax = plt.subplots(figsize=(14, 12))

    if secondary_edges:
        nx.draw_networkx_edges(
            graph,
            positions,
            edgelist=secondary_edges,
            width=[
                edge_width_map[frozenset(edge)]
                for edge in secondary_edges
            ],
            edge_color="#9AA0A6",
            alpha=0.30,
            ax=ax,
        )

    if backbone.number_of_edges() > 0:
        nx.draw_networkx_edges(
            graph,
            positions,
            edgelist=list(backbone.edges()),
            width=[
                edge_width_map[frozenset(edge)]
                for edge in backbone.edges()
            ],
            edge_color="#263238",
            alpha=0.82,
            ax=ax,
        )

    nx.draw_networkx_nodes(
        graph,
        positions,
        nodelist=cluster_order,
        node_size=[node_size_map[node] for node in cluster_order],
        node_color=[colors[node] for node in cluster_order],
        edgecolors="#172033",
        linewidths=1.4,
        alpha=0.98,
        ax=ax,
    )

    nx.draw_networkx_labels(
        graph,
        positions,
        labels={node: node for node in cluster_order},
        font_size=9,
        font_weight="bold",
        ax=ax,
    )

    ax.set_title(
        "Clone-resolved transcriptomic state coupling network\n"
        "Edge width = 0.35 × coupling z-score; "
        "dark edges = maximum-spanning-forest backbone",
        fontsize=15,
        pad=18,
    )
    ax.axis("off")
    fig.tight_layout()

    fig.savefig(
        outdir / "clonal_coupling_network.png",
        dpi=350,
        bbox_inches="tight",
    )
    fig.savefig(
        outdir / "clonal_coupling_network.svg",
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_backbone_only(
    graph: nx.Graph,
    backbone: nx.Graph,
    node_table: pd.DataFrame,
    positions: dict[str, np.ndarray],
    outdir: Path,
) -> None:
    cluster_order = node_table["cluster"].astype(str).tolist()
    colors = node_colors_for_clusters(cluster_order)

    node_size_map = dict(zip(
        node_table["cluster"].astype(str),
        scale_values(node_table["n_total_cells"], 700, 3200),
    ))

    # Use the exact same proportional scale as the full network.
    edge_width_map = proportional_edge_width_map(graph)
    backbone_widths = [
        edge_width_map[frozenset((u, v))]
        for u, v in backbone.edges()
    ]

    fig, ax = plt.subplots(figsize=(14, 12))

    if backbone.number_of_edges() > 0:
        nx.draw_networkx_edges(
            backbone,
            positions,
            edgelist=list(backbone.edges()),
            width=backbone_widths,
            edge_color="#263238",
            alpha=0.85,
            ax=ax,
        )

    nx.draw_networkx_nodes(
        backbone,
        positions,
        nodelist=cluster_order,
        node_size=[node_size_map[node] for node in cluster_order],
        node_color=[colors[node] for node in cluster_order],
        edgecolors="#172033",
        linewidths=1.4,
        alpha=0.98,
        ax=ax,
    )

    nx.draw_networkx_labels(
        backbone,
        positions,
        labels={node: node for node in cluster_order},
        font_size=9,
        font_weight="bold",
        ax=ax,
    )

    ax.set_title(
        "Maximum-spanning-forest clonal coupling backbone\n"
        "Edge width = 0.35 × coupling z-score",
        fontsize=15,
        pad=18,
    )
    ax.axis("off")
    fig.tight_layout()

    fig.savefig(
        outdir / "clonal_coupling_backbone.png",
        dpi=350,
        bbox_inches="tight",
    )
    fig.savefig(
        outdir / "clonal_coupling_backbone.svg",
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_zscore_heatmap(
    z_score: np.ndarray,
    cluster_order: list[str],
    outdir: Path,
) -> None:
    finite = z_score[np.isfinite(z_score)]
    vmax = np.max(np.abs(finite)) if finite.size > 0 else 1.0
    vmax = max(vmax, 1.0)

    fig_size = max(9, 0.34 * len(cluster_order) + 4)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))

    image = ax.imshow(
        z_score,
        cmap="RdBu_r",
        vmin=-vmax,
        vmax=vmax,
        interpolation="nearest",
        aspect="equal",
    )

    ax.set_xticks(np.arange(len(cluster_order)))
    ax.set_yticks(np.arange(len(cluster_order)))
    ax.set_xticklabels(cluster_order, rotation=90, fontsize=8)
    ax.set_yticklabels(cluster_order, fontsize=8)
    ax.set_title("Clonal coupling z-score matrix", pad=14)

    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label("Observed coupling relative to within-embryo null")

    fig.tight_layout()
    fig.savefig(
        outdir / "clonal_coupling_zscore_heatmap.png",
        dpi=350,
        bbox_inches="tight",
    )
    fig.savefig(
        outdir / "clonal_coupling_zscore_heatmap.svg",
        bbox_inches="tight",
    )
    plt.close(fig)


def write_interactive_network(
    graph: nx.Graph,
    backbone: nx.Graph,
    node_table: pd.DataFrame,
    positions: dict[str, np.ndarray],
    outdir: Path,
) -> None:
    if go is None:
        warnings.warn(
            "plotly is not installed; skipping interactive HTML output. "
            "Install it with: python -m pip install plotly"
        )
        return

    cluster_order = node_table["cluster"].astype(str).tolist()
    colors = node_colors_for_clusters(cluster_order)

    node_size_map = dict(zip(
        node_table["cluster"].astype(str),
        scale_values(node_table["n_total_cells"], 18, 50),
    ))

    backbone_edges = {frozenset(edge) for edge in backbone.edges()}

    # Same directly proportional width scale as the static outputs.
    edge_width_map = proportional_edge_width_map(graph)

    figure = go.Figure()

    # Draw each edge separately so width can encode its exact relative z-score.
    for u, v, attrs in graph.edges(data=True):
        x0, y0 = positions[u]
        x1, y1 = positions[v]

        edge_key = frozenset((u, v))
        is_backbone = edge_key in backbone_edges

        figure.add_trace(go.Scatter(
            x=[x0, x1],
            y=[y0, y1],
            mode="lines",
            line={
                # No grey-edge multiplier: both classes use the same width.
                "width": edge_width_map[edge_key],
                "color": (
                    "rgba(38,50,56,0.82)"
                    if is_backbone
                    else "rgba(154,160,166,0.30)"
                ),
            },
            hoverinfo="skip",
            showlegend=False,
        ))

        midpoint_x = (x0 + x1) / 2.0
        midpoint_y = (y0 + y1) / 2.0
        hover = (
            f"{u} ↔ {v}<br>"
            f"z-score: {attrs['z_score']:.3f}<br>"
            f"plotted width: {edge_width_map[edge_key]:.3f}<br>"
            f"log2 enrichment: {attrs['log2_enrichment']:.3f}<br>"
            f"shared clones: {attrs['shared_clones']}<br>"
            f"FDR: {attrs['q_enrichment']:.4g}<br>"
            f"backbone: {'yes' if is_backbone else 'no'}"
        )

        # Transparent midpoint marker provides the edge tooltip.
        figure.add_trace(go.Scatter(
            x=[midpoint_x],
            y=[midpoint_y],
            mode="markers",
            marker={"size": 12, "color": "rgba(0,0,0,0)"},
            hovertemplate=hover + "<extra></extra>",
            showlegend=False,
        ))

    node_rows = node_table.set_index("cluster")
    node_x: list[float] = []
    node_y: list[float] = []
    node_text: list[str] = []
    node_hover: list[str] = []
    node_sizes: list[float] = []
    node_colors: list[str] = []

    for node in cluster_order:
        x, y = positions[node]
        row = node_rows.loc[node]

        node_x.append(float(x))
        node_y.append(float(y))
        node_text.append(node)
        node_sizes.append(node_size_map[node])
        node_colors.append(colors[node])
        node_hover.append(
            f"cluster: {node}<br>"
            f"total cells: {int(row['n_total_cells']):,}<br>"
            f"barcoded cells: {int(row['n_barcoded_cells']):,}<br>"
            f"unique clones: {int(row['n_unique_clones']):,}<br>"
            f"embryos: {int(row['n_embryos']):,}<br>"
            f"network degree: {graph.degree(node)}"
        )

    figure.add_trace(go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers+text",
        text=node_text,
        textposition="top center",
        hovertext=node_hover,
        hoverinfo="text",
        marker={
            "size": node_sizes,
            "color": node_colors,
            "line": {"width": 1.5, "color": "#172033"},
            "opacity": 0.98,
        },
        showlegend=False,
    ))

    figure.update_layout(
        title={
            "text": (
                "Clone-resolved transcriptomic state coupling network"
                "<br><sup>"
                "Edge width = 0.35 × coupling z-score; "
                "dark edges = maximum-spanning-forest backbone; "
                "grey edges = additional significant positive coupling"
                "</sup>"
            ),
            "x": 0.5,
        },
        template="plotly_white",
        hovermode="closest",
        margin={"l": 20, "r": 20, "t": 90, "b": 20},
        xaxis={"visible": False},
        yaxis={"visible": False, "scaleanchor": "x", "scaleratio": 1},
        width=1200,
        height=1000,
    )

    figure.write_html(
        outdir / "clonal_coupling_network_interactive.html",
        include_plotlyjs=True,
        full_html=True,
    )


# =============================================================================
# MAIN
# =============================================================================


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cells, cluster_summary = load_inputs(
        input_path=args.input,
        cluster_summary_path=args.cluster_summary,
    )

    cluster_order = determine_cluster_order(cells, cluster_summary)
    cluster_to_code = {cluster: index for index, cluster in enumerate(cluster_order)}

    clone_categories = pd.Index(cells["clone_uid"].drop_duplicates())
    clone_to_code = {clone: index for index, clone in enumerate(clone_categories)}

    clone_codes = cells["clone_uid"].map(clone_to_code).to_numpy(dtype=np.int32)
    cluster_codes = cells["cluster"].map(cluster_to_code).to_numpy(dtype=np.int32)

    embryo_indices = tuple(
        group.index.to_numpy(dtype=np.int64)
        for _, group in cells.reset_index(drop=True).groupby("embryo", sort=False)
    )

    n_clones = len(clone_categories)
    n_clusters = len(cluster_order)

    observed_counts = counts_from_assignments(
        clone_codes=clone_codes,
        cluster_codes=cluster_codes,
        n_clones=n_clones,
        n_clusters=n_clusters,
    )

    clone_sizes = observed_counts.sum(axis=1).astype(np.float64)
    if np.any(clone_sizes <= 0):
        raise RuntimeError("Encountered an empty clone after matrix construction.")
    inverse_clone_sizes = 1.0 / clone_sizes

    observed_score, shared_clones = coupling_score_from_counts(
        observed_counts,
        inverse_clone_sizes,
    )

    print("\nInput summary")
    print("-------------")
    print(f"Barcode-positive cells: {len(cells):,}")
    print(f"Embryos:               {cells['embryo'].nunique():,}")
    print(f"Embryo-specific clones:{n_clones:>10,}")
    print(f"Clusters:              {n_clusters:>10,}")
    print(f"Permutations:          {args.permutations:>10,}")

    batch_plan = generate_batch_plan(
        n_permutations=args.permutations,
        batch_size=args.batch_size,
        seed=args.seed,
    )

    actual_jobs = effective_n_jobs(args.n_jobs)
    print(f"Parallel workers:      {actual_jobs:>10,}\n")

    batch_results = Parallel(
        n_jobs=args.n_jobs,
        backend="loky",
        verbose=10,
    )(
        delayed(run_permutation_batch)(
            batch_seed=batch_seed,
            n_permutations=batch_count,
            clone_codes=clone_codes,
            cluster_codes=cluster_codes,
            embryo_indices=embryo_indices,
            n_clones=n_clones,
            n_clusters=n_clusters,
            inverse_clone_sizes=inverse_clone_sizes,
            observed_score=observed_score,
        )
        for batch_seed, batch_count in batch_plan
    )

    score_sum = np.zeros_like(observed_score, dtype=np.float64)
    score_sum_sq = np.zeros_like(observed_score, dtype=np.float64)
    greater_equal_count = np.zeros_like(observed_score, dtype=np.int64)
    less_equal_count = np.zeros_like(observed_score, dtype=np.int64)
    completed_permutations = 0

    for result in batch_results:
        batch_sum, batch_sum_sq, batch_ge, batch_le, batch_n = result
        score_sum += batch_sum
        score_sum_sq += batch_sum_sq
        greater_equal_count += batch_ge
        less_equal_count += batch_le
        completed_permutations += batch_n

    if completed_permutations != args.permutations:
        raise RuntimeError(
            f"Expected {args.permutations} permutations but completed "
            f"{completed_permutations}."
        )

    null_mean = score_sum / completed_permutations

    if completed_permutations > 1:
        variance_numerator = (
            score_sum_sq
            - completed_permutations * null_mean * null_mean
        )
        variance_numerator = np.maximum(variance_numerator, 0.0)
        null_variance = variance_numerator / (completed_permutations - 1)
        null_sd = np.sqrt(null_variance)
    else:
        null_sd = np.full_like(null_mean, np.nan)

    z_score = safe_zscore(observed_score, null_mean, null_sd)

    p_enrichment = (
        greater_equal_count.astype(np.float64) + 1.0
    ) / (completed_permutations + 1.0)
    p_depletion = (
        less_equal_count.astype(np.float64) + 1.0
    ) / (completed_permutations + 1.0)

    upper = np.triu_indices(n_clusters, k=1)

    q_enrichment = np.ones_like(p_enrichment)
    q_depletion = np.ones_like(p_depletion)

    q_enrichment_values = bh_fdr(p_enrichment[upper])
    q_depletion_values = bh_fdr(p_depletion[upper])

    q_enrichment[upper] = q_enrichment_values
    q_enrichment[(upper[1], upper[0])] = q_enrichment_values
    q_depletion[upper] = q_depletion_values
    q_depletion[(upper[1], upper[0])] = q_depletion_values

    np.fill_diagonal(z_score, 0.0)
    np.fill_diagonal(p_enrichment, 1.0)
    np.fill_diagonal(p_depletion, 1.0)
    np.fill_diagonal(q_enrichment, 1.0)
    np.fill_diagonal(q_depletion, 1.0)

    pairwise = make_pairwise_table(
        cluster_order=cluster_order,
        observed_score=observed_score,
        shared_clones=shared_clones,
        null_mean=null_mean,
        null_sd=null_sd,
        z_score=z_score,
        p_enrichment=p_enrichment,
        p_depletion=p_depletion,
        q_enrichment=q_enrichment,
        q_depletion=q_depletion,
    )

    significant_edges = pairwise.loc[
        (pairwise["z_score"] > args.min_z)
        & (pairwise["q_enrichment"] <= args.fdr)
        & (pairwise["shared_clones"] >= args.min_shared_clones)
    ].copy()

    significant_edges = significant_edges.sort_values(
        ["z_score", "shared_clones"],
        ascending=[False, False],
    ).reset_index(drop=True)

    node_table = build_node_table(
        cells=cells,
        cluster_order=cluster_order,
        counts=observed_counts,
        cluster_summary=cluster_summary,
    )

    graph, backbone = build_graphs(node_table, significant_edges)
    positions = get_layout(
        graph=graph,
        seed=args.seed,
        layout_k=args.layout_k,
        iterations=args.layout_iterations,
    )

    # Add network-derived node metrics.
    node_table["degree"] = node_table["cluster"].map(dict(graph.degree())).astype(int)
    node_table["weighted_degree"] = node_table["cluster"].map(
        dict(graph.degree(weight="weight"))
    ).astype(float)
    node_table["backbone_degree"] = node_table["cluster"].map(
        dict(backbone.degree())
    ).astype(int)

    # Rebuild graphs so exported node attributes include network metrics.
    graph, backbone = build_graphs(node_table, significant_edges)

    # Save pairwise and node tables.
    pairwise.to_csv(args.outdir / "pairwise_coupling_all.csv", index=False)
    significant_edges.to_csv(
        args.outdir / "pairwise_coupling_significant_positive_edges.csv",
        index=False,
    )
    node_table.to_csv(args.outdir / "clonal_coupling_nodes.csv", index=False)

    # Save complete matrices. Negative values remain here even though they are
    # not included in the positive-edge graph.
    matrix_to_dataframe(observed_score, cluster_order).to_csv(
        args.outdir / "observed_coupling_score_matrix.csv"
    )
    matrix_to_dataframe(null_mean, cluster_order).to_csv(
        args.outdir / "null_mean_coupling_score_matrix.csv"
    )
    matrix_to_dataframe(z_score, cluster_order).to_csv(
        args.outdir / "clonal_coupling_zscore_matrix.csv"
    )
    matrix_to_dataframe(shared_clones, cluster_order).to_csv(
        args.outdir / "shared_clone_count_matrix.csv"
    )
    matrix_to_dataframe(q_enrichment, cluster_order).to_csv(
        args.outdir / "clonal_coupling_enrichment_fdr_matrix.csv"
    )

    epsilon = 1e-9
    log2_enrichment_matrix = np.log2(
        (observed_score + epsilon) / (null_mean + epsilon)
    )
    np.fill_diagonal(log2_enrichment_matrix, 0.0)
    matrix_to_dataframe(log2_enrichment_matrix, cluster_order).to_csv(
        args.outdir / "clonal_coupling_log2_enrichment_matrix.csv"
    )

    # Graph files can be re-opened in Cytoscape or Gephi.
    nx.write_graphml(graph, args.outdir / "clonal_coupling_network.graphml")
    nx.write_graphml(backbone, args.outdir / "clonal_coupling_backbone.graphml")

    # Save layout coordinates for exact reuse elsewhere.
    pd.DataFrame({
        "cluster": list(positions.keys()),
        "x": [float(positions[node][0]) for node in positions],
        "y": [float(positions[node][1]) for node in positions],
    }).to_csv(args.outdir / "network_layout_coordinates.csv", index=False)

    plot_static_network(
        graph=graph,
        backbone=backbone,
        node_table=node_table,
        positions=positions,
        outdir=args.outdir,
    )
    plot_backbone_only(
        graph=graph,
        backbone=backbone,
        node_table=node_table,
        positions=positions,
        outdir=args.outdir,
    )
    plot_zscore_heatmap(
        z_score=z_score,
        cluster_order=cluster_order,
        outdir=args.outdir,
    )
    write_interactive_network(
        graph=graph,
        backbone=backbone,
        node_table=node_table,
        positions=positions,
        outdir=args.outdir,
    )

    n_isolates = nx.number_of_isolates(graph)

    print("\nAnalysis complete")
    print("-----------------")
    print(f"Significant positive edges: {len(significant_edges):,}")
    print(f"Backbone edges:             {backbone.number_of_edges():,}")
    print(f"Isolated nodes:             {n_isolates:,}")
    print(f"Results written to:         {args.outdir.resolve()}\n")

    if significant_edges.empty:
        print(
            "No edges passed the current thresholds. Inspect pairwise_coupling_all.csv "
            "and consider whether --min-shared-clones or --fdr is too strict."
        )


if __name__ == "__main__":
    main()
