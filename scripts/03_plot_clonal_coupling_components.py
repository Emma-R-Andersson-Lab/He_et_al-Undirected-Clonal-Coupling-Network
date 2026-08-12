#!/usr/bin/env python3
"""
Plot one publication-style PNG and SVG for each connected component of a
previously calculated clonal-coupling network.

This script does NOT recalculate coupling statistics and does NOT change which
edges are significant. It reads the standard outputs from
clonal_coupling_network.py:

    pairwise_coupling_significant_positive_edges.csv
    clonal_coupling_nodes.csv

For every connected component containing at least one edge, it produces:

    component_XX_nodesN_edgesM.png
    component_XX_nodesN_edgesM.svg

Meaning of the visualization:
    - Node: transcriptomic cluster
    - Node area: total cluster cell count, scaled globally across all components
    - Node label: cluster name and total cell count
    - Edge: significant positive clonal coupling already accepted upstream
    - Edge width: EDGE_WIDTH_FACTOR * coupling z-score
    - Edge label: z-score and number of shared clones by default
    - Dark edge: maximum-spanning-tree edge within that connected component
    - Grey edge: additional significant positive coupling

The same node-size scale and edge-width scale are used for every component,
which makes separate component figures directly comparable.

An optional --min-shared-clones plotting filter can further restrict the
already-significant edge table. Connected components and backbones are then
recalculated from the retained edges; coupling statistics are not recomputed.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import re
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like, to_rgb
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
import networkx as nx
import numpy as np
import pandas as pd


# =============================================================================
# USER-EDITABLE COLOUR CONFIGURATION
# =============================================================================

# For labels beginning with "c14_", the suffix is looked up here.
# Example: "c14_10+1" -> C14_CLUSTER_COLORS["10+1"]
C14_CLUSTER_COLORS: dict[str, str] = {
    "0": "#41A02C",     # grass green
    "10+1": "#E8590C",  # orange
    "11+9": "#109C6B",  # jade green
    "12": "#5FA8E0",    # light blue
    "13": "#DC257F",    # magenta
    "14": "#BE970D",    # dark gold
    "15": "#8E1F63",    # deep plum
    "16": "#F28EBA",    # pale pink
    "17": "#8C613C",    # brown
    "2+7": "#C22E42",   # deep red
    "3": "#8B5FE0",     # bright violet
    "4": "#6C63AC",     # muted violet
    "5": "#17798F",     # dark cyan
    "6_0": "#4356D0",   # blue
    "6_1": "#0E9AA7",   # teal
    "8": "#FCAF17",     # amber
}

# All other cluster labels are looked up directly here.
INTEGRATED_CLUSTER_COLORS: dict[str, str] = {
    "0": "#DA93E7",
    "1": "#DA70D6",
    "2": "#FD673A",
    "3": "#F28C22",
    "4": "#FF7F50",
    "5": "#7DC5C6",
    "6": "#F28C30",
    "7": "#0ABAB5",
    "8": "#BF40BD",
    "9": "#93C572",
    "10": "#D74826",
    "11": "#FF69B4",
    "12": "#FD673A",
    "13": "#DA70C7",
    "14": "#DA93D2",
    "15": "#FF7F80",
    "16": "#819D4D",
    "17": "#7CFC00",
    "18": "#EFBB60",
    "19": "#7DC5B1",
    "20": "#FFBB30",
    "21": "#D74129",
    "22": "#DA70C7",
    "23": "#673147",
    "24": "#93C572",
    "25": "#0F52BA",
    "26": "#673147",
    "27": "#f7d40a",
    "28": "#0F72BA",
    "29": "#BF40BF",
}

C14_PREFIX = "c14_"
DEFAULT_NODE_COLOR = "#B8BDC7"
BACKBONE_EDGE_COLOR = "#263238"
SECONDARY_EDGE_COLOR = "#9AA0A6"
NODE_BORDER_COLOR = "#172033"
EDGE_WIDTH_FACTOR = 0.35

# Keep text editable in vector graphics.
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.family"] = "DejaVu Sans"


# =============================================================================
# ARGUMENTS
# =============================================================================


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create one PNG and SVG per connected component from an existing "
            "clonal-coupling analysis output directory."
        )
    )
    parser.add_argument(
        "--analysis-dir",
        required=True,
        type=Path,
        help=(
            "Directory containing pairwise_coupling_significant_positive_edges.csv "
            "and clonal_coupling_nodes.csv."
        ),
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory for component PNG/SVG files and summary tables.",
    )
    parser.add_argument(
        "--edges-name",
        default="pairwise_coupling_significant_positive_edges.csv",
        help="Significant-edge CSV filename inside --analysis-dir.",
    )
    parser.add_argument(
        "--nodes-name",
        default="clonal_coupling_nodes.csv",
        help="Node CSV filename inside --analysis-dir.",
    )
    parser.add_argument(
        "--colors-csv",
        type=Path,
        default=None,
        help=(
            "Optional CSV with columns 'cluster' and 'color'. Full cluster-label "
            "matches override the built-in colour dictionaries."
        ),
    )
    parser.add_argument(
        "--node-size-column",
        default="n_total_cells",
        help=(
            "Node-table column used for node area. Falls back to n_barcoded_cells "
            "if unavailable. Default: n_total_cells."
        ),
    )
    parser.add_argument(
        "--min-shared-clones",
        type=int,
        default=0,
        help=(
            "Additional plotting-time filter for edge support. Keep only edges "
            "with at least this many shared clones before connected components "
            "and backbones are recalculated. Default: 0, which keeps every edge "
            "already present in the upstream significant-edge CSV."
        ),
    )
    parser.add_argument(
        "--edge-width-factor",
        type=float,
        default=EDGE_WIDTH_FACTOR,
        help=(
            "Direct multiplicative plotting factor for z-scores. Width = factor × z. "
            "Default: 0.35."
        ),
    )
    parser.add_argument(
        "--min-node-area",
        type=float,
        default=90.0,
        help=(
            "Minimum marker area in points squared for visibility. Default: 90. "
            "Set to 0 for exact area proportionality without a visibility floor."
        ),
    )
    parser.add_argument(
        "--max-node-area",
        type=float,
        default=2200.0,
        help="Marker area assigned to the largest node globally. Default: 2200.",
    )
    parser.add_argument(
        "--edge-labels",
        choices=("z", "z_shared", "full", "none"),
        default="z_shared",
        help=(
            "Edge annotation: z; z and shared clones; full z/shared clones/FDR; "
            "or none. Default: z_shared."
        ),
    )
    parser.add_argument(
        "--edge-label-decimals",
        type=int,
        default=2,
        help="Number of decimals shown for z-scores. Default: 2.",
    )
    parser.add_argument(
        "--node-labels-only",
        action="store_true",
        help=(
            "Keep only transcriptomic cluster names as direct plot annotations. "
            "This hides node-count annotations, backbone edge statistics, title, "
            "subtitle, and footer, while retaining the node-size and edge-width "
            "legends."
        ),
    )
    parser.add_argument(
        "--edge-routing",
        choices=("straight", "curved"),
        default="curved",
        help=(
            "Edge routing style. 'curved' gives every edge a deterministic arc "
            "and draws a white halo at crossings; 'straight' uses ordinary line "
            "segments. Default: curved."
        ),
    )
    parser.add_argument(
        "--edge-curvature",
        type=float,
        default=0.16,
        help=(
            "Base curvature for curved edge routing. Larger values separate nearby "
            "edges more strongly. Default: 0.16."
        ),
    )
    parser.add_argument(
        "--edge-halo-width",
        type=float,
        default=1.2,
        help=(
            "White outline added around curved edges, in points, so crossings remain "
            "visually separable. Set to 0 to disable. Default: 1.2."
        ),
    )
    parser.add_argument(
        "--layout-seed",
        type=int,
        default=1234,
        help="Random seed for reproducible component layouts. Default: 1234.",
    )
    parser.add_argument(
        "--layout-iterations",
        type=int,
        default=4000,
        help="Layout iterations. Default: 4000.",
    )
    parser.add_argument(
        "--layout-mode",
        choices=("community", "spring"),
        default="community",
        help=(
            "Layout strategy. 'community' detects coupling modules and separates "
            "them before arranging nodes within each module; 'spring' uses one "
            "global force-directed layout. Default: community."
        ),
    )
    parser.add_argument(
        "--community-spacing",
        type=float,
        default=4.5,
        help=(
            "Separation between detected communities in community layout. "
            "Larger values spread modules farther apart. Default: 4.5."
        ),
    )
    parser.add_argument(
        "--within-community-scale",
        type=float,
        default=0.95,
        help=(
            "Base radius used to spread nodes inside each detected community. "
            "Larger values spread nodes within modules. Default: 0.95."
        ),
    )
    parser.add_argument(
        "--include-isolates",
        action="store_true",
        help="Also make one-node figures for isolated clusters.",
    )
    parser.add_argument(
        "--png-dpi",
        type=int,
        default=400,
        help="PNG resolution. Default: 400 DPI.",
    )
    parser.add_argument(
        "--transparent",
        action="store_true",
        help="Save PNG and SVG with transparent backgrounds.",
    )
    parser.add_argument(
        "--c14-prefix",
        default=C14_PREFIX,
        help="Prefix identifying c14 subclusters. Default: c14_.",
    )
    return parser.parse_args()


# =============================================================================
# INPUT AND VALIDATION
# =============================================================================


def natural_key(value: str) -> list[object]:
    return [
        int(token) if token.isdigit() else token.lower()
        for token in re.split(r"(\d+)", str(value))
    ]


def load_color_overrides(path: Path | None) -> dict[str, str]:
    if path is None:
        return {}

    table = pd.read_csv(path, dtype=str)
    required = {"cluster", "color"}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(
            f"Colour override CSV is missing column(s): {', '.join(missing)}"
        )

    overrides: dict[str, str] = {}
    for row in table[["cluster", "color"]].dropna().itertuples(index=False):
        cluster = str(row.cluster).strip()
        color = str(row.color).strip()
        if not is_color_like(color):
            raise ValueError(
                f"Invalid colour {color!r} for cluster {cluster!r} in {path}."
            )
        overrides[cluster] = color

    return overrides


def resolve_node_color(
    cluster: str,
    overrides: dict[str, str],
    c14_prefix: str,
) -> tuple[str, bool]:
    """Return (colour, matched_known_mapping)."""
    cluster = str(cluster)

    if cluster in overrides:
        return overrides[cluster], True

    if cluster.startswith(c14_prefix):
        suffix = cluster[len(c14_prefix):]
        if suffix in C14_CLUSTER_COLORS:
            return C14_CLUSTER_COLORS[suffix], True
        return DEFAULT_NODE_COLOR, False

    if cluster in INTEGRATED_CLUSTER_COLORS:
        return INTEGRATED_CLUSTER_COLORS[cluster], True

    return DEFAULT_NODE_COLOR, False


def load_network_tables(
    edges_path: Path,
    nodes_path: Path,
    node_size_column: str,
) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    if not edges_path.exists():
        raise FileNotFoundError(f"Significant-edge file not found: {edges_path}")
    if not nodes_path.exists():
        raise FileNotFoundError(f"Node file not found: {nodes_path}")

    edges = pd.read_csv(edges_path)
    nodes = pd.read_csv(nodes_path)

    required_edges = {
        "cluster_1",
        "cluster_2",
        "z_score",
        "shared_clones",
        "q_enrichment",
    }
    missing_edges = sorted(required_edges - set(edges.columns))
    if missing_edges:
        raise ValueError(
            "Edge CSV is missing required column(s): " + ", ".join(missing_edges)
        )

    if "cluster" not in nodes.columns:
        raise ValueError("Node CSV must contain a 'cluster' column.")

    nodes = nodes.copy()
    edges = edges.copy()

    nodes["cluster"] = nodes["cluster"].astype(str).str.strip()
    edges["cluster_1"] = edges["cluster_1"].astype(str).str.strip()
    edges["cluster_2"] = edges["cluster_2"].astype(str).str.strip()

    if nodes["cluster"].duplicated().any():
        duplicates = nodes.loc[nodes["cluster"].duplicated(), "cluster"].tolist()
        raise ValueError(
            "Node CSV contains duplicate cluster rows, including: "
            + ", ".join(duplicates[:5])
        )

    for column in ("z_score", "shared_clones", "q_enrichment"):
        edges[column] = pd.to_numeric(edges[column], errors="coerce")

    invalid_edges = (
        edges["z_score"].isna()
        | edges["shared_clones"].isna()
        | edges["q_enrichment"].isna()
    )
    if invalid_edges.any():
        warnings.warn(
            f"Dropping {int(invalid_edges.sum()):,} edges with invalid numeric values."
        )
        edges = edges.loc[~invalid_edges].copy()

    # The upstream significant-edge file should contain positive edges only.
    nonpositive = edges["z_score"] <= 0
    if nonpositive.any():
        warnings.warn(
            f"Dropping {int(nonpositive.sum()):,} non-positive edges from the "
            "component visualisation."
        )
        edges = edges.loc[~nonpositive].copy()

    if node_size_column not in nodes.columns:
        if "n_barcoded_cells" not in nodes.columns:
            raise ValueError(
                f"Node-size column {node_size_column!r} is unavailable, and "
                "n_barcoded_cells is also absent."
            )
        warnings.warn(
            f"Node-size column {node_size_column!r} not found; using "
            "'n_barcoded_cells' instead."
        )
        node_size_column = "n_barcoded_cells"

    nodes[node_size_column] = pd.to_numeric(
        nodes[node_size_column], errors="coerce"
    ).fillna(0.0)

    negative_sizes = nodes[node_size_column] < 0
    if negative_sizes.any():
        raise ValueError(
            f"Node-size column {node_size_column!r} contains negative values."
        )

    node_set = set(nodes["cluster"])
    edge_nodes = set(edges["cluster_1"]) | set(edges["cluster_2"])
    missing_nodes = sorted(edge_nodes - node_set, key=natural_key)
    if missing_nodes:
        raise ValueError(
            "Edges reference cluster(s) absent from the node CSV: "
            + ", ".join(missing_nodes[:10])
        )

    return edges, nodes, node_size_column


def build_graph(edges: pd.DataFrame, nodes: pd.DataFrame) -> nx.Graph:
    graph = nx.Graph()

    for row in nodes.itertuples(index=False):
        attrs = row._asdict()
        cluster = str(attrs.pop("cluster"))
        clean_attrs = {
            key: value.item() if isinstance(value, np.generic) else value
            for key, value in attrs.items()
            if pd.notna(value)
        }
        graph.add_node(cluster, **clean_attrs)

    for row in edges.itertuples(index=False):
        attrs = row._asdict()
        u = str(attrs.pop("cluster_1"))
        v = str(attrs.pop("cluster_2"))
        z = float(attrs.get("z_score"))

        clean_attrs = {
            key: value.item() if isinstance(value, np.generic) else value
            for key, value in attrs.items()
            if pd.notna(value)
        }
        clean_attrs["z_score"] = z
        clean_attrs["weight"] = z
        graph.add_edge(u, v, **clean_attrs)

    return graph


# =============================================================================
# GLOBAL VISUAL SCALES
# =============================================================================


def make_node_area_function(
    nodes: pd.DataFrame,
    node_size_column: str,
    min_area: float,
    max_area: float,
):
    if min_area < 0:
        raise ValueError("--min-node-area must be non-negative.")
    if max_area <= 0:
        raise ValueError("--max-node-area must be greater than zero.")
    if min_area > max_area:
        raise ValueError("--min-node-area cannot exceed --max-node-area.")

    values = nodes[node_size_column].to_numpy(dtype=float)
    largest = float(np.max(values)) if values.size else 0.0
    if largest <= 0:
        largest = 1.0

    def area_from_count(count: float) -> float:
        proportional = (max(float(count), 0.0) / largest) * max_area
        return max(min_area, proportional)

    return area_from_count, largest


def edge_width(z_score: float, factor: float, largest_finite_z: float) -> float:
    z = float(z_score)
    if not np.isfinite(z):
        z = largest_finite_z
    if z <= 0:
        raise ValueError(f"Cannot draw a non-positive coupling edge: z={z}")
    return z * factor


def choose_node_size_legend_values(
    nodes: pd.DataFrame,
    node_size_column: str,
) -> list[int]:
    values = nodes[node_size_column].to_numpy(dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    if values.size == 0:
        return []

    candidates = np.quantile(values, [0.25, 0.50, 1.00])
    rounded: list[int] = []
    for value in candidates:
        integer = int(round(float(value)))
        if integer > 0 and integer not in rounded:
            rounded.append(integer)
    return rounded


def choose_edge_width_legend_values(graph: nx.Graph) -> list[float]:
    """Choose globally comparable z-scores for the edge-width legend."""
    values = np.asarray(
        [
            float(attrs["z_score"])
            for _, _, attrs in graph.edges(data=True)
            if np.isfinite(float(attrs["z_score"]))
            and float(attrs["z_score"]) > 0
        ],
        dtype=float,
    )
    if values.size == 0:
        return []

    candidates = np.quantile(values, [0.25, 0.50, 1.00])
    selected: list[float] = []
    seen: set[float] = set()
    for value in candidates:
        rounded = round(float(value), 2)
        if rounded > 0 and rounded not in seen:
            selected.append(rounded)
            seen.add(rounded)
    return selected


# =============================================================================
# LAYOUT AND LABEL HELPERS
# =============================================================================


def component_layout(
    graph: nx.Graph,
    seed: int,
    iterations: int,
    mode: str,
    community_spacing: float,
    within_community_scale: float,
) -> dict[str, np.ndarray]:
    """Create a readable layout without changing any network statistic.

    The community layout first detects densely coupled node groups, positions
    those groups apart from one another, and then runs a separate spring layout
    inside each group. For layout only, z-score weights are log-compressed so a
    very large z-score cannot collapse several nodes onto the same location.
    """
    n_nodes = graph.number_of_nodes()

    if n_nodes == 1:
        node = next(iter(graph.nodes()))
        return {node: np.array([0.0, 0.0])}

    if n_nodes == 2:
        nodes = sorted(graph.nodes(), key=natural_key)
        return {
            nodes[0]: np.array([-1.20, 0.0]),
            nodes[1]: np.array([1.20, 0.0]),
        }

    # Make a layout-only copy. This changes positions only, never edge values,
    # significance, widths, labels, or the maximum-spanning-tree backbone.
    layout_graph = graph.copy()
    for u, v, attrs in layout_graph.edges(data=True):
        z = max(float(attrs.get("z_score", attrs.get("weight", 1.0))), 0.0)
        attrs["layout_weight"] = math.log1p(z)

    def global_spring() -> dict[str, np.ndarray]:
        positions = nx.spring_layout(
            layout_graph,
            seed=seed,
            weight="layout_weight",
            k=3.0 / math.sqrt(n_nodes),
            iterations=iterations,
            scale=2.8,
        )
        return {
            node: np.asarray(position, dtype=float)
            for node, position in positions.items()
        }

    if mode == "spring":
        return global_spring()

    communities = list(
        nx.community.greedy_modularity_communities(
            layout_graph,
            weight="layout_weight",
        )
    )

    # Deterministic ordering helps keep output reproducible.
    communities = sorted(
        communities,
        key=lambda members: (
            -len(members),
            natural_key(min((str(node) for node in members), key=natural_key)),
        ),
    )

    if len(communities) <= 1:
        return global_spring()

    node_to_community: dict[str, int] = {}
    for community_index, members in enumerate(communities):
        for node in members:
            node_to_community[str(node)] = community_index

    community_graph = nx.Graph()
    community_graph.add_nodes_from(range(len(communities)))

    for u, v, attrs in layout_graph.edges(data=True):
        community_u = node_to_community[str(u)]
        community_v = node_to_community[str(v)]
        if community_u == community_v:
            continue

        weight = float(attrs.get("layout_weight", 1.0))
        if community_graph.has_edge(community_u, community_v):
            community_graph[community_u][community_v]["weight"] += weight
        else:
            community_graph.add_edge(community_u, community_v, weight=weight)

    n_communities = len(communities)
    if n_communities == 2:
        community_centres = {
            0: np.array([-community_spacing / 2.0, 0.0]),
            1: np.array([community_spacing / 2.0, 0.0]),
        }
    else:
        community_centres = nx.spring_layout(
            community_graph,
            seed=seed,
            weight="weight",
            k=2.4 / math.sqrt(n_communities),
            iterations=iterations,
            scale=community_spacing,
        )

    positions: dict[str, np.ndarray] = {}

    for community_index, members in enumerate(communities):
        members_sorted = sorted((str(node) for node in members), key=natural_key)
        subgraph = layout_graph.subgraph(members_sorted).copy()
        community_size = len(members_sorted)
        centre = np.asarray(community_centres[community_index], dtype=float)

        if community_size == 1:
            local_positions = {members_sorted[0]: np.array([0.0, 0.0])}
        elif community_size == 2:
            local_positions = {
                members_sorted[0]: np.array([-0.65, 0.0]),
                members_sorted[1]: np.array([0.65, 0.0]),
            }
        else:
            local_scale = within_community_scale * (
                0.85 + 0.22 * math.sqrt(community_size)
            )
            local_positions = nx.spring_layout(
                subgraph,
                seed=seed + 1009 * (community_index + 1),
                weight="layout_weight",
                k=2.7 / math.sqrt(community_size),
                iterations=iterations,
                scale=local_scale,
            )

        for node, local_position in local_positions.items():
            positions[str(node)] = centre + np.asarray(local_position, dtype=float)

    return positions


def display_cluster_label(cluster: str, c14_prefix: str) -> str:
    cluster = str(cluster)
    if cluster.startswith(c14_prefix):
        return f"{c14_prefix}\n{cluster[len(c14_prefix):]}"
    return cluster


def contrasting_text_color(hex_color: str) -> str:
    red, green, blue = to_rgb(hex_color)
    # Relative luminance approximation for readable node labels.
    luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
    return "#111111" if luminance > 0.58 else "#FFFFFF"


def format_edge_label(
    attrs: dict[str, object],
    mode: str,
    decimals: int,
) -> str:
    if mode == "none":
        return ""

    z = float(attrs["z_score"])
    z_text = f"z={z:.{decimals}f}"

    if mode == "z":
        return z_text

    shared = int(float(attrs.get("shared_clones", 0)))
    if mode == "z_shared":
        return f"{z_text}\nclones={shared}"

    q = float(attrs.get("q_enrichment", np.nan))
    q_text = f"q={q:.2g}" if np.isfinite(q) else "q=NA"
    return f"{z_text}\nclones={shared}\n{q_text}"


def sanitise_filename(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    return cleaned.strip("._") or "component"


def deterministic_edge_curvature(
    u: str,
    v: str,
    base_curvature: float,
    is_backbone: bool,
) -> float:
    """Return a stable signed curvature for one undirected edge."""
    first, second = sorted((str(u), str(v)), key=natural_key)
    digest = hashlib.sha1(f"{first}|{second}".encode("utf-8")).digest()
    sign = -1.0 if digest[0] % 2 else 1.0
    level = 1.0 + 0.18 * (digest[1] % 4)
    backbone_multiplier = 0.65 if is_backbone else 1.0
    return sign * base_curvature * level * backbone_multiplier


def curved_midpoint(
    start: np.ndarray,
    end: np.ndarray,
    curvature: float,
) -> np.ndarray:
    """Approximate the visible midpoint of a Matplotlib arc3 edge."""
    start = np.asarray(start, dtype=float)
    end = np.asarray(end, dtype=float)
    midpoint = (start + end) / 2.0
    delta = end - start
    distance = float(np.linalg.norm(delta))
    if distance <= 0:
        return midpoint
    perpendicular = np.array([-delta[1], delta[0]], dtype=float) / distance
    return midpoint + perpendicular * curvature * distance * 0.55


def draw_routed_edges(
    ax: plt.Axes,
    graph: nx.Graph,
    positions: dict[str, np.ndarray],
    edges: list[tuple[str, str]],
    node_areas: dict[str, float],
    color: str,
    alpha: float,
    zorder: float,
    width_factor: float,
    largest_finite_z: float,
    routing: str,
    base_curvature: float,
    halo_width: float,
    backbone_edge_keys: set[frozenset[str]],
) -> dict[frozenset[str], float]:
    """Draw routed edges and return the curvature used for each edge."""
    curvature_by_edge: dict[frozenset[str], float] = {}

    if routing == "straight":
        nx.draw_networkx_edges(
            graph,
            positions,
            edgelist=edges,
            width=[
                edge_width(
                    graph[u][v]["z_score"],
                    width_factor,
                    largest_finite_z,
                )
                for u, v in edges
            ],
            edge_color=color,
            alpha=alpha,
            ax=ax,
        )
        for u, v in edges:
            curvature_by_edge[frozenset((u, v))] = 0.0
        return curvature_by_edge

    for u, v in edges:
        key = frozenset((u, v))
        is_backbone = key in backbone_edge_keys
        curvature = deterministic_edge_curvature(
            u=u,
            v=v,
            base_curvature=base_curvature,
            is_backbone=is_backbone,
        )
        curvature_by_edge[key] = curvature

        width = edge_width(
            graph[u][v]["z_score"],
            width_factor,
            largest_finite_z,
        )
        shrink_a = math.sqrt(max(node_areas[u], 0.0) / math.pi)
        shrink_b = math.sqrt(max(node_areas[v], 0.0) / math.pi)

        if halo_width > 0:
            halo = FancyArrowPatch(
                posA=positions[u],
                posB=positions[v],
                arrowstyle="-",
                connectionstyle=f"arc3,rad={curvature}",
                linewidth=width + 2.0 * halo_width,
                color="white",
                alpha=0.96,
                shrinkA=shrink_a,
                shrinkB=shrink_b,
                capstyle="round",
                joinstyle="round",
                zorder=zorder - 0.10,
            )
            ax.add_patch(halo)

        edge_patch = FancyArrowPatch(
            posA=positions[u],
            posB=positions[v],
            arrowstyle="-",
            connectionstyle=f"arc3,rad={curvature}",
            linewidth=width,
            color=color,
            alpha=alpha,
            shrinkA=shrink_a,
            shrinkB=shrink_b,
            capstyle="round",
            joinstyle="round",
            zorder=zorder,
        )
        ax.add_patch(edge_patch)

    return curvature_by_edge


# =============================================================================
# COMPONENT PLOTTING
# =============================================================================


def plot_component(
    component_id: int,
    component_graph: nx.Graph,
    nodes_indexed: pd.DataFrame,
    node_size_column: str,
    area_from_count,
    largest_finite_z: float,
    size_legend_values: list[int],
    edge_legend_values: list[float],
    colors: dict[str, str],
    args: argparse.Namespace,
) -> dict[str, object]:
    component_nodes = sorted(component_graph.nodes(), key=natural_key)
    component_graph = component_graph.subgraph(component_nodes).copy()
    backbone = nx.maximum_spanning_tree(component_graph, weight="weight")

    positions = component_layout(
        component_graph,
        seed=args.layout_seed + component_id,
        iterations=args.layout_iterations,
        mode=args.layout_mode,
        community_spacing=args.community_spacing,
        within_community_scale=args.within_community_scale,
    )

    backbone_edges = {frozenset(edge) for edge in backbone.edges()}
    secondary_edges = [
        edge
        for edge in component_graph.edges()
        if frozenset(edge) not in backbone_edges
    ]

    node_counts = {
        node: float(nodes_indexed.loc[node, node_size_column])
        for node in component_nodes
    }
    node_areas = {
        node: area_from_count(node_counts[node])
        for node in component_nodes
    }

    n_nodes = component_graph.number_of_nodes()
    n_edges = component_graph.number_of_edges()
    total_cells = int(round(sum(node_counts.values())))

    # Figure size grows sublinearly with component complexity.
    base_side = 7.6 + 1.15 * math.sqrt(max(n_nodes, 1))
    density_bonus = min(4.0, 0.08 * max(n_edges - n_nodes, 0))
    side = min(20.0, base_side + density_bonus)

    fig, ax = plt.subplots(figsize=(side, side))
    fig.patch.set_alpha(0.0 if args.transparent else 1.0)
    ax.set_facecolor("none" if args.transparent else "white")

    curvature_by_edge: dict[frozenset[str], float] = {}

    if secondary_edges:
        curvature_by_edge.update(
            draw_routed_edges(
                ax=ax,
                graph=component_graph,
                positions=positions,
                edges=list(secondary_edges),
                node_areas=node_areas,
                color=SECONDARY_EDGE_COLOR,
                alpha=0.32,
                zorder=1.0,
                width_factor=args.edge_width_factor,
                largest_finite_z=largest_finite_z,
                routing=args.edge_routing,
                base_curvature=args.edge_curvature,
                halo_width=args.edge_halo_width,
                backbone_edge_keys=backbone_edges,
            )
        )

    if backbone.number_of_edges() > 0:
        curvature_by_edge.update(
            draw_routed_edges(
                ax=ax,
                graph=component_graph,
                positions=positions,
                edges=list(backbone.edges()),
                node_areas=node_areas,
                color=BACKBONE_EDGE_COLOR,
                alpha=0.88,
                zorder=2.0,
                width_factor=args.edge_width_factor,
                largest_finite_z=largest_finite_z,
                routing=args.edge_routing,
                base_curvature=args.edge_curvature,
                halo_width=args.edge_halo_width,
                backbone_edge_keys=backbone_edges,
            )
        )

    nx.draw_networkx_nodes(
        component_graph,
        positions,
        nodelist=component_nodes,
        node_size=[node_areas[node] for node in component_nodes],
        node_color=[colors[node] for node in component_nodes],
        edgecolors=NODE_BORDER_COLOR,
        linewidths=1.8,
        alpha=0.98,
        ax=ax,
    )

    # Cluster names are placed inside nodes; total-cell counts are placed just
    # below each node in a small white-backed annotation.
    for node in component_nodes:
        x, y = positions[node]
        node_color = colors[node]
        ax.text(
            x,
            y,
            display_cluster_label(node, args.c14_prefix),
            ha="center",
            va="center",
            fontsize=6.8 if n_nodes <= 18 else 5.8,
            fontweight="bold",
            color=contrasting_text_color(node_color),
            zorder=5,
        )

        if not args.node_labels_only:
            vertical_offset = -(math.sqrt(node_areas[node]) / 2.0 + 4.0)
            ax.annotate(
                f"N={int(round(node_counts[node])):,}",
                xy=(x, y),
                xytext=(0, vertical_offset),
                textcoords="offset points",
                ha="center",
                va="top",
                fontsize=5.5,
                color="#20242A",
                bbox={
                    "boxstyle": "round,pad=0.10",
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.78,
                },
                zorder=6,
            )

    # Draw all significant edges, but annotate statistics only on backbone
    # edges. Labels are shifted onto the corresponding curved route.
    if (
        not args.node_labels_only
        and args.edge_labels != "none"
        and backbone.number_of_edges() > 0
    ):
        for u, v in backbone.edges():
            key = frozenset((u, v))
            curvature = curvature_by_edge.get(key, 0.0)
            label_position = curved_midpoint(
                positions[u],
                positions[v],
                curvature=curvature,
            )
            ax.text(
                float(label_position[0]),
                float(label_position[1]),
                format_edge_label(
                    component_graph[u][v],
                    mode=args.edge_labels,
                    decimals=args.edge_label_decimals,
                ),
                ha="center",
                va="center",
                fontsize=4.6,
                color="#252A31",
                bbox={
                    "boxstyle": "round,pad=0.08",
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.76,
                },
                zorder=7,
            )

    finite_z = [
        float(attrs["z_score"])
        for _, _, attrs in component_graph.edges(data=True)
        if np.isfinite(float(attrs["z_score"]))
    ]
    z_min = min(finite_z) if finite_z else float("nan")
    z_max = max(finite_z) if finite_z else float("nan")

    if not args.node_labels_only:
        fig.suptitle(
            f"Clonal coupling network — connected component {component_id:02d}",
            fontsize=13,
            fontweight="bold",
            y=0.982,
        )
        subtitle = (
            f"{n_nodes} clusters  |  {n_edges} significant positive couplings  |  "
            f"total cells={total_cells:,}"
        )
        if finite_z:
            subtitle += f"  |  z range={z_min:.2f}–{z_max:.2f}"
        fig.text(0.5, 0.947, subtitle, ha="center", va="center", fontsize=8.5)

    edge_type_handles: list[object] = [
        Line2D(
            [0],
            [0],
            color=BACKBONE_EDGE_COLOR,
            linewidth=3.0,
            label="Maximum-spanning-tree backbone",
        ),
        Line2D(
            [0],
            [0],
            color=SECONDARY_EDGE_COLOR,
            linewidth=3.0,
            alpha=0.55,
            label="Additional significant coupling",
        ),
    ]

    edge_width_handles: list[object] = []
    for z_value in edge_legend_values:
        edge_width_handles.append(
            Line2D(
                [0],
                [0],
                color="#5B626B",
                linewidth=edge_width(
                    z_value,
                    args.edge_width_factor,
                    largest_finite_z,
                ),
                solid_capstyle="round",
                label=f"z={z_value:g}",
            )
        )

    node_size_handles: list[object] = []
    for count in size_legend_values:
        node_size_handles.append(
            Line2D(
                [0],
                [0],
                linestyle="none",
                marker="o",
                markersize=math.sqrt(area_from_count(count)),
                markerfacecolor="#D8DCE2",
                markeredgecolor=NODE_BORDER_COLOR,
                markeredgewidth=1.2,
                label=f"{count:,} cells",
            )
        )

    edge_legend = fig.legend(
        handles=edge_type_handles + edge_width_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.105),
        ncol=max(1, len(edge_type_handles + edge_width_handles)),
        frameon=False,
        fontsize=7.0,
        handlelength=3.2,
        handleheight=1.4,
        columnspacing=1.4,
        title="Edge type and z-score width scale",
        title_fontsize=7.2,
    )
    edge_legend.set_zorder(10)

    if node_size_handles:
        node_legend = fig.legend(
            handles=node_size_handles,
            loc="lower center",
            bbox_to_anchor=(0.5, 0.018),
            ncol=len(node_size_handles),
            frameon=False,
            fontsize=7.0,
            handlelength=2.6,
            handleheight=3.8,
            columnspacing=2.0,
        )
        node_legend.set_zorder(10)

    if not args.node_labels_only:
        footer = (
            f"Edge width = {args.edge_width_factor:g} × z-score for every component. "
            f"Node area uses {node_size_column!r} on one global scale. "
            "Dark/grey colour does not alter edge thickness. "
            f"Edge routing={args.edge_routing}."
        )
        if args.min_node_area > 0:
            footer += " A minimum node area is applied only for visibility."
        fig.text(0.5, 0.018, footer, ha="center", va="bottom", fontsize=6.6)

    ax.margins(0.30)
    ax.set_aspect("equal")
    ax.axis("off")
    if args.node_labels_only:
        fig.subplots_adjust(top=0.97, bottom=0.28, left=0.04, right=0.96)
    else:
        fig.subplots_adjust(top=0.91, bottom=0.28, left=0.04, right=0.96)

    file_base = sanitise_filename(
        f"component_{component_id:02d}_nodes{n_nodes}_edges{n_edges}"
    )
    png_path = args.outdir / f"{file_base}.png"
    svg_path = args.outdir / f"{file_base}.svg"

    save_kwargs = {
        "bbox_inches": "tight",
        "facecolor": "none" if args.transparent else "white",
        "transparent": args.transparent,
    }
    fig.savefig(png_path, dpi=args.png_dpi, **save_kwargs)
    fig.savefig(svg_path, **save_kwargs)
    plt.close(fig)

    return {
        "component_id": component_id,
        "n_nodes": n_nodes,
        "n_edges": n_edges,
        "n_backbone_edges": backbone.number_of_edges(),
        "total_cells": total_cells,
        "minimum_z": z_min,
        "maximum_z": z_max,
        "nodes": ";".join(component_nodes),
        "png": png_path.name,
        "svg": svg_path.name,
    }


# =============================================================================
# MAIN
# =============================================================================


def main() -> None:
    args = parse_args()

    if args.min_shared_clones < 0:
        raise ValueError("--min-shared-clones must be zero or greater.")
    if args.edge_width_factor <= 0:
        raise ValueError("--edge-width-factor must be greater than zero.")
    if args.edge_label_decimals < 0:
        raise ValueError("--edge-label-decimals must be non-negative.")
    if args.edge_curvature < 0:
        raise ValueError("--edge-curvature must be zero or greater.")
    if args.edge_halo_width < 0:
        raise ValueError("--edge-halo-width must be zero or greater.")
    if args.layout_iterations <= 0:
        raise ValueError("--layout-iterations must be greater than zero.")
    if args.community_spacing <= 0:
        raise ValueError("--community-spacing must be greater than zero.")
    if args.within_community_scale <= 0:
        raise ValueError("--within-community-scale must be greater than zero.")

    args.outdir.mkdir(parents=True, exist_ok=True)

    edges_path = args.analysis_dir / args.edges_name
    nodes_path = args.analysis_dir / args.nodes_name

    edges, nodes, node_size_column = load_network_tables(
        edges_path=edges_path,
        nodes_path=nodes_path,
        node_size_column=args.node_size_column,
    )

    input_edge_count = len(edges)
    if args.min_shared_clones > 0:
        edges = edges.loc[
            edges["shared_clones"] >= args.min_shared_clones
        ].copy()
        edges = edges.reset_index(drop=True)

    graph = build_graph(edges, nodes)

    overrides = load_color_overrides(args.colors_csv)
    colors: dict[str, str] = {}
    unmatched: list[str] = []
    for node in graph.nodes():
        color, matched = resolve_node_color(
            node,
            overrides=overrides,
            c14_prefix=args.c14_prefix,
        )
        colors[node] = color
        if not matched:
            unmatched.append(node)

    if unmatched:
        unmatched_table = pd.DataFrame(
            {
                "cluster": sorted(unmatched, key=natural_key),
                "assigned_fallback_color": DEFAULT_NODE_COLOR,
            }
        )
        unmatched_path = args.outdir / "unmatched_cluster_colors.csv"
        unmatched_table.to_csv(unmatched_path, index=False)
        warnings.warn(
            f"{len(unmatched):,} cluster(s) had no colour mapping and were drawn "
            f"in {DEFAULT_NODE_COLOR}. See {unmatched_path}."
        )

    area_from_count, largest_node_count = make_node_area_function(
        nodes=nodes,
        node_size_column=node_size_column,
        min_area=args.min_node_area,
        max_area=args.max_node_area,
    )

    finite_positive_z = [
        float(attrs["z_score"])
        for _, _, attrs in graph.edges(data=True)
        if np.isfinite(float(attrs["z_score"]))
        and float(attrs["z_score"]) > 0
    ]
    largest_finite_z = max(finite_positive_z, default=1.0)

    components = list(nx.connected_components(graph))
    components.sort(
        key=lambda component: (
            -len(component),
            -graph.subgraph(component).number_of_edges(),
            [natural_key(node) for node in sorted(component, key=natural_key)],
        )
    )

    isolates = sorted(list(nx.isolates(graph)), key=natural_key)
    if isolates:
        isolate_rows = nodes.loc[nodes["cluster"].isin(isolates)].copy()
        isolate_rows.to_csv(args.outdir / "isolated_nodes.csv", index=False)

    plotted_components: list[set[str]] = []
    for component in components:
        subgraph = graph.subgraph(component)
        if subgraph.number_of_edges() == 0 and not args.include_isolates:
            continue
        plotted_components.append(set(component))

    plotted_node_names = set().union(*plotted_components) if plotted_components else set()
    plotted_nodes_table = nodes.loc[nodes["cluster"].isin(plotted_node_names)].copy()
    size_legend_values = choose_node_size_legend_values(
        plotted_nodes_table if not plotted_nodes_table.empty else nodes,
        node_size_column=node_size_column,
    )
    edge_legend_values = choose_edge_width_legend_values(graph)

    nodes_indexed = nodes.set_index("cluster", drop=False)
    manifest_rows: list[dict[str, object]] = []

    for component_id, component in enumerate(plotted_components, start=1):
        component_graph = graph.subgraph(component).copy()
        manifest_rows.append(
            plot_component(
                component_id=component_id,
                component_graph=component_graph,
                nodes_indexed=nodes_indexed,
                node_size_column=node_size_column,
                area_from_count=area_from_count,
                largest_finite_z=largest_finite_z,
                size_legend_values=size_legend_values,
                edge_legend_values=edge_legend_values,
                colors=colors,
                args=args,
            )
        )

    manifest = pd.DataFrame(manifest_rows)
    manifest.to_csv(args.outdir / "connected_component_manifest.csv", index=False)

    print("\nConnected-component plotting complete")
    print("------------------------------------")
    print(f"Input nodes:                 {graph.number_of_nodes():,}")
    print(f"Upstream significant edges:  {input_edge_count:,}")
    print(f"Minimum shared clones:       {args.min_shared_clones:,}")
    print(f"Edges after clone filter:    {graph.number_of_edges():,}")
    print(f"Connected components total:  {len(components):,}")
    print(f"Components plotted:          {len(plotted_components):,}")
    print(f"Isolated nodes:              {len(isolates):,}")
    print(f"Node-size column:            {node_size_column}")
    print(f"Largest node count:          {largest_node_count:,.0f}")
    print(f"Edge-width rule:             width = {args.edge_width_factor:g} × z")
    print(f"Edge routing:                {args.edge_routing}")
    print(f"Results written to:          {args.outdir.resolve()}\n")

    if not plotted_components:
        print(
            "No connected component with at least one edge was available to plot. "
            "Use --include-isolates to draw one-node figures."
        )


if __name__ == "__main__":
    main()
