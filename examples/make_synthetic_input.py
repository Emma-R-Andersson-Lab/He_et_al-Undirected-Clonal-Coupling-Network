#!/usr/bin/env python3
"""
Generate a small synthetic dataset in the format produced by
scripts/01_export_lineage_input.R.

This exists so the pipeline can be run, tested and demonstrated without access
to the real sequencing data. It is a positive control, not a simulation of any
real biology: clones are drawn from three overlapping "fate modules", so the
analysis step should recover exactly the within-module cluster pairs as
significant edges and leave the remaining clusters unconnected.

    python examples/make_synthetic_input.py --outdir example_input

Writes lineage_cells.csv and cluster_summary.csv.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# Clusters 0-2, 3-5 and 6-8 each form a coupled module.
# Clusters 9-11 receive clones drawn independently of module, so they should
# end up with no significant edges.
MODULES = {
    0: [0, 1, 2],
    1: [3, 4, 5],
    2: [6, 7, 8],
}
MODULE_WEIGHTS = [0.40, 0.35, 0.25]

N_CLUSTERS = 12
N_EMBRYOS = 4
N_CLONES_PER_EMBRYO = 150
MODULE_ENRICHMENT = 25.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path("example_input"),
        help="Directory for the generated CSV files. Default: example_input.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed. Default: 42.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    rows: list[tuple[str, str, str, str, str]] = []

    for embryo_index in range(N_EMBRYOS):
        embryo = f"E{embryo_index + 1}"

        # Uneven cluster abundances: this is exactly the confound that the
        # within-embryo permutation null has to absorb.
        abundance = rng.dirichlet(np.ones(N_CLUSTERS) * 0.8)

        for clone_index in range(N_CLONES_PER_EMBRYO):
            n_cells = int(rng.integers(1, 12))

            module = rng.choice(len(MODULES), p=MODULE_WEIGHTS)
            allowed = np.zeros(N_CLUSTERS, dtype=bool)
            allowed[MODULES[int(module)]] = True

            probability = np.where(allowed, abundance * MODULE_ENRICHMENT, abundance)
            probability = probability / probability.sum()

            clusters = rng.choice(N_CLUSTERS, n_cells, p=probability)

            clone_id = f"clone{clone_index}"
            for cell_index, cluster in enumerate(clusters):
                rows.append((
                    f"{embryo}_{clone_id}_{cell_index}",
                    embryo,
                    clone_id,
                    f"{embryo}::{clone_id}",
                    str(cluster),
                ))

    lineage_cells = pd.DataFrame(
        rows,
        columns=["cell_id", "embryo", "clone_id", "clone_uid", "cluster"],
    )

    # Node size uses n_total_cells, i.e. all clustered cells including the
    # barcode-negative ones. Simulate a barcoding efficiency of ~35-65%.
    barcoded = (
        lineage_cells.groupby("cluster")
        .agg(
            n_barcoded_cells=("cell_id", "size"),
            n_unique_clones=("clone_uid", "nunique"),
            n_embryos=("embryo", "nunique"),
        )
        .reset_index()
    )

    efficiency = rng.uniform(0.35, 0.65, len(barcoded))
    barcoded["n_total_cells"] = np.ceil(
        barcoded["n_barcoded_cells"] / efficiency
    ).astype(int)
    barcoded["barcode_fraction"] = (
        barcoded["n_barcoded_cells"] / barcoded["n_total_cells"]
    )

    barcoded["cluster"] = barcoded["cluster"].astype(str)
    barcoded = barcoded.sort_values(
        "cluster", key=lambda s: s.astype(int)
    ).reset_index(drop=True)
    barcoded.insert(1, "cluster_order", np.arange(1, len(barcoded) + 1))

    cluster_summary = barcoded[[
        "cluster",
        "cluster_order",
        "n_total_cells",
        "n_barcoded_cells",
        "n_unique_clones",
        "n_embryos",
        "barcode_fraction",
    ]]

    lineage_cells.to_csv(args.outdir / "lineage_cells.csv", index=False)
    cluster_summary.to_csv(args.outdir / "cluster_summary.csv", index=False)

    expected = sorted(
        f"{a}-{b}"
        for members in MODULES.values()
        for i, a in enumerate(members)
        for b in members[i + 1:]
    )

    print(f"Wrote {args.outdir.resolve()}")
    print(f"  barcode-positive cells : {len(lineage_cells):,}")
    print(f"  embryo-specific clones : {lineage_cells.clone_uid.nunique():,}")
    print(f"  clusters               : {lineage_cells.cluster.nunique()}")
    print(f"  embryos                : {lineage_cells.embryo.nunique()}")
    print()
    print(f"Expected significant edges ({len(expected)} within-module pairs):")
    print("  " + ", ".join(expected))
    print("Clusters 9, 10 and 11 should remain unconnected.")


if __name__ == "__main__":
    main()
