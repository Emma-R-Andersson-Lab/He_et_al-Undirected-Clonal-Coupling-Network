#!/usr/bin/env Rscript
# =============================================================================
# 01_export_lineage_input.R
#
# Step 1 of the undirected clonal-coupling network pipeline.
#
# Exports the cell-level tables that the Python analysis step
# (02_clonal_coupling_network.py) consumes, starting from a Seurat object that
# carries three pieces of metadata per cell:
#
#     - a transcriptomic cluster label,
#     - a lineage barcode / clone identifier,
#     - a biological replicate (embryo) identifier.
#
# Written files (all in --outdir):
#
#     lineage_cells.csv         one row per barcode-positive cell; the only
#                               file strictly required downstream
#     cluster_summary.csv       per-cluster node annotation, including
#                               n_total_cells over ALL clustered cells, which
#                               is what sizes the network nodes
#     clone_summary.csv         per-clone size and cluster occupancy
#     clone_cluster_counts.csv  the clone-by-cluster count matrix in long form
#     embryo_summary.csv        per-replicate totals
#
# The same barcode sequence observed in two different embryos is treated as two
# independent clones. This is enforced here by building clone_uid as
# "<embryo>::<clone_id>", and it is re-validated by the Python step.
#
# Usage
# -----
# Interactively, with a Seurat object already in the session:
#
#     source("scripts/01_export_lineage_input.R")
#     export_lineage_input(
#       obj         = my_seurat_object,
#       cluster_col = "res0_6_valve_subset_c14_refined",
#       clone_col   = "cloneid_prefixed",
#       embryo_col  = "embryo",
#       outdir      = "input"
#     )
#
# From the command line, against a saved object:
#
#     Rscript scripts/01_export_lineage_input.R \
#       --rds seurat_object.rds \
#       --cluster-col res0_6_valve_subset_c14_refined \
#       --clone-col cloneid_prefixed \
#       --embryo-col embryo \
#       --outdir input
#
# Set --cluster-col IDENT to use the object's active identities instead of a
# metadata column.
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(tidyr)
})

# readr::write_csv is used when available because it does not quote fields
# unnecessarily, but it is not a hard requirement.
write_csv_compat <- function(x, path) {
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_csv(x, path, na = "")
  } else {
    utils::write.csv(x, path, row.names = FALSE, na = "")
  }
  invisible(path)
}


# =============================================================================
# DEFAULTS
# =============================================================================

# Values in the clone column that mean "this cell carries no usable barcode".
# Do not add "0" unless 0 genuinely encodes an unbarcoded cell in your object.
DEFAULT_INVALID_CLONE_VALUES <- c(
  "", "NA", "NaN", "nan", "None", "none",
  "unassigned", "Unassigned", "no_clone", "NoClone"
)

# Values treated as missing in the embryo and cluster columns.
DEFAULT_INVALID_LABEL_VALUES <- c("", "NA", "NaN", "nan")


# =============================================================================
# EXPORTER
# =============================================================================

#' Export clonal-coupling input tables from a Seurat object.
#'
#' @param obj A Seurat object.
#' @param cluster_col Metadata column holding the transcriptomic cluster label,
#'   or the literal string "IDENT" to use Idents(obj).
#' @param clone_col Metadata column holding the clone / lineage barcode ID.
#' @param embryo_col Metadata column holding the biological replicate ID.
#' @param outdir Directory to write the CSV files into. Created if absent.
#' @param invalid_clone_values Character vector of clone-column values to treat
#'   as missing.
#' @param verbose Print a summary report when TRUE.
#'
#' @return Invisibly, a named list of the tables that were written.
export_lineage_input <- function(obj,
                                 cluster_col,
                                 clone_col,
                                 embryo_col,
                                 outdir,
                                 invalid_clone_values = DEFAULT_INVALID_CLONE_VALUES,
                                 verbose = TRUE) {

  # ---------------------------------------------------------------------------
  # Validation
  # ---------------------------------------------------------------------------
  if (!inherits(obj, "Seurat")) {
    stop("`obj` must be a Seurat object.", call. = FALSE)
  }

  meta <- obj[[]] %>%
    rownames_to_column("cell_id")

  required_cols <- c(clone_col, embryo_col)
  if (cluster_col != "IDENT") {
    required_cols <- c(required_cols, cluster_col)
  }

  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required metadata column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (cluster_col == "IDENT") {
    cluster_values <- as.character(Idents(obj)[meta$cell_id])
  } else {
    cluster_values <- as.character(meta[[cluster_col]])
  }

  # ---------------------------------------------------------------------------
  # Cell-level table
  # ---------------------------------------------------------------------------
  all_cells <- tibble(
    cell_id  = as.character(meta$cell_id),
    embryo   = trimws(as.character(meta[[embryo_col]])),
    clone_id = trimws(as.character(meta[[clone_col]])),
    cluster  = trimws(cluster_values)
  ) %>%
    mutate(
      clone_id = if_else(
        is.na(clone_id) | clone_id %in% invalid_clone_values,
        NA_character_,
        clone_id
      ),
      embryo = if_else(
        is.na(embryo) | embryo %in% DEFAULT_INVALID_LABEL_VALUES,
        NA_character_,
        embryo
      ),
      cluster = if_else(
        is.na(cluster) | cluster %in% DEFAULT_INVALID_LABEL_VALUES,
        NA_character_,
        cluster
      )
    )

  if (anyDuplicated(all_cells$cell_id) > 0) {
    stop(
      "Duplicate cell IDs were found. The exporter requires exactly one row ",
      "per cell.",
      call. = FALSE
    )
  }

  # Every clustered cell contributes to node size, barcoded or not.
  all_clustered_cells <- all_cells %>%
    filter(!is.na(cluster))

  # Only barcode-positive cells with a replicate and a cluster enter the
  # coupling statistics.
  lineage_cells <- all_cells %>%
    filter(
      !is.na(cluster),
      !is.na(embryo),
      !is.na(clone_id)
    ) %>%
    mutate(
      # The same barcode in different embryos is a different clone.
      clone_uid = paste(embryo, clone_id, sep = "::")
    ) %>%
    select(cell_id, embryo, clone_id, clone_uid, cluster)

  if (nrow(lineage_cells) == 0) {
    stop("No valid barcode-positive cells remain after filtering.", call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Cluster ordering
  #
  # A factor's level order is meaningful (it usually encodes the biological
  # ordering the authors chose), so it is preserved and handed to the Python
  # step via cluster_summary.csv. Clusters absent from the levels are appended
  # in sorted order.
  # ---------------------------------------------------------------------------
  observed_clusters <- unique(all_clustered_cells$cluster)

  cluster_levels <- NULL
  if (cluster_col != "IDENT" && is.factor(meta[[cluster_col]])) {
    cluster_levels <- levels(meta[[cluster_col]])
  } else if (cluster_col == "IDENT" && is.factor(Idents(obj))) {
    cluster_levels <- levels(Idents(obj))
  }

  if (is.null(cluster_levels)) {
    cluster_order <- sort(observed_clusters)
  } else {
    cluster_order <- c(
      cluster_levels[cluster_levels %in% observed_clusters],
      sort(setdiff(observed_clusters, cluster_levels))
    )
  }

  cluster_order_df <- tibble(
    cluster = cluster_order,
    cluster_order = seq_along(cluster_order)
  )

  # ---------------------------------------------------------------------------
  # Summary tables
  # ---------------------------------------------------------------------------
  cluster_total <- all_clustered_cells %>%
    count(cluster, name = "n_total_cells")

  cluster_barcoded <- lineage_cells %>%
    group_by(cluster) %>%
    summarise(
      n_barcoded_cells = n(),
      n_unique_clones  = n_distinct(clone_uid),
      n_embryos        = n_distinct(embryo),
      .groups = "drop"
    )

  cluster_summary <- cluster_order_df %>%
    left_join(cluster_total, by = "cluster") %>%
    left_join(cluster_barcoded, by = "cluster") %>%
    mutate(
      across(
        c(n_total_cells, n_barcoded_cells, n_unique_clones, n_embryos),
        ~ tidyr::replace_na(.x, 0L)
      ),
      barcode_fraction = if_else(
        n_total_cells > 0,
        n_barcoded_cells / n_total_cells,
        NA_real_
      )
    ) %>%
    arrange(cluster_order)

  clone_summary <- lineage_cells %>%
    group_by(embryo, clone_id, clone_uid) %>%
    summarise(
      clone_size_cells = n(),
      n_clusters = n_distinct(cluster),
      clusters   = paste(sort(unique(cluster)), collapse = "|"),
      .groups = "drop"
    ) %>%
    arrange(embryo, desc(clone_size_cells), clone_uid)

  clone_cluster_counts <- lineage_cells %>%
    count(embryo, clone_id, clone_uid, cluster, name = "n_cells") %>%
    arrange(embryo, clone_uid, cluster)

  embryo_summary <- lineage_cells %>%
    group_by(embryo) %>%
    summarise(
      n_barcoded_cells = n(),
      n_unique_clones  = n_distinct(clone_uid),
      n_clusters       = n_distinct(cluster),
      .groups = "drop"
    ) %>%
    arrange(embryo)

  # ---------------------------------------------------------------------------
  # Write
  # ---------------------------------------------------------------------------
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  tables <- list(
    lineage_cells        = lineage_cells,
    cluster_summary      = cluster_summary,
    clone_summary        = clone_summary,
    clone_cluster_counts = clone_cluster_counts,
    embryo_summary       = embryo_summary
  )

  for (name in names(tables)) {
    write_csv_compat(tables[[name]], file.path(outdir, paste0(name, ".csv")))
  }

  # ---------------------------------------------------------------------------
  # Report
  # ---------------------------------------------------------------------------
  if (verbose) {
    cat("\nExport complete.\n")
    cat("Output directory: ", normalizePath(outdir), "\n", sep = "")
    cat("Total cells with a cluster: ",
        format(nrow(all_clustered_cells), big.mark = ","), "\n", sep = "")
    cat("Barcode-positive cells exported: ",
        format(nrow(lineage_cells), big.mark = ","), "\n", sep = "")
    cat("Unique embryo-specific clones: ",
        format(n_distinct(lineage_cells$clone_uid), big.mark = ","), "\n", sep = "")
    cat("Transcriptomic clusters: ",
        n_distinct(lineage_cells$cluster), "\n", sep = "")
    cat("Embryos: ", n_distinct(lineage_cells$embryo), "\n\n", sep = "")

    print(embryo_summary)
    print(cluster_summary)
  }

  invisible(tables)
}


# =============================================================================
# COMMAND-LINE ENTRY POINT
# =============================================================================

parse_cli <- function(args) {
  opts <- list(
    rds         = NULL,
    cluster_col = NULL,
    clone_col   = NULL,
    embryo_col  = "embryo",
    outdir      = "input"
  )

  key_map <- c(
    "--rds"         = "rds",
    "--cluster-col" = "cluster_col",
    "--clone-col"   = "clone_col",
    "--embryo-col"  = "embryo_col",
    "--outdir"      = "outdir"
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!key %in% names(key_map)) {
      stop("Unknown argument: ", key, call. = FALSE)
    }
    if (i + 1 > length(args)) {
      stop("Argument ", key, " requires a value.", call. = FALSE)
    }
    opts[[key_map[[key]]]] <- args[[i + 1]]
    i <- i + 2
  }

  missing <- names(opts)[vapply(opts, is.null, logical(1))]
  if (length(missing) > 0) {
    stop(
      "Missing required argument(s): --",
      paste(gsub("_", "-", missing), collapse = ", --"),
      call. = FALSE
    )
  }

  opts
}

if (sys.nframe() == 0 && !interactive()) {
  cli <- parse_cli(commandArgs(trailingOnly = TRUE))

  if (!file.exists(cli$rds)) {
    stop("Seurat object not found: ", cli$rds, call. = FALSE)
  }

  message("Reading Seurat object from ", cli$rds, " ...")
  seurat_object <- readRDS(cli$rds)

  export_lineage_input(
    obj         = seurat_object,
    cluster_col = cli$cluster_col,
    clone_col   = cli$clone_col,
    embryo_col  = cli$embryo_col,
    outdir      = cli$outdir
  )
}
