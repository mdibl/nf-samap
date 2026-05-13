#!/usr/bin/env Rscript

# Author : Markus Sujansky
# Date   : 2026-05-07
# Version: 2.0.6
# Purpose: Build Shared Feature Space Loupe File

# --------------------------------------------------
# Libraries
suppressPackageStartupMessages({
    library(argparse)
    library(loupeR)
    library(Matrix)
})

# --------------------------------------------------
get_args <- function() {
    parser <- ArgumentParser(
        description = "Build a shared feature space Loupe file from SAMap output",
        formatter_class = "argparse.ArgumentDefaultsHelpFormatter"
    )

    parser$add_argument("--id1",      required=TRUE, help="First species label")
    parser$add_argument("--id2",      required=TRUE, help="Second species label")
    parser$add_argument("--counts",   required=TRUE, help="Path to counts .mtx file")
    parser$add_argument("--barcodes", required=TRUE, help="Path to barcodes .csv file")
    parser$add_argument("--features", required=TRUE, help="Path to features .csv file")
    parser$add_argument("--umap",     required=TRUE, help="Path to UMAP .csv file")
    parser$add_argument("--meta",     required=TRUE, help="Path to metadata .csv file")

    args <- parser$parse_args()
    return(args)
}

# --------------------------------------------------
main <- function() {
    args <- get_args()
    prefix <- paste0(args$id1, "_", args$id2)

    tryCatch(
        loupeR::setup(),
        error = function(e) cat("[INFO] loupeR already set up\n")
    )

    cat("[INFO] Loading count matrix\n")
    counts   <- readMM(args$counts)
    barcodes <- read.csv(args$barcodes, header=FALSE)[[1]]
    features <- read.csv(args$features, header=FALSE)[[1]]
    rownames(counts) <- features
    colnames(counts) <- barcodes
    counts <- as(counts, "dgCMatrix")
    cat("[INFO] Count matrix dims:", dim(counts), "\n")

    cat("[INFO] Loading metadata\n")
    meta <- read.csv(args$meta, row.names=1)

    cat("[INFO] Loading UMAP\n")
    umap <- read.csv(args$umap, row.names=1)

    cat("[INFO] Creating Loupe file\n")
    create_loupe(
        count_mat   = counts,
        clusters    = list(
            cell_type        = setNames(as.factor(meta$cell_type_labeled), rownames(meta)),
            species          = setNames(as.factor(meta$species),           rownames(meta)),
            alignment_family = setNames(as.factor(meta$alignment_family),  rownames(meta))
        ),
        projections = list(
            SAMap_UMAP = as.matrix(umap)
        ),
        output_name = prefix
    )
    cat("[INFO] Done. Loupe file written to", paste0(prefix, ".cloupe"), "\n")
}

# --------------------------------------------------
main()