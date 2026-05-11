#!/usr/bin/env Rscript

# Author : Markus Sujansky
# Date   : 2026-05-07
# Version: 2.0.3
# Purpose: Build Shared Feature Space Loupe File

# --------------------------------------------------
# Libraries
suppressPackageStartupMessages({
    library(argparse)
    library(loupeR)
    library(anndata)
})

# --------------------------------------------------
# Define and parse arguments
get_args <- function() {
    parser <- ArgumentParser(
        description = "Build a shared feature space Loupe file from SAMap output",
        formatter_class = "argparse.ArgumentDefaultsHelpFormatter"
    )
    
    parser$add_argument(
        "--id1",
        required = TRUE,
        help = "First species label"
    )
    parser$add_argument(
        "--id2",
        required = TRUE,
        help = "Second species label"
    )
    parser$add_argument(
        "--h5ad",
        required = TRUE,
        help = "Path to the combined h5ad file produced by create_loupe_input.py"
    )

    args <- parser$parse_args()
    return(args)
}

# --------------------------------------------------

main <- function() {
    args <- get_args()

    tryCatch(
        loupeR::setup(),
        error = function(e) cat("[INFO] loupeR already set up\n")
    )

    cat("[INFO] Loading combined h5ad from", args$h5ad, "\n")
    ad <- read_h5ad(args$h5ad)
    cat("[INFO] ad$X class:", class(ad$X), "\n")
    cat("[INFO] ad$X type:", typeof(ad$X), "\n")

    cat("[INFO] Creating Loupe file", "\n")
    barcodes <- if (!is.null(rownames(ad$obs))) rownames(ad$obs) else ad$obs_names

    count_mat <- as(ad$X, "dgCMatrix")


    create_loupe(
        count_mat   = t(count_mat),
        clusters    = list(
            cell_type = setNames(as.factor(ad$obs$cell_type_labeled), barcodes),
            species   = setNames(as.factor(ad$obs$species),           barcodes)
        ),
        projections = list(
            SAMap_UMAP = as.matrix(ad$obsm[['X_umap']])
        ),
        output_name = paste0(args$id1, "_", args$id2, "_Loupe")
    )
    cat("[INFO] Done. Loupe file written to output directory", "\n")
}

# --------------------------------------------------
main()