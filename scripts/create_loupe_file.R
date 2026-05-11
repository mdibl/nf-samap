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

    cat("[INFO] Creating Loupe file", "\n")
    create_loupe(
        count_mat   = t(ad$X),
        clusters    = list(
            cell_type = setNames(as.factor(ad$obs$cell_type_labeled), rownames(ad$obs)),
            species   = setNames(as.factor(ad$obs$species),           rownames(ad$obs))
        ),
        projections = list(
            SAMap_UMAP = ad$obsm[['X_umap']]
        ),
        output_name = paste0(args$id1, "_", args$id2, "_Loupe.cloupe")
    )
    cat("[INFO] Done. Loupe file written to output directory", "\n")
}

# --------------------------------------------------
main()