#!/usr/bin/env Rscript

# --------------------------------------------------
# Libraries
suppressPackageStartupMessages({
  library(argparse)
  library(readr)
  library(Seurat)
  library(Matrix)
})

# --------------------------------------------------
# Define and parse arguments (analogous to Python get_args())
get_args <- function() {
  parser <- ArgumentParser(
    description = "Build a SAMAP object from a directory of SAMs and a sample sheet",
    formatter_class = "argparse.ArgumentDefaultsHelpFormatter"
  )

  parser$add_argument(
    "-s", "--so",
    required = TRUE,
    help = "Path to Seurat Object"
  )

  parser$add_argument(
    "-i", "--id",
    required = TRUE,
    help = "ID for the Seurat Object"
  )

  parser$add_argument(
    "-a", "--anno",
    required = TRUE,
    help = "Annotation layer to be extracted from the Seurat Object"
  )


  args <- parser$parse_args()

  # Create a structured object similar to a Python NamedTuple
  Args <- list(
    so = args$so,
    id = args$id,
    anno = args$anno
  )

  return(Args)
}

# MAIN -------------------------------------------------- MAIN

if (sys.nframe() == 0) {
  args <- get_args()
}
if (!file.exists(args$so)) {
    stop(paste0("Seurat Object not found for ID: ", args$id))
}
seurat_obj <- readRDS(args$so)
print(length(seurat_obj@meta.data[["old.ident"]]))
print(length(seurat_obj@meta.data[["orig.ident"]]))
print(length(Cells(seurat_obj)))
print(length(seurat_obj@meta.data[[Anno]]))


ID <- args$id
Anno <- args$anno

DefaultAssay(seurat_obj) <- "RNA"

#".x" part of the AnnData Object, swapped rows and columns to fit correct dimensions
Counts <- as(t(seurat_obj@assays[["RNA"]]@layers[["counts"]]), "dgCMatrix")

#".var" part of the AnnData Object, maybe make into df? Look at desired format
Feats <- Features(seurat_obj)

##Cell Metadata, including original sample identity and Harmony clustering label (resolution of choice)
obs <- data.frame(
  cell_id = Cells(seurat_obj),
  sample  = seurat_obj@meta.data[["old.ident"]],
  value   = seurat_obj@meta.data[[Anno]]
)

# Rename the column dynamically
colnames(obs)[which(names(obs) == "value")] <- Anno


Counts <- as(Counts, "CsparseMatrix")
writeMM(Counts, file = paste0(ID,"_Counts.mtx"))

write.csv(obs, paste0(ID,"_Obs.csv"), row.names = FALSE)

write.table(data.frame(values = unlist(Feats)),
            paste0(ID,"_Feats.csv"),
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = ",")