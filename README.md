# nf-samap: Cross-species transcriptome alignment pipeline

A [Nextflow](https://www.nextflow.io/) pipeline wrapping [SAMap](https://github.com/atarashansky/SAMap) to perform cross-species single-cell RNA-seq alignment. Given expression data (Seurat objects or AnnData files) and reference transcriptomes/proteomes for two or more species, the pipeline runs reciprocal BLAST to build gene-homology maps and then uses SAMap to compute gene-level mapping scores, pairwise differential expression analysis, and optional [10x Loupe Browser](https://www.10xgenomics.com/support/software/loupe-browser) visualization files.

---

## Prerequisites

- [Nextflow](https://www.nextflow.io/docs/latest/install.html) (≥ 22.10)
- [Docker](https://docs.docker.com/get-docker/)
- The pipeline's Docker images (see step 1 below)

---

## Setup

### 1. Build Docker images

The pipeline uses five custom images. Build and push all of them with:

```bash
make docker
```

Or build them individually:

```bash
docker build --platform=linux/amd64 -f Dockerfile.preprocessing -t docker.io/mdiblbiocore/preprocessing:latest .
docker build -f Dockerfile.samap -t docker.io/mdiblbiocore/samap:latest .
docker build --platform=linux/amd64 -f Dockerfile.blast -t docker.io/mdiblbiocore/samap-blast:latest .
docker build --platform=linux/amd64 -f Dockerfile.postanalysis -t docker.io/mdiblbiocore/postanalysis:latest .
docker build --platform=linux/amd64 -f Dockerfile.loupe -t docker.io/mdiblbiocore/loupe:latest .
```

> If you are pulling from an existing registry instead of building locally, you can skip this step as long as Docker can reach the `mdiblbiocore` images listed in `nextflow.config`.

### 2. Prepare your input data

You need three things per species:

| File | Description |
|------|-------------|
| Expression data (`.RDS` or `.h5ad`) | Single-cell expression data. Seurat objects (`.rds`, `.rdata`, `.rda`) and AnnData files (`.h5ad`, `.h5`, `.loom`) are both accepted. Must contain a metadata column with cell-type annotations. |
| Transcriptome or proteome FASTA (`.fasta`) | Reference sequences used for reciprocal BLAST. Can be nucleotide (transcriptome) or amino acid (proteome). |
| A row in `sample_sheet.csv` | Metadata linking the above files and providing pipeline configuration for each species. |

A typical directory layout:

```
sample_sheet.csv
data/
├── transcriptomes/
│   ├── planarian.fasta
│   └── hydra.fasta
├── planarian.h5ad
└── hydra.h5ad
```

### 3. Create your sample sheet

The sample sheet is a CSV that describes each species. Each row is one species.

```csv
id,so,fasta,annotation,mapping_dict
00,data/planarian.h5ad,data/transcriptomes/planarian.fasta,cluster,
01,data/hydra.h5ad,data/transcriptomes/hydra.fasta,Cluster,
```

#### Required columns

| Column | Description |
|--------|-------------|
| `id` | A unique **two-character** alphanumeric identifier for the species (e.g. `00`, `pl`, `hy`). Used to prefix all output files for that species. |
| `so` | Path to the expression data file. Accepts Seurat objects (`.rds`/`.rdata`/`.rda`) or AnnData files (`.h5ad`/`.h5`/`.loom`). The pipeline detects the format automatically and routes through the appropriate preprocessing steps. |
| `annotation` | The name of the metadata column in the expression object that contains cell-type labels (e.g. `cluster`, `cell_type`, `tissue`). This is used for grouping cells in post-analysis. |

#### Optional columns

| Column | Description | When to use |
|--------|-------------|-------------|
| `fasta` | Path to the transcriptome or proteome FASTA file. Used as the database for reciprocal BLAST between species pairs. Nucleotide and amino acid FASTAs can be mixed across species — the pipeline selects the correct BLAST program automatically (`tblastx`, `blastx`, `tblastn`, or `blastp`). | Required unless you are providing precomputed BLAST maps via `maps_dir`, in which case this column can be omitted entirely. |
| `mapping_dict` | Path to a two-column text file mapping FASTA sequence IDs to the gene symbols used in the expression object. | Required when the identifiers in your FASTA file (e.g. transcript IDs like `TRINITY_DN1234_c0_g1_i1`) differ from the gene names in your data (e.g. `slit-1`). Each line should contain one FASTA ID and one gene symbol, comma-separated. |
| `type` | Sequence type — `nucl` for nucleotide transcriptomes, `prot` for amino acid proteomes. | Optional; the pipeline can infer this automatically in most cases, but providing it explicitly avoids ambiguity when file extensions are non-standard. |

### 4. Create a params file

Pipeline parameters should be provided via a JSON params file rather than on the command line. Create a file (e.g. `params.json`) in your working directory:

```json
{
    "sample_sheet": "sample_sheet.csv",
    "outdir": "out",
    "run_id": "my_experiment"
}
```

Only include parameters you want to override — anything omitted falls back to the defaults in `nextflow.config`.

### 5. Run the pipeline

```bash
nextflow run main.nf -params-file params.json
```

To resume a failed run without re-running completed steps:

```bash
nextflow run main.nf -params-file params.json -resume
```

---

## Parameters

These go in your `params.json` and are passed to Nextflow with `-params-file params.json`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `sample_sheet` | `sample_sheet.csv` | Path to the sample sheet CSV. |
| `outdir` | `out` | Directory where all outputs are written. |
| `run_id` | *(timestamp)* | Label for this run, used to prefix log files. If not provided, a `yyyyMMdd_HHmmss` timestamp is used. |
| `maps_dir` | `null` | Path to a directory of **precomputed BLAST maps**. If provided, the BLAST step is skipped entirely and these maps are used directly. Useful for re-running analysis with different parameters without repeating the (slow) BLAST step. |
| `var_genes` | `3000` | Number of highly variable genes selected per species during SAM preprocessing. |
| `create_pairwise_loupe` | `'true'` | Generate a `.cloupe` file for each species pair, viewable in 10x Loupe Browser. Requires agreeing to the 10x EULA (see `eula10x`). |
| `create_all_loupe` | `'false'` | Generate a single `.cloupe` file combining all species. Requires agreeing to the 10x EULA (see `eula10x`). |
| `eula10x` | `'Agree'` | Must be set to `"Agree"` to enable Loupe file generation. Set to any other value to disable Loupe output and skip the EULA gate. |

---

## Compute profiles

Pass `-profile <name>` to select a resource profile:

| Profile | Description |
|---------|-------------|
| `docker` | Default. Enables Docker with per-process CPU/memory allocations tuned for a workstation or small server (64 CPUs for BLAST, 32 for SAMap). |
| `worm` | Adjusted allocations for large multi-species datasets. |
| `cluster` | Conservative settings (8 CPUs, 8 GB) for HPC login/testing nodes. |
| `test` | Minimal resources (1 CPU, 8 GB, 10-minute timeout) for quick smoke tests. |
| `arm` | Adds `--platform=linux/amd64` to Docker run options for Apple Silicon Macs. |

Example:

```bash
nextflow run main.nf -profile docker
```

---

## Using precomputed BLAST maps

BLAST is the slowest step. If you have already run the pipeline once (or have maps from another source), you can skip it by setting `maps_dir` in your params file:

```json
{
    "sample_sheet": "sample_sheet.csv",
    "outdir": "out",
    "maps_dir": "out/blast/maps/"
}
```

```bash
nextflow run main.nf -params-file params.json
```

The maps directory should contain the reciprocal BLAST output files generated by `scripts/map_genes.sh` (one file per species pair, in both directions).

When using precomputed maps, the `fasta` column is not required in your sample sheet and can be omitted.

---

## Output structure

All outputs land under `--outdir` (default: `out/`):

```
out/
├── preprocess_seurat_object/       # Intermediate: counts matrix, obs, and feature tables (Seurat input only)
├── preprocess_anndata_object/      # Intermediate: initialized .h5ad files (Seurat input only)
├── preprocess_sam_object/          # Intermediate: preprocessed .h5ad files with variable genes selected
├── blast/                          # Reciprocal BLAST maps (one subdirectory per species pair)
├── load_sams/                      # Pickled SAM objects
├── build_samap/                    # Pickled SAMAP object (pre-algorithm)
├── run_samap/                      # Pickled SAMAP results (post-algorithm)
└── Analysis/
    ├── AllSpecies/
    │   └── Loupe/                  # All-species .cloupe file (if create_all_loupe = true)
    └── {id1}-{id2}/                # One directory per species pair
        ├── AnalysisResults/        # GenePairs.csv, pms_cluster_alignment_scores.csv, Grouping_Analysis/
        ├── AlignmentFamilies/      # Per-family alignment confusion tables
        └── Loupe/                  # Pairwise .cloupe file (if create_pairwise_loupe = true)
```

---

## Pipeline stages

| Stage | Module | Description |
|-------|--------|-------------|
| 1 | `PREPROCESS_SEURAT_OBJECT` | Extracts counts matrix, cell metadata (obs), and gene metadata (feats) from each Seurat `.RDS` file using R. Skipped if input is already an AnnData file. |
| 2 | `PREPROCESS_ANNDATA_OBJECT` | Converts the extracted Seurat data into an initialized `.h5ad` (AnnData) file for each species. Skipped if input is already an AnnData file. |
| 3 | `PREPROCESS_SAM_OBJECT` | Preprocesses the AnnData object for SAMap: filters, normalizes, and selects the top `var_genes` highly variable genes. Produces a `_preprocessed.h5ad` used downstream. |
| 4 | `RUN_BLAST_PAIR` | For every unordered species pair, runs reciprocal BLAST using the appropriate program for the FASTA types. Skipped if `--maps_dir` is set. |
| 5 | `LOAD_SAMS` | Loads each preprocessed `.h5ad` into a SAM (Self-Assembling Manifold) object. |
| 6 | `BUILD_SAMAP` | Combines all SAM objects with the BLAST maps to construct a SAMAP object. |
| 7 | `RUN_SAMAP` | Runs the SAMap algorithm to compute gene-to-gene mapping scores across species. |
| 8 | `SUMMARY_SAMAP` | For each species pair: extracts top gene-pair mappings, computes pairwise cell-cluster alignment scores, and generates summary visualizations. |
| 9 | `CONNECTED_DE` | Identifies alignment families and runs differential expression analysis within each family across both species. |
| 10 | `ADDITIONAL_ANALYSIS` | Integrates gene-pair scores with DE results into final per-pair output tables. |
| 11 | `CREATE_LOUPE_INPUT` | Merges raw counts matrices and SAMap UMAP embeddings for each species pair into Loupe-ready inputs. Runs if `create_pairwise_loupe = 'true'`. |
| 12 | `CREATE_LOUPE_FILE` | Builds a pairwise `.cloupe` file from the merged counts and UMAP. Runs if `create_pairwise_loupe = 'true'`. |
| 13 | `CREATE_ALL_LOUPE_INPUT` | Merges raw counts and UMAP embeddings for all species into Loupe-ready inputs. Runs if `create_all_loupe = 'true'`. |
| 14 | `CREATE_ALL_LOUPE_FILE` | Builds a single all-species `.cloupe` file. Runs if `create_all_loupe = 'true'`. |

---

## Loupe Browser output

When `create_pairwise_loupe = 'true'` (the default), the pipeline generates a `.cloupe` file for each species pair. These can be opened in [10x Loupe Browser](https://www.10xgenomics.com/support/software/loupe-browser) and contain:

- SAMap UMAP embedding with cells from both species
- Cell-type labels per species
- Alignment family assignments per cell

An all-species Loupe file can also be generated by setting `create_all_loupe = 'true'`.

**EULA:** Loupe file creation uses the 10x `loupeR` package, which requires accepting the [10x Genomics End-User License Agreement](https://www.10xgenomics.com/legal/end-user-software-license-agreement). The pipeline enforces this via the `eula10x` parameter — set it to any value other than `"Agree"` to disable Loupe output entirely.

---

## Troubleshooting

**Pipeline fails at BLAST with "command not found"**  
Make sure Docker is running and the `mdiblbiocore/samap-blast` image is available locally. Run `make docker` or pull the image manually.

**Gene names not matching between FASTA and expression object**  
If your FASTA uses transcript IDs (e.g. Trinity IDs) while your data uses gene symbols, provide a `mapping_dict` file for the affected species in the sample sheet.

**Out-of-memory errors**  
BLAST and SAMap are memory-intensive. Use `-profile worm` for larger datasets, or override resources in your params file and pass `-profile docker`.

**Resuming after a partial run**  
Always add `-resume` to avoid re-running completed steps:

```bash
nextflow run main.nf -params-file params.json -resume
```

**Inspect logs**  
Per-module log files are written alongside their outputs (e.g. `out/run_samap/*.log`).

**Loupe file generation fails with EULA error**  
The `eula10x` parameter must be set to `"Agree"` in your params file. If you do not want to generate Loupe files, set `create_pairwise_loupe = 'false'` and `create_all_loupe = 'false'`.

---

## Links and acknowledgements

- [SAMap repository](https://github.com/atarashansky/SAMap)
- [SAMap paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8139856/)
- [SAMap Docker image](https://hub.docker.com/r/avianalter/samap)
- [BLAST Docker image](https://hub.docker.com/r/staphb/blast)
- [10x Loupe Browser](https://www.10xgenomics.com/support/software/loupe-browser)

---

## Authors

**Markus Sujansky** (current maintainer) — [MDIBL Bioinformatics Core](https://mdibl.org/)

**Ryan Sonderman** (original author) — [@RyanSonder](https://github.com/RyanSonder)

**Riley Grindle** — [@Riley-Grindle](https://github.com/Riley-Grindle)

This pipeline is licensed under the MIT License. See [LICENSE](LICENSE) for details.
