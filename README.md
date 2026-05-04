# nf-samap: Cross-species transcriptome alignment pipeline

A [Nextflow](https://www.nextflow.io/) pipeline wrapping [SAMap](https://github.com/atarashansky/SAMap) to perform cross-species single-cell RNA-seq alignment. Given expression data (Seurat objects) and reference transcriptomes/proteomes for two or more species, the pipeline runs reciprocal BLAST to build gene-homology maps and then uses SAMap to compute gene-level mapping scores and pairwise differential expression analysis across species.

---

## Prerequisites

- [Nextflow](https://www.nextflow.io/docs/latest/install.html) (≥ 22.10)
- [Docker](https://docs.docker.com/get-docker/)
- The pipeline's Docker images (see step 1 below)

---

## Setup

### 1. Build Docker images

The pipeline uses four custom images. Build and push all of them with:

```bash
make docker
```

Or build them individually:

```bash
docker build --platform=linux/amd64 -f Dockerfile.preprocessing -t docker.io/mdiblbiocore/preprocessing:latest .
docker build -f Dockerfile.samap -t docker.io/mdiblbiocore/samap:latest .
docker build --platform=linux/amd64 -f Dockerfile.blast -t docker.io/mdiblbiocore/samap-blast:latest .
docker build --platform=linux/amd64 -f Dockerfile.postanalysis -t docker.io/mdiblbiocore/postanalysis:latest .
```

> If you are pulling from an existing registry instead of building locally, you can skip this step as long as Docker can reach the `mdiblbiocore` images listed in `nextflow.config`.

### 2. Prepare your input data

You need three things per species:

| File | Description |
|------|-------------|
| Seurat object (`.RDS`) | Single-cell expression data. Must contain a metadata column with cell-type annotations. |
| Transcriptome or proteome FASTA (`.fasta`) | Reference sequences used for reciprocal BLAST. Can be nucleotide (transcriptome) or amino acid (proteome). |
| A row in `sample_sheet.csv` | Metadata linking the above files and providing pipeline configuration for each species. |

A typical directory layout:

```
sample_sheet.csv
data/
├── transcriptomes/
│   ├── planarian.fasta
│   └── hydra.fasta
├── planarian.RDS
└── hydra.RDS
```

### 3. Create your sample sheet

The sample sheet is a CSV that describes each species. Each row is one species.

```csv
id,so,fasta,annotation,mapping_dict
00,data/planarian.RDS,data/transcriptomes/planarian.fasta,cluster,
01,data/hydra.RDS,data/transcriptomes/hydra.fasta,Cluster,
```

#### Required columns

| Column | Description |
|--------|-------------|
| `id` | A unique **two-character** alphanumeric identifier for the species (e.g. `00`, `pl`, `hy`). Used to prefix all output files for that species. |
| `so` | Path to the Seurat object (`.RDS`). The pipeline extracts the counts matrix, cell metadata, and gene features from this file. |
| `fasta` | Path to the transcriptome or proteome FASTA file. Used as the database for reciprocal BLAST between species pairs. Nucleotide and amino acid FASTAs can be mixed across species — the pipeline selects the correct BLAST program automatically (`tblastx`, `blastx`, `tblastn`, or `blastp`). |
| `annotation` | The name of the metadata column in the Seurat object that contains cell-type labels (e.g. `cluster`, `cell_type`, `tissue`). This is used for grouping cells in post-analysis. |

#### Optional columns

| Column | Description | When to use |
|--------|-------------|-------------|
| `mapping_dict` | Path to a two-column text file mapping FASTA sequence IDs to the gene symbols used in the Seurat object. | Required when the identifiers in your FASTA file (e.g. transcript IDs like `TRINITY_DN1234_c0_g1_i1`) differ from the gene names in your Seurat object (e.g. `slit-1`). Each line should contain one FASTA ID and one gene symbol, comma-separated. |
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

---

## Output structure

All outputs land under `--outdir` (default: `out/`):

```
out/
├── preprocess_seurat_object/       # Intermediate: counts matrix, obs, and feature tables
├── preprocess_anndata_object/      # Intermediate: .h5ad files built from Seurat objects
├── blast/                          # Reciprocal BLAST maps (one subdirectory per species pair)
├── load_sams/                      # Pickled SAM objects
├── build_samap/                    # Pickled SAMAP object (pre-algorithm)
├── run_samap/                      # Pickled SAMAP results (post-algorithm)
└── Analysis/
    └── AnalysisResults/
        └── {id1}-{id2}/            # One directory per species pair
            ├── GenePairs.csv       # Top gene-pair mappings with scores
            ├── pms_cluster_alignment_scores.csv  # Pairwise cell-cluster alignment scores
            └── Grouping_Analysis/  # Per-alignment-family differential expression tables
```

> **Note:** Visualization outputs (Sankey diagrams, scatter plots) are not currently produced — that module is disabled pending improvements. Summary statistics and differential expression results are the primary analysis outputs.

---

## Pipeline stages

| Stage | Module | Description |
|-------|--------|-------------|
| 1 | `PREPROCESS_SEURAT_OBJECT` | Extracts counts matrix, cell metadata (obs), and gene metadata (feats) from each Seurat `.RDS` file using R. |
| 2 | `PREPROCESS_ANNDATA_OBJECT` | Converts the extracted data into an `.h5ad` (AnnData) file for each species. |
| 3 | `RUN_BLAST_PAIR` | For every unordered species pair, runs reciprocal BLAST using the appropriate program for the FASTA types. Skipped if `--maps_dir` is set. |
| 4 | `LOAD_SAMS` | Loads each `.h5ad` into a SAM (Self-Assembling Manifold) object. |
| 5 | `BUILD_SAMAP` | Combines all SAM objects with the BLAST maps to construct a SAMAP object. |
| 6 | `RUN_SAMAP` | Runs the SAMap algorithm to compute gene-to-gene mapping scores across species. |
| 7 | `SUMMARY_SAMAP` | For each species pair: extracts top gene-pair mappings and computes pairwise cell-cluster alignment scores. |
| 8 | `CONNECTED_DE` | Identifies alignment families and runs differential expression analysis within each family across both species. |
| 9 | `ADDITIONAL_ANALYSIS` | Integrates gene-pair scores with DE results into final per-pair output tables. |

---

## Troubleshooting

**Pipeline fails at BLAST with "command not found"**  
Make sure Docker is running and the `mdiblbiocore/samap-blast` image is available locally. Run `make docker` or pull the image manually.

**Gene names not matching between FASTA and Seurat object**  
If your FASTA uses transcript IDs (e.g. Trinity IDs) while your Seurat object uses gene symbols, provide a `mapping_dict` file for the affected species in the sample sheet.

**Out-of-memory errors**  
BLAST and SAMap are memory-intensive. Use `-profile worm` for larger datasets, or override resources in your params file and pass `-profile docker`.

**Resuming after a partial run**  
Always add `-resume` to avoid re-running completed steps:

```bash
nextflow run main.nf -params-file params.json -resume
```

**Inspect logs**  
Per-module log files are written alongside their outputs (e.g. `out/run_samap/*.log`).

---

## Links and acknowledgements

- [SAMap repository](https://github.com/atarashansky/SAMap)
- [SAMap paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8139856/)
- [SAMap Docker image](https://hub.docker.com/r/avianalter/samap)
- [BLAST Docker image](https://hub.docker.com/r/staphb/blast)

---

## Authors

**Markus Sujansky** (current maintainer) — [MDIBL Bioinformatics Core](https://mdibl.org/)

**Ryan Sonderman** (original author) — [@RyanSonder](https://github.com/RyanSonder)

**Riley Grindle** — [@Riley-Grindle](https://github.com/Riley-Grindle)

This pipeline is licensed under the MIT License. See [LICENSE](LICENSE) for details.
