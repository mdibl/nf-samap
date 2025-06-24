# 🧬 SAMap Pipeline - Cross-species transcriptome alignment using BLAST and SAMap

> This project wraps SAMap (https://github.com/atarashansky/SAMap) in a NextFlow pipeline. 

> Current version: `v1.1.0`

---
# 🚀 Quickstart

This project uses Makefile to simplify many of the necessary actions. Each step can be done manually or with a make target. 

I will soon add a make target to clone the example data from SAMap and format it correctly.

## 1. Build the custom docker images
```bash
docker build -f Dockerfile.samap -t pipeline/samap:latest .
docker build -f Dockerfile.blast -t pipeline/samap-blast:latest .
```
-or-
```bash
make docker
```

## 2. Run the pipeline
```bash
nextflow run main.nf --with-docker
```
-or-
```bash
make run
```

---
# 📂 Input Files

The pipeline expects the following input files to be present:

```
sample_sheet.csv
*.fasta
*.h5ad
```

An example tree:
```
sample_sheet.csv
data/
├── transcriptomes/
│   ├── hydra.fasta
│   ├── planarian.fasta
│   └── schistosome.fasta
├── hydra.h5ad
├── planarian.h5ad
└── schistosome.h5ad
```

---
# 📜 Sample Sheet Format

The sample sheet dictates metadata about each sample. Samples will not be put through the pipeline unless they are present and correctly described in the sample sheet. An example `sample_sheet.csv` might look like:

```csv
id,h5ad,fasta,annotation
00,data/planarian.h5ad,data/transcriptomes/planarian_transcriptome.fasta,cluster
01,data/hydra_mod.h5ad,data/transcriptomes/hydra_transcriptome.fasta,Cluster
02,data/schistosome.h5ad,data/transcriptomes/schistosome_proteome.fasta,tissue
```

---
# ⚙️ Parameters

| Parameter | Requirement | Description | Default |
|-----------|-------------|-------------|---------|
| `run_id` | Optional | Custom run ID | `null` |
| `sample_sheet` | Optional | Sample sheet describing sample metadata | `'sample_sheet.csv'` |
| `data_dir` | Optional | Path to directory containing sample data | `'data'` |
| `maps_dir` | Optional | Path to directory of precomputed BLAST maps | `null` |
| `results_dir` | Optional | Path to directory where results are stored | `'results'` |

---
# 🏁 Output Files

Results are stored in `results/{run_id}/`.

| Path | Description |
|------|-------------|
| {run_id}_sample_sheet.csv | Processed sample sheet |
| csv/hms.csv | Highest mapping scores |
| csv/pms.csv | Pairwise mapping scores |
| plots/chord.html | Chord plot |
| plots/sankey.html | Sankey plot |
| plots/scatter.png | Scatterplot |
| samap_objects/samap_results.pkl | Pickled SAMAP object after running SAMap |
| samap_objects/samap.pkl | Pickled SAMAP object before running SAMap |
| sams/* | Pickled SAM objects named according to the 2-char hash assigned to their sample |
| logs/* | Logfile output for each module |

---
# 🧱 Module Overview

### 1. PREPROCESS

Reads the sample_sheet.csv, classifies transcriptomes based on input FASTA files, and assigns unique two-character IDs. Outputs an enriched sample sheet with metadata used downstream.

### 2. RUN_BLAST_PAIR

For each unique unordered species pair, performs a reciprocal BLAST to generate mapping files. Skipped if --use_precomputed_blast is true.

### 3. LOAD_SAMS

Loads input .h5ad files and constructs SAM objects required for SAMap. Outputs pickled SAM objects.

### 4. BUILD_SAMAP

Combines the SAM objects and reciprocal BLAST maps to build a SAMAP object.

### 5. RUN_SAMAP

Runs the SAMap algorithm on the built object to calculate pairwise gene mapping scores.

### 6. VISUALIZE_SAMAP

Generates outputs such as Sankey diagrams, scatter plots, and CSV summaries of the alignment results for downstream analysis or interpretation.

---
# 🔗 Links and Acknowledgements

- SAMap Repository:   https://github.com/atarashansky/SAMap
- SAMap Paper:        https://pmc.ncbi.nlm.nih.gov/articles/PMC8139856/
- SAMap Docker Image: https://hub.docker.com/r/avianalter/samap
- BLAST Docker Image: https://hub.docker.com/r/staphb/blast

---
# 👤 Authors and Licenses

**Ryan Sonderman**

GitHub: [@RyanSonder](https://github.com/RyanSonder)

**Riley Grindle**

GitHub: [@Riley-Grindle](https://github.com/Riley-Grindle)

This pipeline is licensed under the MIT License. See the [LICENSE](LICENSE) file for full details.

---
# 📋 To Do List

- [ ] **Visualization Improvements**
    - Sort the Sankey diagram for better interpretability
    - Add a legend for input organisms in the Sankey and scatter plots
    - Improve coloring in scatter plots for distinct organism/group visualization

- [ ] **Reproducibility & Reporting**
    - Add version reporting for SAMap in logs and outputs
    - Ensure Docker image versioning is clear and consistent (`pipeline/samap:v1.0.0`)

- [ ] **Data Accesibility**
    - Add a way to easily clone the data from the original SAMap repo

- [ ] **nf-core requirements**
    - [x] Acknowledgements
    - [ ] Continuous integration testing
    - [x] Community owned
    - [ ] Docker support (no latest)
    - [ ] Bundled documentation
    - [ ] Use nf-core git branches
    - [ ] Identity and branding
    - [ ] GitHub keywords
    - [ ] Pass lint tests
    - [ ] Minimum inputs
    - [x] MIT License
    - [x] Nextflow
    - [ ] Standardised parameters
    - [ ] Research object crate
    - [ ] Semantic versioning
    - [x] Single command 
    - [ ] Use the template
    - [ ] Workflow name
    - [x] Workflow size

- [ ] **nf-core recommendations**
    - [ ] Publication credit
    - [ ] Testing
    - [ ] Bioconda
    - [ ] Build with community
    - [ ] Cloud compatible
    - [ ] Custom containers
    - [ ] DOIs
    - [ ] Fiel formats
