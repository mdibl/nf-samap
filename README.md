# Ryan Sonderman's SAMap Pipeline

## Source

This pipeline is based on SAMap, which can be found at https://github.com/atarashansky/SAMap.

## Citation

Tarashansky, Alexander J., et al. "Mapping single-cell atlases throughout Metazoa unravels cell type evolution." Elife 10 (2021): e66747.

## About

This project uses a Nextflow pipeline to run SAMap, using a modified docker image from `docker.io/avianalter/samap`. 

Currently, the pipeline is designed to take a JSON config describing the locations of your input files and orthology maps. The expected directory structure is:

```
data/
├── hydro.h5ad
├── planarian.h5ad
├── schistosome.h5ad
└── maps
    ├── hypl
    │   ├── hy_to_pl.txt
    │   └── pl_to_hy.txt
    ├── hysc
    │   ├── hy_to_sc.txt
    │   └── sc_to_hy.txt
    └── plsc
        ├── pl_to_sc.txt
        └── sc_to_pl.txt
```

## Usage

### 1. Build the Docker Image

```bash
make docker
```

This will update the docker image with any changes to scripts, patches, etc.

### 2. Prepare Your Config

Create a `config.json` file describing your species and map locations. Example that follows the previously mentioned tree:

```json
{
  "maps": "data/maps/",
  "species": {
    "pl": "data/planarian.h5ad",
    "hy": "data/hydra.h5ad",
    "sc": "data/schistosome.h5ad"
  }
}
```

### 3. Run the Pipeline

```bash
make run
```

This will run the Nextflow pipeline using Docker.

## Development

- To open a shell inside the Docker container with your workspace mounted:

  ```bash
  make docker-shell
  ```

## Output

The pipeline will generate:

- `samap_obj.pkl`: The main SAMap results object (pickle)
- more visualization stuff (to be implemented)

## Notes

- All output files are written with the correct user permissions when using the provided Makefile or Nextflow pipeline.
- If you modify the scripts, rebuild the Docker image with `make docker`.
- For custom species or obs keys, provide a JSON mapping as described in the documentation.

## License

See [LICENSE](LICENSE) for details.

---


## TODO

- [ ] **Helper Scripts**
    - Bash script to build the keys JSON for the Python visualization script
    - Add support for SAM format I/O (SAM → JSON → SAMap)

- [ ] **Visualization Improvements**
    - Sort the Sankey diagram for better interpretability
    - Add a legend for input organisms in the Sankey and scatter plots
    - Improve coloring in scatter plots for distinct organism/group visualization

- [ ] **Reproducibility & Reporting**
    - Add version reporting for SAMap in logs and outputs
    - Ensure Docker image versioning is clear and consistent (`ryansonder/samap:v1.0.0`)