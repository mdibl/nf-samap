# my-samap

## Instructions from Joel:

Here is a rough overview:

    1. Acquire SAMap software.  Run the test examples to make sure we get the expected, check for computational requirements, etc
    2. Ideally (but probably optional), wrap the workflow that was presented as a jupyter notebook in the github repo into a nextflow workflow instead.  (Riley and Ryan both have experience with this, and are your best resources for how to proceed)
    3. Run the axolotl mutliome data through the multiome pipeline on 10X resources.  [alternative, or potential additional task: run the same data through the mutliome option of the nf-core/scrnaseq pipeline]
    4. Run the mouse multiome data through the same pipeline
    5. Extract only the expression data from the multiome sets and convert them into the hd5 format expected by SAMAP
    6. Run the SAMap pipeline on the combined mouse and axolotl data- and interpret and visualize the output 

Assuming success on what was written above, we will likely want to follow up by making SAMap run on the data that Talia is working on, comparative retina samples from zebrafish, mouse, and chicken.


TODO
redesign the python inputs to work with a data/ path assuming the following example tree structure:
data/
├── hydra.h5ad
├── maps
│   ├── hypl
│   │   ├── hy_to_pl.txt
│   │   └── pl_to_hy.txt
│   ├── hysc
│   │   ├── hy_to_sc.txt
│   │   └── sc_to_hy.txt
│   └── plsc
│       ├── pl_to_sc.txt
│       └── sc_to_pl.txt
├── planarian.h5ad
└── samap_output.pkl

Samap will run on this data directory and output samap_output.pkl, which is a pickled samap object.