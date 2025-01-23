# AWIGen2Microbiome

Repository for the AWI-Gen 2 microbiome project workflows, analysis, and 
results. 

Please note that this repository aims to document our analysis approach and is
not intended for a full reproduction of the results, since the metadata tables 
are needed to fully run many of the scripts. Given the nature of the project, 
the metadata tables are available from EGA only after approval from 
the relevant data access committee. 

### Reference

Please refer to this publication:

> [Maghini, Ovoduran _et al._ Expanding the Gut Microbiome Atlas of 
Africa __Nature__ 2025](https://www.nature.com/articles/s41586-024-08485-8)

### Repository organization

- `./data/` Folder to hold classification and metadata tables.
- `./files/` Folder to hold derived data tables after running the scripts.
- `./figures/` Folder to hold figures after running the scripts.
- `./src/` R code for the analysis
- `./cluster_scripts/` Other analysis and data-wrangling scripts that were run 
on our cluster
- `./workflows/` Legacy folder for the nextflow processing workflows

## Contact and questions

If you have any questions or comments, please feel free to 
[open an issue](https://github.com/bhattlab/AWIGen2Microbiome/issues/new) or
contact the corresponding author on the above-mentioned publication.
