# RecycNecromass

Analysis scripts for the aerobic recycled-necromass study.

## Script organization

- `scripts/physiology/`: physiology, ammonium, and inorganic nitrogen analyses.
- `scripts/16S/`: QIIME 2/phyloseq import, diversity, dbRDA, differential abundance, and aerobic ordinations.
- `scripts/metabolites/`: targeted and untargeted metabolite ordinations, heatmaps, and dbRDA.
- `scripts/community_and_networks/`: genus–environment correlation, community assembly, cross-layer 16S/KO integration, and functional redundancy.
- `scripts/workflow/`: local QIIME 2, PICRUSt2, and metadata-manifest helpers.

The scripts are provided without raw data, generated figures, or local software environments. Before running them, set the project directory and update any input-file paths for the local data layout. Some original analysis scripts still contain legacy absolute paths and should be parameterized before use on another machine.

## Reproducibility notes

Please record the R/Python versions and package environments, document the input files used for each figure, and keep raw sequencing and metabolomics data outside this repository unless redistribution is permitted by the data agreements.
