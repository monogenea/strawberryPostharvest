## Strawberry Post-harvest MS

This repository contains the datasets and R scripts necessary to reproduce the core analyses reported in XXXX et al. 2020.

# Instructions

Run the three R scripts in the following order, after installing the packages `readxl`, `stringr`, `pcaMethods` (Bioconductor), `pheatmap`, `limma` (Bioconductor), `xlsx`, `ellipse` and `RColorBrewer`:

1) `0-preproc.R`, to pre-process the metabolomics data from Tables S1-3
2) `1-heatmapPCA.R`, to produce the heatmap and PCA figures from the pre-processed data
3) `2-limma.R`, to model the pre-processed data using the `limma` R package

The directory `data/` contains the pre-processed data as well as a set of results from the modeling approach. Read the file `sessionInfo` to access exact R package versions and hardware configuration.
