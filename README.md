# The Circadian Clock in the Ageing Plant

This repository contains the code to run the analysis and create figures for this project.

## Installation

This project uses [`renv`](https://rstudio.github.io/renv/articles/renv.html). `renv` should automatically set up and download required packages when you first load this project.

Follow [this guide](https://rstudio.github.io/renv/articles/collaborating.html) for more instructions.

## Folder structure

-   `data`
    -   `chip_seq` = data of targets of core circadian clock genes (CCA1, ELF3, ELF4, LUX, LHY, PRRs)
    -   `gene_lists` = collections of genes from other sources, such as circadian clock genes, P450 enzymes, and photoperiodic pathway genes
    -   `other_datasets` = gene expression data from other publications (currently just Mockler lab)
    -   `quants` = outputs from Salmon from my RNA-seq samples
-   `R` = the code
-   `outputs` = final figures / supplementary figures for the paper plus any additional data outputs
-   `pre_compute` = (You should create this as an empty directory after initialising the repo!) Outputs from my R code that require quite a lot of compute time or call external databases (so they should not be run every time the code is called)
-   `renv` = the folder controlled by renv to create a reproducible R environment with specific versions of packages
