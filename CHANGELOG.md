# nf-core/eisca: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.7.0 - [2026-07-16]
- add a differential expression analysis of pseudobulk samples using package PyDeseq2 

## v2.6.0 - [2026-05-20]
- add a new module for CellBender to remove empty droplets.

## v2.5.0 - [2025-11-19]
- add the scvi integration method to the clustering module using scvi-tools
- add a new module annotation_scvi for cell-type prediction using scvi-tools
- add a new module dea_scvi for differential expression analysis using scvi-tools

## v2.0.0 - [2025-07-01]
- add cell-cell communication analysis with tool cellchat
- allow genearting PDF figures

## v1.8.0 - [2024-12-11]
- add support for smart-seq2 data analysis

## v1.6.0 - [2024-11-05]
- add tertiary analysis: Differential expression analysis
- fix plot size issue in report


## v1.4.0 - [2024-10-22]
- add tertiary analysis: Cell type annotation
- add 'analysis parameters' for analysis sections in report 

## v1.1.0 - [2024-10-04]
- add module 'Merging/integration of samples' for secondary analysis

## v1.0.0 - [2024-09-19]

Initial release of EarlhamInst/eisca, created with the [nf-core](https://nf-co.re/) template.
The pipeline including following analyses:
1. Primary analysis:
   - Quality control of raw reads
   - Mapping reads and quantification
   - Convert count matrix into Anndata and Securat objects
2. Secondary analysis
   - Single-cell quality control
   - Cell filtering
   - Clustering analysis

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
