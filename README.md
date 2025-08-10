# Processing Alevin Bootstraps for Seurat Integration

This repository contains an R workflow and supporting code for processing Salmon Alevin bootstrap replicates to account for technical uncertainty before integrating single-cell RNA-seq datasets in Seurat. [See the workflow here](http://seantbresnahan.com/Overdispersion-Correction-for-Alevin).

When quantifying single-cell transcript expression, fragment assignment ambiguity and multi-mapping reads can introduce **inferential variance** in estimated counts.  
Salmon Alevin can produce per-cell bootstrap replicates, analogous to transcript-level bootstraps in bulk RNA-seq, which capture this uncertainty.  

If left uncorrected, this overdispersion can inflate variance estimates, biasing downstream analyses such as differential expression, clustering, and trajectory inference.

Overdispersion modeling is widely used in bulk RNA-seq (e.g., `edgeR::catchSalmon()`), but here the approach is adapted for single-cell data by treating cells, rather than bulk samples, as the unit over which bootstrap variance is estimated.  

The goal is to compute a moderated per-transcript overdispersion estimate and use it to adjust counts prior to downstream modeling, improving biological signal detection.

---

## Contents

- **Example Salmon Alevin command** for transcriptrating per-cell transcript-level bootstraps, highlighting required parameters:
  - `--numCellBootstraps` — number of replicates per cell.
  - `--dumpFeatures` — output transcript-level bootstraps (`quants_boot_mat.gz`).
  - `--tgMap` — transcript-to-self (rather than transcript-to-gene) mapping table; required for transcript-level quantification.

- **R functions to:**
  - Read Alevin transcript-level counts and bootstrap replicates.
  - Compute bootstrap-derived per-transcript overdispersion estimates following the `edgeR::catchSalmon()` methodology.
  - Adjust counts by dividing out overdispersion, as in *Baldoni et al.* (2024), before Seurat integration.

- **Annotated example workflow showing:**
  - Import of Alevin outputs.
  - Overdispersion correction.
  - Creation and merging of corrected Seurat objects.

---
