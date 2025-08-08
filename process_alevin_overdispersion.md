---
title: "Overdispersion Correction for Alevin"
bibliography: Alevin_OD.bib
csl: nature.csl
output:
  html_document:
    toc: true
    toc_depth: 2
    toc_float: true
---

# Processing Alevin Bootstraps for Seurat Integration

Accurate quantification of gene expression from single-cell RNA-seq depends not only on the raw counts, 
but also on accounting for **technical uncertainty** in those counts. Salmon Alevin[@srivastava_bayesian_2020; @patro_salmon_2017; @srivastava_alevin_2019] can output bootstrap 
replicates that capture inferential uncertainty due to multi-mapping reads and fragment assignment, 
analogous to transcript-level bootstraps in bulk RNA-seq. This uncertainty manifests as **overdispersion**, 
which, if uncorrected, can inflate variability estimates and distort downstream analyses such as differential 
expression or trajectory inference.[@baldoni_dividing_2024]

While overdispersion modeling is well established in bulk RNA-seq (e.g., in `edgeR::catchSalmon()`), 
these same principles apply to single-cell data, Here, each **cell** is the unit over which bootstrap 
variance is pooled, rather than the sample. In both cases, the goal is to obtain a moderated per-feature 
overdispersion estimate and use it to **adjust counts** prior to downstream modeling.[@baldoni_dividing_2024]

To perform this procedure, Alevin must be run with parameters that produce per-cell bootstrap replicates 
at the **gene** level. For example:

```bash
salmon alevin \
  -l ISR \
  -1 example_1.fastq.gz \
  -2 example_2.fastq.gz \
  --chromiumV3 \
  -i index/transcripts \
  -p 10 \
  --whitelist index/3M-february-2018.txt \
  --numCellBootstraps 20 \
  --dumpFeatures \
  -o quants/example \
  --tgMap index/tx2g.tsv
```

The key options are:

- `--numCellBootstraps N` — number of bootstrap replicates to generate per cell.
- `--dumpFeatures` — ensures that bootstraps are output at the feature (gene) level, producing `quants_boot_mat.gz`.
- `--tgMap` — transcript-to-gene mapping; required for gene-level quantification.

The remainder of this document walks through an annotated R workflow that:

- Reads **Salmon Alevin** gene-level outputs (counts and bootstrap replicates).
- Computes per-gene overdispersion following the logic of **`edgeR::catchSalmon()`**.
- Divides out the overdispersion from the counts [@baldoni_dividing_2024] prior to creating and merging **Seurat**[@Butler2018Integrating] objects.


---

## 1) Setup

```r
library(Seurat)
library(tximport)
library(Matrix)
library(eds)
library(jsonlite)

# timestamped logger
.log <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", sprintf(...))
```

**Note:** Alevin’s count and bootstrap matrices are stored in the **EDS** format. The `eds` package provides `readEDS()` to load them efficiently.

---

## 2) Helper: read EDS counts (features × cells)

```r
# Helper: read an EDS as features x cells, with logging
read_eds_gc <- function(path, n_feat, n_cells, label) {
  .log("[read] %s expecting features=%d cells=%d (file-native order)", label, n_feat, n_cells)
  m <- eds::readEDS(path, numOfGenes = n_feat, numOfOriginalCells = n_cells)  # returns features x cells
  .log("[read] %s got %d x %d (features x cells)", label, nrow(m), ncol(m))
  if (nrow(m) != n_feat || ncol(m) != n_cells) {
    stop(sprintf("[read] %s dimension mismatch: expected %dx%d, got %dx%d",
                 label, n_feat, n_cells, nrow(m), ncol(m)))
  }
  m
}
```

Alevin’s **counts** (`quants_mat.gz`) are stored as **features × cells** (genes × barcodes). We assert the dimensions match the accompanying `quants_mat_cols.txt` (genes) and `quants_mat_rows.txt` (cells).

---

## 3) Overdispersion from Alevin bootstraps (gene-level)

```r
# Compute overdispersion from raw Alevin bootstrap replicates (gene-level)
compute_overdisp_from_boot_gene <- function(alevin_dir, n_boot, verbose=TRUE) {
  f_boot  <- file.path(alevin_dir, "quants_boot_mat.gz")
  f_cols  <- file.path(alevin_dir, "quants_mat_cols.txt")  # genes (G)
  f_rows  <- file.path(alevin_dir, "quants_mat_rows.txt")  # cells (C)
  stopifnot(file.exists(f_boot), file.exists(f_cols), file.exists(f_rows))
  
  genes <- scan(f_cols, what="character", quiet=TRUE)
  cells <- scan(f_rows, what="character", quiet=TRUE)
  G <- length(genes); C <- length(cells)
  .log("[OD] Expecting boot with C=%d, G=%d, N=%d", C, G, n_boot)
  
  # readEDS returns G x (C*N)
  boot <- eds::readEDS(f_boot, numOfGenes = G, numOfOriginalCells = C * n_boot)
  .log("[OD] boot dims = %d x %d (genes x C*N)", nrow(boot), ncol(boot))
  stopifnot(nrow(boot) == G, ncol(boot) == C * n_boot)
  
  OverDisp_sum <- numeric(G)
  DF <- integer(G)
  eps <- 1e-8
  
  for (ci in seq_len(C)) {
    cols <- ((ci - 1L) * n_boot + 1L):(ci * n_boot)  # N columns for this cell
    Bc <- boot[, cols, drop = FALSE]                  # G x N
    mu <- Matrix::rowMeans(Bc)                        # length G
    s2 <- Matrix::rowSums((Bc - mu)^2) / (n_boot - 1L)
    
    keep <- (mu > 0)
    if (any(keep)) {
      OverDisp_sum[keep] <- OverDisp_sum[keep] + (n_boot - 1L) * (s2[keep] / pmax(mu[keep], eps))
      DF[keep] <- DF[keep] + (n_boot - 1L)
    }
    if (verbose && (ci %% 1000 == 0 || ci == C)) .log("[OD] pooled %d/%d cells", ci, C)
  }
  
  i <- (DF > 0L)
  OverDisp <- rep_len(NA_real_, G)
  if (any(i)) {
    OverDisp[i] <- OverDisp_sum[i] / DF[i]
    DFMedian <- stats::median(DF[i]); DFPrior <- 3
    OverDispPrior <- stats::median(OverDisp[i]) / stats::qf(0.5, df1=DFMedian, df2=DFPrior)
    if (!is.finite(OverDispPrior) || OverDispPrior < 1) OverDispPrior <- 1
    OverDisp[i] <- (DFPrior * OverDispPrior + DF[i] * OverDisp[i]) / (DFPrior + DF[i])
    OverDisp <- pmax(OverDisp, 1)
    OverDisp[!i] <- OverDispPrior
  }
  .log("[OD] OD range: min=%.3f, med=%.3f, max=%.3f",
       suppressWarnings(min(OverDisp, na.rm=TRUE)),
       suppressWarnings(median(OverDisp, na.rm=TRUE)),
       suppressWarnings(max(OverDisp, na.rm=TRUE)))
  
  list(OverDisp=OverDisp, DF=DF, features=genes)
}
```

---

## 4) Comparison to `edgeR::catchSalmon()`

Alevin output is **gene-level**, `quants_boot_mat.gz` has shape **\(G \times (C \cdot N)\)**; each cell contributes **N contiguous columns** (one per bootstrap).

For a gene \( g \) and cell \( c \), let \( b_{g,c,1},\dots,b_{g,c,B} \) be the \( B \) bootstrap replicate counts. Define:

- **Mean** \( \mu_{g,c} = \frac{1}{B} \sum_{k=1}^{B} b_{g,c,k} \).

- **Unbiased variance** \( s^{2}_{g,c} = \frac{1}{B-1} \sum_{k=1}^{B} (b_{g,c,k} - \mu_{g,c})^{2} \).

The **per-cell contribution** to the overdispersion numerator is
\[
\mathrm{OD}^{\mathrm{sum}}_{g,c} \;=\; (B-1)\,\frac{s^{2}_{g,c}}{\mu_{g,c}} \quad\text{for}\;\mu_{g,c} > 0.
\]

Pooling across all cells \( c=1,\dots,C \), we accumulate
\[
\mathrm{OD}^{\mathrm{sum}}_{g} \;=\; \sum_{c=1}^{C} \mathrm{OD}^{\mathrm{sum}}_{g,c}, 
\qquad
\mathrm{DF}_{g} \;=\; \sum_{c=1}^{C} (B-1).
\]

The **raw** overdispersion estimate is the per-gene average:
\[
\widehat{\mathrm{OD}}_{g} \;=\; \frac{\mathrm{OD}^{\mathrm{sum}}_{g}}{\mathrm{DF}_{g}}.
\]

### Moderation (same as `edgeR::catchSalmon()`)

Let \( \mathrm{DF}_{\text{median}} \) be the median of \( \mathrm{DF}_{g} \) over genes with \( \mathrm{DF}_{g} > 0 \). Let the prior degrees of freedom be \( \mathrm{DF}_{\text{prior}} = 3 \). Define the prior overdispersion as:
\[
\mathrm{OD}_{\text{prior}} \;=\; \frac{\operatorname{median}\!\left(\widehat{\mathrm{OD}}_{g}\right)}
                                      {F^{-1}_{\,\text{df1}=\mathrm{DF}_{\text{median}},\;\text{df2}=\mathrm{DF}_{\text{prior}}}(0.5)},
\]
where \(F^{-1}\) is the inverse CDF of the F distribution.

Then **moderate** each gene:
\[
\mathrm{OD}^{\text{mod}}_{g} \;=\;
\frac{\mathrm{DF}_{\text{prior}}\cdot \mathrm{OD}_{\text{prior}} + \mathrm{DF}_{g}\cdot \widehat{\mathrm{OD}}_{g}}
     {\mathrm{DF}_{\text{prior}} + \mathrm{DF}_{g}},
\qquad
\mathrm{OD}^{\text{mod}}_{g} \leftarrow \max\!\big(\mathrm{OD}^{\text{mod}}_{g}, 1\big).
\]

For genes with \( \mathrm{DF}_{g} = 0 \), set \( \mathrm{OD}^{\text{mod}}_{g} = \mathrm{OD}_{\text{prior}} \).  
These are **exactly the same** moderation and clamping rules used by `edgeR::catchSalmon()`.

### Correction applied to counts

Let \( \mathrm{counts}_{g,c} \) be the original count for gene \( g \) in cell \( c \). The **corrected** count is:
\[
\mathrm{counts}^{\text{corr}}_{g,c} \;=\; \frac{\mathrm{counts}_{g,c}}{\mathrm{OD}^{\text{mod}}_{g}}.
\]

This is implemented by left-multiplying the counts matrix by a diagonal matrix with entries \( 1 / \mathrm{OD}^{\text{mod}}_{g} \).

---

## 5) Sample processing function

```r
# Function to build a Seurat object with counts / overdispersion
process_sample <- function(sample_id, base_dir, n_boot) {
  t0 <- Sys.time()
  alevin_dir <- file.path(base_dir, sample_id, "alevin")
  .log("[%s] Starting in %s", sample_id, alevin_dir)
  
  f_counts <- file.path(alevin_dir, "quants_mat.gz")
  f_cols   <- file.path(alevin_dir, "quants_mat_cols.txt") # genes
  f_rows   <- file.path(alevin_dir, "quants_mat_rows.txt") # cells
  
  stopifnot(file.exists(f_counts), file.exists(f_cols), file.exists(f_rows))
  
  feats <- scan(f_cols, what="character", quiet=TRUE)
  cells <- scan(f_rows, what="character", quiet=TRUE)
  .log("[%s] Features=%d (genes) Cells=%d", sample_id, length(feats), length(cells))
  
  counts <- read_eds_gc(f_counts,
                        n_feat = length(feats),
                        n_cells = length(cells),
                        label = sprintf("%s:counts", sample_id))
  
  stopifnot(nrow(counts) == length(feats), ncol(counts) == length(cells))
  rownames(counts) <- feats
  colnames(counts) <- paste0(cells, "_", sample_id)
  
  # gene-level OD from boots, then counts/OD
  od <- compute_overdisp_from_boot_gene(alevin_dir, n_boot = n_boot, verbose = TRUE)
  
  stopifnot(identical(rownames(counts), od$features))
  inv_od <- 1 / od$OverDisp
  inv_od[!is.finite(inv_od)] <- 1
  .log("[%s] Scaling counts by 1/OD", sample_id)
  counts_corr <- Matrix::Diagonal(x = inv_od) %*% counts
  
  # reattach gene + cell names in case the sparse op dropped them
  dimnames(counts_corr) <- dimnames(counts)
  
  seu <- Seurat::CreateSeuratObject(counts = counts_corr, project = sample_id)
  DefaultAssay(seu) <- "RNA"
  seu$sample <- sample_id
  .log("[%s] Done in %.1fs | %d genes x %d cells", sample_id,
       as.numeric(difftime(Sys.time(), t0, "secs")), nrow(seu), ncol(seu))
  seu
}
```

---

## 6) Example for processing and merging multiple samples

```r
dir_out <- "/rsrch5/home/epi/bhattacharya_lab/data/Placenta_scRNAseq/analysis"
dir_data <- "/rsrch5/home/epi/bhattacharya_lab/data/Placenta_scRNAseq/salmon/quants"
SRA <- c("SRR14134833","SRR14134834","SRR14134836","SRR15193608")

files <- paste(dir_data,SRA,"alevin/quants_mat.gz",sep="/")

objs <- lapply(SRA, function(s) process_sample(s, base_dir = dir_data, n_boot = 20))

.log("Merging %d Seurat objects", length(objs))
combined <- Reduce(function(x, y) merge(x, y), objs)
.log("Merged object: %d features x %d cells", nrow(combined), ncol(combined))
```

---

## References