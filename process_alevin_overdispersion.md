---
title: "Correcting for overdispersion from read-to-transcript mapping uncertainty in scRNA-seq"
bibliography: Alevin_OD.bib
csl: nature.csl
output:
  html_document:
    toc: true
    toc_depth: 2
    toc_float: true
---

# Processing Alevin bootstraps for Seurat integration

Accurate quantification of transcript expression from single-cell RNA-seq depends not only on the raw counts, but also on accounting for **technical uncertainty** in those counts. Salmon Alevin can output bootstrap replicates that capture inferential uncertainty due to fragments mapping to multiple transcripts of the same gene.[@srivastava_bayesian_2020; @patro_salmon_2017; @srivastava_alevin_2019] This uncertainty manifests as **overdispersion**, which, if uncorrected, can inflate variability estimates and distort downstream analyses (differential expression, clustering, trajectory inference).[@baldoni_dividing_2024]

While overdispersion modeling is well established in bulk RNA-seq (e.g., `edgeR::catchSalmon()`), here we adapt the approach for single-cell data by treating **cells** as the unit over which bootstrap variance is estimated. The goal is to compute a moderated per-transcript overdispersion estimate and use it to **adjust counts** prior to downstream modeling with Seurat[@Butler2018Integrating], improving biological signal detection.[@baldoni_dividing_2024]

## Dataset for example analysis

For this tutorial, we use scRNA-seq from **four healthy term placentas**,[@lu-culligan2021maternal] specifically: SRR14134833, SRR14134834, SRR14134836, and SRR15193608. Libraries consisting of approximately 4.8M cells across these four samples were prepared using the 10x Genomics Single Cell 3′ Reagent Kit v3 and sequenced on an Illumina NovaSeq platform.

Quantification in this tutorial uses our **long-read defined, placenta-specific transcriptome reference**.[@Bresnahan_placenta] Villous placenta RNA was used to prepare ONT cDNA libraries (PCS111), with raw signals basecalled by Guppy. Full-length cDNA reads were oriented and trimmed with Pychopper2, then corrected for microindels using FMLRC2 with matched short-read placenta RNA-seq. Corrected reads were aligned to hg38 using minimap2; primary alignments across all samples were concatenated and assembled de novo with ESPRESSO. The assembly was characterized with SQANTI3 against GENCODE v45 to assign structural classes, and high-confidence isoforms were identified using multiple orthogonal support signals including long-read coverage, short-read splice-junction support, and proximity of TSS/TTS to placenta CAGE-seq, DNase-seq, and SAPAS sites.

## Alevin parameters

Run Alevin to quantify transcript expression and produce per-cell transcript-level bootstrap replicates:

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
  --tgMap index/tx2g.tsv # tx-to-tx identity table
```

Key options:

- `--numCellBootstraps N` — number of bootstrap replicates per cell.
- `--dumpFeatures` — outputs bootstraps at the feature (transcript) level, producing `quants_boot_mat.gz`.
- `--tgMap` — transcript-to-self (rather than transcript-to-gene) mapping table; required for transcript-level quantification.

[You can find the official Alevin documentation here.](https://salmon.readthedocs.io/en/latest/alevin.html)

---

# 1) Setup

```r
suppressPackageStartupMessages({
  library(Seurat)
  library(edgeR)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(eds)
  library(jsonlite)
  library(tximport)
})

# timestamped logger
.log <- function(...) message(format(Sys.time(), "%H:%M:%S"), " | ", sprintf(...))
```

**Notes**: Seurat v5 prefers `layer = "counts"` instead of the deprecated `slot` argument in `GetAssayData()`. Alevin’s EDS matrices are efficiently read with `eds::readEDS()`.

---

# 2) Helper: read EDS counts (features × cells)

```r
# Helper: read an EDS as features x cells, with logging
read_eds_gc <- function(path, n_feat, n_cells, label) {
  .log("[read] %s expecting features=%d cells=%d (file-native order)", label, n_feat, n_cells)
  m <- eds::readEDS(path, numOftranscripts = n_feat, numOfOriginalCells = n_cells)  # returns features x cells
  .log("[read] %s got %d x %d (features x cells)", label, nrow(m), ncol(m))
  if (nrow(m) != n_feat || ncol(m) != n_cells) {
    stop(sprintf("[read] %s dimension mismatch: expected %dx%d, got %dx%d",
                 label, n_feat, n_cells, nrow(m), ncol(m)))
  }
  m
}
```

---

# 3) Fast parallel overdispersion from Alevin bootstraps

This implementation streams transcript-by-cell bootstrap counts, accumulates per-transcript overdispersion across cells in blocks, and then applies the same moderation/clamping used in catchSalmon-style approaches.

```r
compute_overdisp_from_boot_transcript_fast_par <- function(
    alevin_dir,
    n_boot,
    block_cells = 128L,
    n_cores = max(1L, parallel::detectCores() - 1L),
    omp_threads = 1L,
    verbose = TRUE
){
  suppressPackageStartupMessages({ library(Matrix) })
  
  f_boot <- file.path(alevin_dir, "quants_boot_mat.gz")
  f_cols <- file.path(alevin_dir, "quants_mat_cols.txt")
  f_rows <- file.path(alevin_dir, "quants_mat_rows.txt")
  stopifnot(file.exists(f_boot), file.exists(f_cols), file.exists(f_rows))
  
  transcripts <- scan(f_cols, what = "character", quiet = TRUE)
  cells <- scan(f_rows, what = "character", quiet = TRUE)
  G <- length(transcripts); C <- length(cells)
  .log("[OD] Expecting boot with C=%d, G=%d, N=%d", C, G, n_boot)
  
  boot <- eds::readEDS(f_boot, numOftranscripts = G, numOfOriginalCells = C * n_boot)
  if (!inherits(boot, "dgCMatrix")) boot <- as(boot, "dgCMatrix")
  .log("[OD] boot dims = %d x %d (transcripts x C*N)", nrow(boot), ncol(boot))
  
  # Keep per-worker threading low to avoid oversubscription
  Sys.setenv(OMP_NUM_THREADS = as.integer(omp_threads))
  
  nbm1 <- as.integer(n_boot - 1L)
  eps  <- 1e-8
  
  # (bs*N) x bs grouping matrix: sums N bootstrap columns per cell
  make_Gmat <- function(bs) {
    Matrix::sparseMatrix(
      i = as.integer(seq_len(bs * n_boot)),
      j = as.integer(rep(seq_len(bs), each = n_boot)),
      x = 1,
      dims = c(bs * n_boot, bs)
    )
  }
  Gmat_full <- make_Gmat(block_cells)
  
  square_sparse <- function(M) { M2 <- M; M2@x <- M2@x^2; M2 }
  
  starts <- seq.int(1L, C, by = block_cells)
  blocks <- lapply(starts, function(s) {
    e <- min(C, s + block_cells - 1L)
    list(start = as.integer(s), end = as.integer(e))
  })
  
  # Worker over a block of cells
  block_fn <- function(bl) {
    start <- bl$start; end <- bl$end
    bs <- as.integer(end - start + 1L)
    
    j1 <- (start - 1L) * n_boot + 1L
    j2 <- end * n_boot
    if (j1 > ncol(boot) || j2 > ncol(boot)) stop("Block column indices out of bounds")
    
    Bblk <- boot[, j1:j2, drop = FALSE]                 # G x (bs*N), sparse
    Gmat <- if (bs == block_cells) Gmat_full else make_Gmat(bs)
    
    # Sum and sumsq across N boot reps per cell (fast sparse %*%)
    S  <- Bblk %*% Gmat                                  # G x bs
    SS <- square_sparse(Bblk) %*% Gmat                   # G x bs
    
    # Coerce to base matrices for stable elementwise ops / rowSums
    S  <- as.matrix(S)
    SS <- as.matrix(SS)
    
    mu <- S / n_boot
    s2 <- (SS - n_boot * (mu * mu)) / nbm1
    
    keep <- (mu > 0)                                    # logical G x bs
    
    # Contribution: nbm1 * s2 / pmax(mu, eps) where keep; 0 elsewhere
    R <- nbm1 * (s2 / pmax(mu, eps))
    R[!keep] <- 0
    
    list(
      OD_sum = rowSums(R),                  # length G (numeric)
      DF_sum = nbm1 * rowSums(keep),        # length G (integer-ish)
      upto   = as.integer(end)
    )
  }
  
  # Parallel over blocks (forked workers share 'boot' on Unix)
  out_list <- parallel::mclapply(blocks, block_fn, mc.cores = n_cores)
  
  OverDisp_sum <- Reduce(`+`, lapply(out_list, `[[`, "OD_sum"))
  DF           <- Reduce(`+`, lapply(out_list, `[[`, "DF_sum"))
  
  if (verbose) .log("[OD] pooled %d/%d cells", max(vapply(out_list, `[[`, 1L, "upto")), C)
  
  # Shrinkage + floor (same as your original logic)
  i <- (DF > 0L)
  OverDisp <- rep_len(NA_real_, G)
  if (any(i)) {
    OverDisp[i] <- OverDisp_sum[i] / DF[i]
    DFMedian <- stats::median(DF[i]); DFPrior <- 3L
    OverDispPrior <- stats::median(OverDisp[i]) / stats::qf(0.5, df1 = DFMedian, df2 = DFPrior)
    if (!is.finite(OverDispPrior) || OverDispPrior < 1) OverDispPrior <- 1
    OverDisp[i] <- (DFPrior * OverDispPrior + DF[i] * OverDisp[i]) / (DFPrior + DF[i])
    OverDisp <- pmax(OverDisp, 1)
    OverDisp[!i] <- OverDispPrior
  }
  
  .log("[OD] OD range: min=%.3f, med=%.3f, max=%.3f",
       suppressWarnings(min(OverDisp, na.rm = TRUE)),
       suppressWarnings(median(OverDisp, na.rm = TRUE)),
       suppressWarnings(max(OverDisp, na.rm = TRUE)))
  
  list(OverDisp = OverDisp, DF = DF, features = transcripts)
}
```

---

# 4) Build Seurat objects with raw and corrected assays

Creates two assays per sample: the raw counts in `RNA` and the overdispersion-corrected counts in `RNA_corr`. Feature-level `OverDisp` and `DF` are stored in each assay’s `Misc("feature_meta")`.

```r
process_sample <- function(sample_id, base_dir, n_boot,
                           block_cells = 128L,
                           n_cores = max(1L, parallel::detectCores() - 1L),
                           omp_threads = 1L) {
  t0 <- Sys.time()
  alevin_dir <- file.path(base_dir, sample_id, "alevin")
  .log("[%s] Starting in %s", sample_id, alevin_dir)
  
  f_counts <- file.path(alevin_dir, "quants_mat.gz")
  f_cols   <- file.path(alevin_dir, "quants_mat_cols.txt")
  f_rows   <- file.path(alevin_dir, "quants_mat_rows.txt")
  stopifnot(file.exists(f_counts), file.exists(f_cols), file.exists(f_rows))
  
  feats <- scan(f_cols, what="character", quiet=TRUE)
  cells <- scan(f_rows, what="character", quiet=TRUE)
  .log("[%s] Features=%d (transcripts) Cells=%d", sample_id, length(feats), length(cells))
  
  counts <- read_eds_gc(
    f_counts,
    n_feat  = length(feats),
    n_cells = length(cells),
    label   = sprintf("%s:counts", sample_id)
  )
  stopifnot(nrow(counts) == length(feats), ncol(counts) == length(cells))
  rownames(counts) <- feats
  colnames(counts) <- paste0(cells, "_", sample_id)
  
  # Parallel transcript-level OD from boots
  od <- compute_overdisp_from_boot_transcript_fast_par(
    alevin_dir  = alevin_dir,
    n_boot      = n_boot,
    block_cells = block_cells,
    n_cores     = n_cores,
    omp_threads = omp_threads,
    verbose     = TRUE
  )
  
  stopifnot(identical(rownames(counts), od$features))
  inv_od <- 1 / od$OverDisp
  inv_od[!is.finite(inv_od)] <- 1
  .log("[%s] Scaling counts by 1/OD", sample_id)
  counts_corr <- Matrix::Diagonal(x = inv_od) %*% counts
  dimnames(counts_corr) <- dimnames(counts)
  
  # Build Seurat with two assays (raw + corrected)
  seu <- Seurat::CreateSeuratObject(counts = counts, project = sample_id, assay = "RNA")
  seu[["RNA_corr"]] <- Seurat::CreateAssayObject(counts = counts_corr)
  
  # Store feature-level OD/DF in Misc (Seurat v5-safe)
  fm <- data.frame(OverDisp = as.numeric(od$OverDisp),
                   DF       = as.integer(od$DF),
                   row.names = feats)
  SeuratObject::Misc(seu[["RNA"]],      "feature_meta") <- fm
  SeuratObject::Misc(seu[["RNA_corr"]], "feature_meta") <- fm
  
  # Sample label; choose default assay
  DefaultAssay(seu) <- "RNA_corr"
  seu$sample <- sample_id
  
  .log("[%s] Done in %.1fs | %d transcripts x %d cells", sample_id,
       as.numeric(difftime(Sys.time(), t0, "secs")), nrow(seu), ncol(seu))
  seu
}
```

---

# 5) Run across samples and merge

```r
dir_out <- "/rsrch5/home/epi/bhattacharya_lab/data/Placenta_scRNAseq/analysis"
dir_data <- "/rsrch5/home/epi/bhattacharya_lab/data/Placenta_scRNAseq/salmon/quants"
SRA <- c("SRR14134833","SRR14134834","SRR14134836","SRR15193608")
files <- paste(dir_data,SRA,"alevin/quants_mat.gz",sep="/")
file.exists(files)

objs <- lapply(
  SRA,
  function(s) process_sample(
    sample_id   = s,
    base_dir    = dir_data,
    n_boot      = 20,
    block_cells = 128,  # safe block size for big G, tune upward if RAM allows
    n_cores     = 8,    # leave ~2 cores headroom for I/O, GC, OS
    omp_threads = 1     # one BLAS/OpenMP thread per worker to avoid oversubscription
  )
)

.log("Merging %d Seurat objects", length(objs))
combined <- Reduce(function(x, y) merge(x, y), objs)
.log("Merged object: %d features x %d cells", nrow(combined), ncol(combined))
```

---

# 6) Between-sample BCV (pseudo-bulk per sample) and plots

We compute between-sample BCV by first pseudo-bulking counts by a `sample_id` (or other grouping), then running edgeR’s dispersion estimation. Plots show tagwise BCV vs average log-CPM with common dispersion (red) and trend (blue).

```r
# Sum cell counts to pseudo-bulk by a metadata column (e.g., sample_id)
pseudobulk_by <- function(
    seurat_obj,
    sample_col = "sample_id",
    assay = "RNA",
    layer = "counts",
    min_cells_per_sample = 1
){
  .log("[pseudobulk_by] assay=%s layer=%s sample_col=%s", assay, layer, sample_col)
  
  # --- metadata checks ---
  md <- seurat_obj[[]]
  if (!sample_col %in% names(md)) {
    stop(sprintf("Metadata column '%s' not found. Available: %s",
                 sample_col, paste(head(names(md), 25), collapse = ", ")))
  }
  aobj <- seurat_obj[[assay]]
  if (is.null(aobj)) stop(sprintf("Assay '%s' not found in object.", assay))
  
  # aggregate one matrix X (genes x cells) by group ids
  agg_one <- function(X, ids, min_cells) {
    if (!inherits(X, "dgCMatrix")) X <- as(Matrix::Matrix(X, sparse = TRUE), "dgCMatrix")
    # build group sizes (ignore NAs)
    tab <- table(ids, useNA = "no")
    keep_groups <- names(tab)[tab >= min_cells]
    if (length(keep_groups) == 0L) {
      return(Matrix::Matrix(0, nrow = nrow(X), ncol = 0, sparse = TRUE,
                            dimnames = list(rownames(X), NULL)))
    }
    keep_cells <- ids %in% keep_groups
    if (!any(keep_cells)) {
      return(Matrix::Matrix(0, nrow = nrow(X), ncol = 0, sparse = TRUE,
                            dimnames = list(rownames(X), NULL)))
    }
    X <- X[, keep_cells, drop = FALSE]
    grp <- factor(ids[keep_cells], levels = keep_groups)
    
    # cells x groups indicator (no contrasts)
    jj <- as.integer(grp)
    M <- Matrix::sparseMatrix(
      i = seq_along(jj),
      j = jj,
      x = 1,
      dims = c(length(jj), length(levels(grp))),
      dimnames = list(NULL, levels(grp))
    )
    
    PB <- X %*% M  # genes x groups (dgCMatrix)
    if (inherits(PB, "dgCMatrix")) PB@x <- as.numeric(PB@x)  # ensure numeric
    PB
  }
  
  # Helper: robustly map ids for a given matrix
  map_ids <- function(X, sample_col) {
    # 1) Best: FetchData on exactly these cells
    ids <- tryCatch(
      SeuratObject::FetchData(seurat_obj, vars = sample_col, cells = colnames(X))[[1]],
      error = function(e) rep(NA_character_, ncol(X))
    )
    # 2) Fallback: rowname match into metadata
    if (all(is.na(ids))) {
      rn <- rownames(md)
      if (!is.null(rn) && length(rn)) {
        idx <- match(colnames(X), rn)
        ids <- md[[sample_col]][idx]
      }
    }
    ids
  }
  
  # --- Seurat v5 (Assay5) with layers ---
  if (inherits(aobj, "Assay5")) {
    lys <- SeuratObject::Layers(aobj)
    if (length(lys) == 0L) stop(sprintf("Assay '%s' has no layers.", assay))
    
    # If user requested a specific layer (or there is only one), use it
    if (length(lys) == 1L || (layer %in% lys && length(lys) > 1L)) {
      use_layer <- if (layer %in% lys) layer else lys[[1]]
      X <- SeuratObject::GetAssayData(aobj, layer = use_layer)
      ids <- map_ids(X, sample_col)
      
      # If mapping still failed (all NA), as a last resort use a constant from layer name
      if (all(is.na(ids))) {
        # Derive something from the layer string (e.g., "counts" or "counts.SRR...")
        fallback <- sub("^counts\\.?","", use_layer)
        if (fallback == "") fallback <- use_layer
        ids <- rep(fallback, ncol(X))
        .log("[pseudobulk_by] Using layer-derived sample id '%s' for single-layer mapping fallback.", fallback)
      }
      
      PB <- agg_one(X, ids, min_cells_per_sample)
      .log("[pseudobulk_by] Output pseudo-bulk: %d genes x %d samples", nrow(PB), ncol(PB))
      return(PB)
    }
    
    # Multiple layers: aggregate each, then sum by sample name
    PB_list <- vector("list", length(lys))
    names(PB_list) <- lys
    for (i in seq_along(lys)) {
      ly <- lys[i]
      X  <- SeuratObject::GetAssayData(aobj, layer = ly)
      
      ids_try <- map_ids(X, sample_col)
      if (all(is.na(ids_try))) {
        # Fallback to layer-derived id (works with counts.SRR... style)
        sample_from_layer <- sub("^counts\\.?","", ly)
        if (sample_from_layer == "") sample_from_layer <- ly
        ids_try <- rep(sample_from_layer, ncol(X))
        .log("[pseudobulk_by] Layer '%s': using layer-derived sample id '%s' (metadata mapping failed).",
             ly, sample_from_layer)
      }
      PB_list[[i]] <- agg_one(X, ids_try, min_cells_per_sample)
    }
    
    PB_list <- Filter(function(M) ncol(M) > 0, PB_list)
    if (length(PB_list) == 0L) {
      assay_genes <- rownames(aobj)
      .log("[pseudobulk_by] No pseudo-bulk columns produced; try lowering min_cells_per_sample.")
      return(Matrix::Matrix(0, nrow = length(assay_genes), ncol = 0, sparse = TRUE,
                            dimnames = list(assay_genes, NULL)))
    }
    
    common_genes <- rownames(aobj)  # shared across layers in Assay5
    all_samples  <- sort(unique(unlist(lapply(PB_list, colnames))))
    PB_sum <- Matrix::Matrix(0, nrow = length(common_genes), ncol = length(all_samples), sparse = TRUE,
                             dimnames = list(common_genes, all_samples))
    for (PB in PB_list) {
      PB <- PB[common_genes, , drop = FALSE]
      if (inherits(PB, "dgCMatrix")) PB@x <- as.numeric(PB@x)
      PB_sum[, colnames(PB)] <- PB_sum[, colnames(PB), drop = FALSE] + PB
    }
    .log("[pseudobulk_by] Output pseudo-bulk: %d genes x %d samples", nrow(PB_sum), ncol(PB_sum))
    return(PB_sum)
  }
  
  # --- Seurat v4 path (slot-based) ---
  X <- SeuratObject::GetAssayData(seurat_obj, assay = assay, slot = layer)
  # Prefer rowname match into metadata (robust to order)
  ids <- {
    rn <- rownames(md)
    if (!is.null(rn) && length(rn)) md[colnames(X), sample_col, drop = TRUE] else md[[sample_col]][colnames(X)]
  }
  PB <- agg_one(X, ids, min_cells_per_sample)
  .log("[pseudobulk_by] Output pseudo-bulk: %d genes x %d samples", nrow(PB), ncol(PB))
  PB
}

# edgeR BCV on a genes x samples count matrix
bcv_from_matrix <- function(
    PB, min_samples = 2, min_counts = 10, min_samples_expr = 1, robust = TRUE
){
  .log("[bcv_from_matrix] Starting: %d genes x %d samples", nrow(PB), ncol(PB))
  if (!inherits(PB, "dgCMatrix")) PB <- as(Matrix::Matrix(PB, sparse = TRUE), "dgCMatrix")
  
  # gene filter: enough counts overall and expressed in >= min_samples_expr samples
  keep_g <- (Matrix::rowSums(PB) >= min_counts) &
    (Matrix::rowSums(PB > 0) >= min_samples_expr)
  .log("[bcv_from_matrix] Genes passing filter: %d / %d", sum(keep_g), length(keep_g))
  PB <- PB[keep_g, , drop = FALSE]
  
  # sample filter: non-zero libs
  lib <- Matrix::colSums(PB)
  keep_s <- lib > 0
  .log("[bcv_from_matrix] Samples nonzero lib: %d / %d", sum(keep_s), length(keep_s))
  PB <- PB[, keep_s, drop = FALSE]
  if (ncol(PB) < max(2L, min_samples) || nrow(PB) < 2L)
    stop("[bcv_from_matrix] Not enough genes/samples after filtering")
  
  dge <- DGEList(counts = PB, group = factor(rep(1, ncol(PB))))
  dge <- calcNormFactors(dge, method = "TMMwsp")
  design <- matrix(1, nrow = ncol(PB), ncol = 1)
  dge <- estimateDisp(dge, design = design, robust = robust)
  .log("[bcv_from_matrix] Median tagwise BCV=%.4f (n=%d genes, %d samples)",
       median(sqrt(dge$tagwise.dispersion), na.rm = TRUE),
       nrow(dge$counts), ncol(dge$counts))
  
  list(
    ave_log_cpm = dge$AveLogCPM,
    bcv_tagwise = sqrt(dge$tagwise.dispersion),
    bcv_common  = sqrt(dge$common.dispersion),
    bcv_trend   = sqrt(dge$trended.dispersion),
    n_genes     = nrow(dge$counts),
    n_samples   = ncol(dge$counts),
    genes       = rownames(dge$counts),
    samples     = colnames(dge$counts)
  )
}

plot_bcv <- function(bcv, title = "") {
  .log("[plot_bcv] %s: plotting %d genes", title, length(bcv$bcv_tagwise))
  df <- data.frame(ave = bcv$ave_log_cpm, bcv = bcv$bcv_tagwise)
  ggplot(df, aes(x = ave, y = bcv)) +
    geom_point(size = 0.3, alpha = 0.6) +
    geom_hline(yintercept = bcv$bcv_common, color = "red", linewidth = 0.8) +
    geom_line(
      data = data.frame(
        ave = sort(bcv$ave_log_cpm),
        bcv_trend = bcv$bcv_trend[order(bcv$ave_log_cpm)]
      ),
      aes(x = ave, y = bcv_trend),
      color = "blue", linewidth = 0.8
    ) +
    labs(x = "Average log CPM", y = "BCV (sqrt NB dispersion)", title = title) +
    theme_classic()
}

# Compare pre vs post *between-sample* BCV (pseudo-bulk per sample)
compare_bcv_between_samples <- function(
    combined,
    sample_col = "sample_id",
    pre_assay = "RNA",
    post_assay = "RNA_corr",
    layer = "counts",
    min_cells_per_sample = 1,
    min_counts = 10,
    min_samples_expr = 1,
    robust = TRUE
){
  .log("[compare_bcv_between_samples] sample_col=%s | %s vs %s",
       sample_col, pre_assay, post_assay)
  
  PBpre  <- pseudobulk_by(combined, sample_col, assay = pre_assay,  layer = layer,
                          min_cells_per_sample = min_cells_per_sample)
  PBpost <- pseudobulk_by(combined, sample_col, assay = post_assay, layer = layer,
                          min_cells_per_sample = min_cells_per_sample)
  
  if (ncol(PBpre) == 0L)
    stop(sprintf("[compare_bcv_between_samples] No pseudo-bulk samples from '%s'. Check layer names, mapping of '%s', or lower min_cells_per_sample.", pre_assay, sample_col))
  if (ncol(PBpost) == 0L)
    stop(sprintf("[compare_bcv_between_samples] No pseudo-bulk samples from '%s'. Check layer names, mapping of '%s', or lower min_cells_per_sample.", post_assay, sample_col))
  
  # Align samples and genes
  common_samples <- intersect(colnames(PBpre), colnames(PBpost))
  PBpre  <- PBpre[,  common_samples, drop = FALSE]
  PBpost <- PBpost[, common_samples, drop = FALSE]
  common_genes <- intersect(rownames(PBpre), rownames(PBpost))
  PBpre  <- PBpre[ common_genes, , drop = FALSE]
  PBpost <- PBpost[common_genes, , drop = FALSE]
  .log("[compare_bcv_between_samples] Using %d genes x %d samples", nrow(PBpre), ncol(PBpre))
  
  bcv_pre  <- bcv_from_matrix(PBpre,  min_counts = min_counts,
                              min_samples = 2, min_samples_expr = min_samples_expr,
                              robust = robust)
  bcv_post <- bcv_from_matrix(PBpost, min_counts = min_counts,
                              min_samples = 2, min_samples_expr = min_samples_expr,
                              robust = robust)
  
  p1 <- plot_bcv(bcv_pre,  sprintf("Between-sample BCV (pre: %s)", pre_assay))
  p2 <- plot_bcv(bcv_post, sprintf("Between-sample BCV (post: %s)", post_assay))
  combo <- p1 | p2
  
  common_g <- intersect(bcv_pre$genes, bcv_post$genes)
  d <- bcv_pre$bcv_tagwise[match(common_g, bcv_pre$genes)] -
    bcv_post$bcv_tagwise[match(common_g, bcv_post$genes)]
  summary <- list(
    n_samples = ncol(PBpre),
    n_genes   = length(common_g),
    median_BCV_before = median(bcv_pre$bcv_tagwise[match(common_g, bcv_pre$genes)], na.rm = TRUE),
    median_BCV_after  = median(bcv_post$bcv_tagwise[match(common_g, bcv_post$genes)], na.rm = TRUE),
    median_delta      = median(d, na.rm = TRUE)
  )
  .log("[compare_bcv_between_samples] Median BCV pre=%.4f post=%.4f Δ=%.4f",
       summary$median_BCV_before, summary$median_BCV_after, summary$median_delta)
  
  list(plots = combo, pre = bcv_pre, post = bcv_post, summary = summary)
}
```

![](plot1.png)

---

# References
