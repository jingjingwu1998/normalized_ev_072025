#!/usr/bin/env Rscript
## Overlay merged SAME/DIFF regions with genes from TxDb (NO GTF)
## Outputs: overlay_results/genes_{TAG}_{same|diff}.csv

suppressPackageStartupMessages({
  library(GenomicFeatures)              # genes()
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)                 # Entrez -> Symbol
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
  library(AnnotationDbi)
  library(GenomeInfoDb)                 # seqlevelsStyle
  if (requireNamespace("rtracklayer", quietly = TRUE)) library(rtracklayer)
})

## ---------------------------
## USER SETTINGS
## ---------------------------
input_dir   <- "/Users/80030577/Desktop/HiC_analysis/EV_normalization/gene_annotation_ready_to_run_8.11.2025"
overlay_dir <- file.path(input_dir, "overlay_results")
dir.create(overlay_dir, showWarnings = FALSE, recursive = TRUE)

tags <- c("rbl_gcbc_vs_nbc_gcbc", "rbl_lcl_vs_nbc_lcl", "rbl_lcl_vs_rbl_nbc")

pad_kb        <- 10
upstream_kb   <- NULL
downstream_kb <- NULL
strand_aware  <- FALSE

## ---------------------------
## Helpers
## ---------------------------
coerce_to_gr <- function(x) {
  if (inherits(x, "GRanges"))     return(x)
  if (inherits(x, "GRangesList")) return(unlist(x, use.names = FALSE))
  if (is.list(x)) {
    parts <- lapply(x, coerce_to_gr)
    parts <- parts[vapply(parts, function(g) inherits(g, "GRanges") && length(g) > 0, logical(1))]
    if (length(parts) == 0L) return(GenomicRanges::GRanges())
    return(do.call(c, parts))
  }
  GenomicRanges::GRanges()
}

load_regions <- function(input_dir, kind = c("same","diff"), tag) {
  kind <- match.arg(kind)
  base <- file.path(input_dir, paste0("regions_", kind, "_", tag))
  paths <- c(paste0(base, ".rds"), paste0(base, ".bed"), paste0(base, ".csv"))
  exists_vec <- file.exists(paths)
  if (!any(exists_vec)) {
    message("[load_regions] Missing: ", base, "(.rds/.bed/.csv). Returning empty GRanges.")
    return(GRanges())
  }
  if (exists_vec[1]) {
    obj <- readRDS(paths[1]); return(coerce_to_gr(obj))
  } else if (exists_vec[2]) {
    if (!requireNamespace("rtracklayer", quietly = TRUE))
      stop("Found BED but rtracklayer not installed; install it or provide .rds/.csv.")
    return(coerce_to_gr(rtracklayer::import(paths[2], format = "BED")))
  } else {
    df <- read.csv(paths[3], stringsAsFactors = FALSE)
    req <- c("seqnames","start","end")
    if (!all(req %in% names(df))) stop("CSV missing seqnames/start/end: ", paths[3])
    return(GRanges(seqnames=df$seqnames, ranges=IRanges(start=as.integer(df$start), end=as.integer(df$end))))
  }
}

pad_genes <- function(genes_gr, pad_kb = 0, upstream_kb = NULL, downstream_kb = NULL, strand_aware = FALSE) {
  genes_gr <- coerce_to_gr(genes_gr)
  if (length(genes_gr) == 0L) return(genes_gr)
  if (!is.null(upstream_kb) || !is.null(downstream_kb)) {
    up <- if (is.null(upstream_kb)) 0L else as.integer(round(upstream_kb*1000))
    dn <- if (is.null(downstream_kb)) 0L else as.integer(round(downstream_kb*1000))
    if (strand_aware) {
      on_plus <- as.character(strand(genes_gr)) %in% c("+","*")
      start(genes_gr)[on_plus]  <- pmax(1L, start(genes_gr)[on_plus] - up)
      end(genes_gr)[on_plus]    <- end(genes_gr)[on_plus] + dn
      start(genes_gr)[!on_plus] <- pmax(1L, start(genes_gr)[!on_plus] - dn)
      end(genes_gr)[!on_plus]   <- end(genes_gr)[!on_plus] + up
    } else {
      start(genes_gr) <- pmax(1L, start(genes_gr) - up)
      end(genes_gr)   <- end(genes_gr) + dn
    }
  } else if (!is.null(pad_kb) && pad_kb > 0) {
    pad <- as.integer(round(pad_kb*1000))
    start(genes_gr) <- pmax(1L, start(genes_gr) - pad)
    end(genes_gr)   <- end(genes_gr) + pad
  }
  genes_gr
}

# Add gene symbols from Entrez IDs
add_symbols <- function(genes_gr) {
  ids <- as.character(mcols(genes_gr)$gene_id)
  syms <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = ids, keytype = "ENTREZID", column = "SYMBOL",
    multiVals = function(x) if (length(x)) x[1] else NA_character_
  )
  mcols(genes_gr)$gene_symbol <- unname(syms[ids])
  genes_gr
}

# Harmonize seqlevels style (e.g., "chr" vs no "chr")
match_seqstyle <- function(query, target) {
  try({
    style <- seqlevelsStyle(target)[1]
    seqlevelsStyle(query) <- style
  }, silent = TRUE)
  query
}

overlay_once <- function(regions, genes_pad, out_csv) {
  regions <- coerce_to_gr(regions)
  if (length(regions) == 0L) { write.csv(data.frame(), out_csv, row.names = FALSE); return(invisible(0L)) }
  
  # Harmonize chromosome naming styles
  regions <- match_seqstyle(regions, genes_pad)
  
  ov <- findOverlaps(regions, genes_pad, ignore.strand = TRUE)
  if (length(ov) == 0L) { write.csv(data.frame(), out_csv, row.names = FALSE); return(invisible(0L)) }
  
  ridx <- queryHits(ov); gidx <- subjectHits(ov)
  reg_df  <- as.data.frame(regions[ridx])
  gene_df <- as.data.frame(genes_pad[gidx])
  
  keep_reg  <- intersect(c("seqnames","start","end","width","sign",
                           "deltaEV1_mean","deltaEV1_min","deltaEV1_max",
                           "deltaEV2_mean","deltaEV2_min","deltaEV2_max","n_bins"),
                         colnames(reg_df))
  # TxDb knownGene has 'gene_id' (Entrez); we added 'gene_symbol'
  keep_gene <- intersect(c("seqnames","start","end","width","gene_id","gene_symbol"),
                         colnames(gene_df))
  
  out <- cbind(
    setNames(reg_df[keep_reg],   paste0("region_", keep_reg)),
    setNames(gene_df[keep_gene], paste0("gene_",   keep_gene))
  )
  write.csv(out, out_csv, row.names = FALSE)
  invisible(length(ov))
}

## ---------------------------
## Main (TxDb only)
## ---------------------------
cat("== Using TxDb.Hsapiens.UCSC.hg19.knownGene (no GTF) ==\n")
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

genes_gr  <- genes(txdb)        # GRanges with gene_id (Entrez)
genes_gr  <- add_symbols(genes_gr)
genes_pad <- pad_genes(genes_gr, pad_kb = pad_kb,
                       upstream_kb = upstream_kb, downstream_kb = downstream_kb,
                       strand_aware = strand_aware)

for (tag in tags) {
  cat(sprintf("\n== Overlay: %s ==\n", tag))
  same_path <- file.path(overlay_dir, paste0("genes_", tag, "_same.csv"))
  diff_path <- file.path(overlay_dir, paste0("genes_", tag, "_diff.csv"))
  
  regions_same <- load_regions(input_dir, "same", tag)
  regions_diff <- load_regions(input_dir, "diff", tag)
  
  cat(sprintf("  Regions: SAME=%d, DIFF=%d\n", length(regions_same), length(regions_diff)))
  n_same <- overlay_once(regions_same, genes_pad, same_path)
  n_diff <- overlay_once(regions_diff,  genes_pad, diff_path)
  cat(sprintf("  Overlaps written: SAME=%d -> %s\n", n_same, same_path))
  cat(sprintf("                     DIFF=%d -> %s\n", n_diff, diff_path))
}

cat("\n✅ Done. Overlay CSVs are in: ", overlay_dir, "\n")
