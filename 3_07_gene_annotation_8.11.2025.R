
## ====== Libraries ======
library(GenomicRanges)
library(GenomeInfoDb)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

## ===== 1) Build promoter set (hg19) =====
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genetxdb  = genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
mcols(genetxdb)$gene_id <- names(genetxdb)
genome(genetxdb) <- "hg19"

U <- 2500L; D <- 500L   # promoter window: -2500/+500
# U <- 0L; D <- 1L   # promoter window: -2500/+500

prom_gene <- promoters(genetxdb, upstream = U, downstream = D)
mcols(prom_gene)$gene_id <- mcols(genetxdb)$gene_id
prom_gene <- trim(prom_gene)
genome(prom_gene) <- "hg19"

## ===== 2) readers/helpers =====

# find first non-empty (non-whitespace) line
.first_non_empty <- function(path){
  lines <- readLines(path, warn = FALSE)
  nz <- which(nzchar(trimws(lines)))
  if (!length(nz)) stop("Empty file: ", path)
  list(lines = lines, first = nz[1])
}

# Reader: skips leading blanks, auto-detects delimiter, handles quoted headers
smart_read_regions <- function(path){
  info <- .first_non_empty(path)
  sep <- if (grepl("\t", info$lines[info$first])) "\t" else ","
  df <- read.table(path,
                   header = TRUE,
                   sep = sep,
                   quote = "\"'",
                   comment.char = "",
                   fill = TRUE,
                   check.names = FALSE,
                   strip.white = TRUE,
                   stringsAsFactors = FALSE,
                   skip = info$first - 1L)
  # sanitize column names (strip quotes/space/punct; lowercase)
  norm <- function(x){
    x <- gsub("^\\s+|\\s+$", "", x)
    x <- gsub("^\"|\"$", "", x); x <- gsub("^'|'$", "", x); x <- gsub("`", "", x)
    x <- gsub("[[:space:]]+", "", x)
    tolower(x)
  }
  names(df) <- norm(names(df))
  df
}

# Convert to GRanges; attach hg19 lengths, then trim BEFORE validation
to_regions_gr <- function(csv_path, template_seqinfo){
  df <- smart_read_regions(csv_path)
  
  # accept your names or common aliases
  pick <- function(cands){
    ix <- match(cands, names(df))
    ix <- ix[!is.na(ix)]
    if (length(ix)) ix[1] else NA_integer_
  }
  i_chr   <- pick(c("region_seqnames","chr","chrom","seqnames","chromosome"))
  i_start <- pick(c("region_start","start","chromstart","start_bp","startpos"))
  i_end   <- pick(c("region_end","end","chromend","end_bp","endpos"))
  
  if (any(is.na(c(i_chr, i_start, i_end)))) {
    stop("Could not find chr/start/end in: ", csv_path,
         "\nColumns seen: ", paste(names(df), collapse=", "))
  }
  
  chr <- as.character(df[[i_chr]])
  st  <- suppressWarnings(as.integer(df[[i_start]]))
  en  <- suppressWarnings(as.integer(df[[i_end]]))
  
  keep <- !is.na(chr) & !is.na(st) & !is.na(en)
  if (!all(keep)) {
    warning("Dropping ", sum(!keep), " rows with missing chr/start/end in ", basename(csv_path))
    chr <- chr[keep]; st <- st[keep]; en <- en[keep]
  }
  
  gr <- GRanges(seqnames = chr, ranges = IRanges(start = st, end = en))
  genome(gr) <- "hg19"
  
  # harmonize styles, restrict to common seqlevels
  seqlevelsStyle(gr) <- seqlevelsStyle(template_seqinfo)[1]
  common <- intersect(seqlevels(gr), seqlevels(template_seqinfo))
  gr <- keepSeqlevels(gr, common, pruning.mode = "coarse")
  
  # attach lengths, then TRIM BEFORE validity check to avoid warnings
  suppressWarnings(seqlengths(gr) <- seqlengths(template_seqinfo)[seqlevels(gr)])
  gr <- trim(gr)  # clip to chromosome bounds
  
  gr
}


# Overlap (bins unstranded) with promoter set; tidy output
overlap_regions_promoters <- function(regions_gr, promoters_gr){
  seqlevelsStyle(promoters_gr) <- seqlevelsStyle(regions_gr)[1]
  common <- intersect(seqlevels(regions_gr), seqlevels(promoters_gr))
  regions_gr   <- keepSeqlevels(regions_gr,  common, pruning.mode = "coarse")
  promoters_gr <- keepSeqlevels(promoters_gr, common, pruning.mode = "coarse")
  
  hits <- findOverlaps(regions_gr, promoters_gr, ignore.strand = TRUE)
  # number of region–promoter pairs
  if (!length(hits)) return(data.frame())
  
  ov <- pintersect(regions_gr[queryHits(hits)], promoters_gr[subjectHits(hits)])
  data.frame(
    region_chr   = as.character(seqnames(regions_gr))[queryHits(hits)],
    region_start = start(regions_gr)[queryHits(hits)],
    region_end   = end(regions_gr)[queryHits(hits)],
    promoter_chr = as.character(seqnames(promoters_gr))[subjectHits(hits)],
    promoter_start = start(promoters_gr)[subjectHits(hits)],
    promoter_end   = end(promoters_gr)[subjectHits(hits)],
    gene_id     = mcols(promoters_gr)$gene_id[subjectHits(hits)],
    overlap_bp  = width(ov),
    frac_region_covered   = width(ov) / width(regions_gr)[queryHits(hits)],
    frac_promoter_covered = width(ov) / width(promoters_gr)[subjectHits(hits)]
  )
}

## ===== 3) Your inputs & run =====
# change file names
input_files <- c(
  "/Users/80030577/Desktop/HiC_analysis/EV_normalization/gene_annotation_ready_to_run_8.11.2025/regions_same_rbl_lcl_vs_nbc_lcl.csv"
)

for (f in input_files) {
  message("Processing: ", f)
  regions_gr <- to_regions_gr(f, template_seqinfo = seqinfo(prom_gene))
  res <- overlap_regions_promoters(regions_gr, prom_gene)
  out_file <- sub("\\.csv$", "_overlap_promoters.csv", f, ignore.case = TRUE)
  write.csv(res, out_file, row.names = FALSE)
  message("  Saved: ", out_file, " (", nrow(res), " overlaps)")
}

message("Done! Overlap files written for all inputs.")
