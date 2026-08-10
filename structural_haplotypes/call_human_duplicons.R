#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
})

option_list <- list(
  make_option("--manifest", type = "character", help = "TSV with haplotype, paf and repeatmasker columns."),
  make_option("--gene-track", type = "character", dest = "gene_track", help = "Representative PGA gene track BED4."),
  make_option("--output-prefix", type = "character", dest = "output_prefix", help = "Output prefix."),
  make_option("--min-alignment-length", type = "numeric", default = 2000, dest = "min_alignment_length"),
  make_option("--merge-distance", type = "numeric", default = 1000, dest = "merge_distance"),
  make_option("--candidate-merge-distance", type = "numeric", default = 1000, dest = "candidate_merge_distance"),
  make_option("--diagonal-buffer", type = "numeric", default = 5000, dest = "diagonal_buffer"),
  make_option("--min-gene-coverage", type = "numeric", default = 1.0, dest = "min_gene_coverage")
)
opt <- parse_args(OptionParser(option_list = option_list))

required <- c("manifest", "gene_track", "output_prefix")
missing <- required[vapply(required, function(x) is.null(opt[[x]]) || opt[[x]] == "", logical(1))]
if (length(missing) > 0) {
  stop("Missing required option(s): ", paste(missing, collapse = ", "), call. = FALSE)
}
if (opt$min_gene_coverage <= 0 || opt$min_gene_coverage > 1) {
  stop("--min-gene-coverage must be in (0, 1].", call. = FALSE)
}

normalize_haplotype <- function(x) {
  x <- as.character(x)
  x <- str_split_fixed(x, "::", 2)[, 1]
  if_else(str_starts(x, "apr"), str_replace(x, "\\.hap([12])$", ".\\1"), x)
}

read_paf_file <- function(paf_file, haplotype) {
  paf <- suppressWarnings(read.table(
    paf_file, sep = "\t", header = FALSE, fill = TRUE, quote = "",
    comment.char = "", stringsAsFactors = FALSE
  ))
  if (ncol(paf) < 12 || nrow(paf) == 0) return(tibble())

  paf <- paf[, 1:12]
  colnames(paf) <- c(
    "queryID", "queryLen", "queryStart", "queryEnd", "strand",
    "refID", "refLen", "refStart", "refEnd", "numResidueMatches",
    "lenAln", "mapQ"
  )

  paf <- paf %>%
    mutate(
      across(c(queryLen, queryStart, queryEnd, refLen, refStart, refEnd,
               numResidueMatches, lenAln, mapQ), as.numeric),
      haplotype = haplotype,
      paf_row = row_number()
    ) %>%
    filter(queryID == refID)

  neg <- which(paf$strand == "-")
  if (length(neg) > 0) {
    tmp <- paf$queryStart[neg]
    paf$queryStart[neg] <- paf$queryEnd[neg]
    paf$queryEnd[neg] <- tmp
  }

  paf %>%
    mutate(
      qmin = pmin(queryStart, queryEnd),
      qmax = pmax(queryStart, queryEnd),
      rmin = pmin(refStart, refEnd),
      rmax = pmax(refStart, refEnd),
      identity = numResidueMatches / lenAln,
      orig_ids = as.character(paf_row)
    )
}

merge_collinear_blocks <- function(aln, merge_dist) {
  if (nrow(aln) == 0) return(aln)
  aln <- as.data.frame(aln, stringsAsFactors = FALSE)
  changed <- TRUE

  while (changed) {
    changed <- FALSE
    merged_list <- list()
    used <- rep(FALSE, nrow(aln))

    for (i in seq_len(nrow(aln))) {
      if (used[i]) next
      curr <- aln[i, ]
      used[i] <- TRUE
      curr_offset <- ifelse(
        curr$strand == "+",
        curr$refStart - curr$queryStart,
        curr$refStart + curr$queryStart
      )

      if (i < nrow(aln)) {
        for (j in (i + 1):nrow(aln)) {
          if (used[j]) next
          cand <- aln[j, ]
          same_axis <- curr$refID == cand$refID && curr$queryID == cand$queryID && curr$strand == cand$strand
          if (!same_axis) next

          cand_offset <- ifelse(
            cand$strand == "+",
            cand$refStart - cand$queryStart,
            cand$refStart + cand$queryStart
          )
          if (abs(curr_offset - cand_offset) > merge_dist) next

          r_dist <- max(curr$refStart, cand$refStart) - min(curr$refEnd, cand$refEnd)
          q_min_curr <- min(curr$queryStart, curr$queryEnd)
          q_max_curr <- max(curr$queryStart, curr$queryEnd)
          q_min_cand <- min(cand$queryStart, cand$queryEnd)
          q_max_cand <- max(cand$queryStart, cand$queryEnd)
          q_dist <- max(q_min_curr, q_min_cand) - min(q_max_curr, q_max_cand)

          dist_ok <- r_dist <= merge_dist && r_dist >= -merge_dist && q_dist <= merge_dist && q_dist >= -merge_dist
          cand_encl <- cand$refStart >= curr$refStart && cand$refEnd <= curr$refEnd &&
            q_min_cand >= q_min_curr && q_max_cand <= q_max_curr
          curr_encl <- curr$refStart >= cand$refStart && curr$refEnd <= cand$refEnd &&
            q_min_curr >= q_min_cand && q_max_curr <= q_max_cand

          if (dist_ok && !cand_encl && !curr_encl) {
            curr$refStart <- min(curr$refStart, cand$refStart)
            curr$refEnd <- max(curr$refEnd, cand$refEnd)
            if (curr$strand == "+") {
              curr$queryStart <- min(curr$queryStart, cand$queryStart)
              curr$queryEnd <- max(curr$queryEnd, cand$queryEnd)
            } else {
              curr$queryStart <- max(curr$queryStart, cand$queryStart)
              curr$queryEnd <- min(curr$queryEnd, cand$queryEnd)
            }
            curr$numResidueMatches <- curr$numResidueMatches + cand$numResidueMatches
            curr$lenAln <- curr$lenAln + cand$lenAln
            curr$orig_ids <- paste(curr$orig_ids, cand$orig_ids, sep = ",")
            used[j] <- TRUE
            changed <- TRUE
          }
        }
      }
      merged_list[[length(merged_list) + 1]] <- curr
    }
    aln <- bind_rows(merged_list)
  }

  aln %>%
    mutate(
      qmin = pmin(queryStart, queryEnd),
      qmax = pmax(queryStart, queryEnd),
      rmin = pmin(refStart, refEnd),
      rmax = pmax(refStart, refEnd),
      identity = numResidueMatches / lenAln
    )
}

build_shadow_candidates <- function(paf, haplotype, min_len, diagonal_buffer, merge_dist) {
  if (nrow(paf) == 0) return(tibble())

  merged <- merge_collinear_blocks(paf, merge_dist) %>%
    mutate(
      diagonal_offset = abs(refStart - queryStart),
      macro_block_id = paste0(haplotype, ":MB", row_number())
    ) %>%
    filter(lenAln > min_len, diagonal_offset >= diagonal_buffer)

  if (nrow(merged) == 0) return(tibble())

  x_shadow <- merged %>%
    transmute(
      haplotype,
      x_block_id = macro_block_id,
      x_start = rmin,
      x_end = rmax,
      x_identity = identity,
      x_paf_rows = orig_ids
    )
  y_shadow <- merged %>%
    transmute(
      y_block_id = macro_block_id,
      y_start = qmin,
      y_end = qmax,
      y_identity = identity,
      y_paf_rows = orig_ids
    )

  tidyr::expand_grid(x_shadow, y_shadow) %>%
    mutate(
      candidate_start = pmax(x_start, y_start),
      candidate_end = pmin(x_end, y_end),
      candidate_length = candidate_end - candidate_start,
      mean_identity = (x_identity + y_identity) / 2
    ) %>%
    filter(candidate_length > 0)
}

read_repeatmasker <- function(path, haplotype) {
  rm <- tryCatch(
    suppressWarnings(read.table(path, skip = 3, fill = TRUE, stringsAsFactors = FALSE, quote = "", comment.char = "")),
    error = function(e) NULL
  )
  if (is.null(rm) || nrow(rm) == 0 || ncol(rm) < 11) return(tibble())

  rm %>%
    transmute(
      haplotype = haplotype,
      seqid = V5,
      te_start = pmin(as.numeric(V6), as.numeric(V7)),
      te_end = pmax(as.numeric(V6), as.numeric(V7)),
      strand = V9,
      repeat_name = V10,
      class_family = V11
    ) %>%
    filter(!is.na(te_start), !is.na(te_end))
}

prepare_genes <- function(gene_track) {
  read.table(
    gene_track, sep = "\t", header = FALSE, quote = "", comment.char = "",
    col.names = c("contig", "gene_start", "gene_end", "gene_name"),
    stringsAsFactors = FALSE
  ) %>%
    mutate(haplotype = normalize_haplotype(contig)) %>%
    arrange(haplotype, gene_start, gene_end) %>%
    group_by(haplotype) %>%
    mutate(
      gene_index = row_number(),
      prev_gene_end = lag(gene_end),
      next_gene_start = lead(gene_start),
      gene_id = paste(haplotype, gene_index, gene_name, sep = "|")
    ) %>%
    ungroup()
}

make_cluster_rows <- function(cands, merge_dist, min_cov) {
  if (nrow(cands) == 0) return(tibble())
  sorted <- cands %>% arrange(candidate_start, candidate_end, candidate_rank)
  gene_len <- sorted$gene_end[1] - sorted$gene_start[1]
  rows <- list()

  emit <- function(base, idx, start, end) {
    overlap <- pmax(0, pmin(end, base$gene_end[1]) - pmax(start, base$gene_start[1]))
    rows[[length(rows) + 1]] <<- tibble(
      haplotype = base$haplotype[1],
      gene_id = base$gene_id[1],
      gene_name = base$gene_name[1],
      gene_index = base$gene_index[1],
      gene_start = base$gene_start[1],
      gene_end = base$gene_end[1],
      prev_gene_end = base$prev_gene_end[1],
      next_gene_start = base$next_gene_start[1],
      cluster_start = start,
      cluster_end = end,
      cluster_length = end - start,
      cluster_coverage = if (gene_len > 0) overlap / gene_len else 0,
      cluster_overlaps_adjacent_gene =
        (!is.na(base$prev_gene_end[1]) && start < base$prev_gene_end[1]) ||
        (!is.na(base$next_gene_start[1]) && end > base$next_gene_start[1]),
      mean_identity = mean(base$mean_identity[idx], na.rm = TRUE),
      source_candidate_count = length(idx),
      source_candidate_ids = paste(base$primitive_candidate_id[idx], collapse = ","),
      x_block_id = paste(unique(base$x_block_id[idx]), collapse = ";"),
      y_block_id = paste(unique(base$y_block_id[idx]), collapse = ";"),
      x_paf_rows = paste(unique(unlist(str_split(base$x_paf_rows[idx], ","))), collapse = ","),
      y_paf_rows = paste(unique(unlist(str_split(base$y_paf_rows[idx], ","))), collapse = ",")
    )
  }

  for (idx in seq_len(nrow(sorted))) {
    emit(sorted, idx, sorted$candidate_start[idx], sorted$candidate_end[idx])
  }

  emit_short_merges <- function(base) {
    if (nrow(base) < 2) return(invisible(NULL))
    for (start_idx in seq_len(nrow(base))) {
      start <- base$candidate_start[start_idx]
      end <- base$candidate_end[start_idx]
      idx <- start_idx
      singleton_cov <- if (gene_len > 0) pmax(0, pmin(end, base$gene_end[1]) - pmax(start, base$gene_start[1])) / gene_len else 0
      if (singleton_cov + 1e-9 >= min_cov || start_idx == nrow(base)) next

      for (end_idx in (start_idx + 1):nrow(base)) {
        if (base$candidate_start[end_idx] > end + merge_dist) break
        start <- min(start, base$candidate_start[end_idx])
        end <- max(end, base$candidate_end[end_idx])
        idx <- c(idx, end_idx)
        cov <- if (gene_len > 0) pmax(0, pmin(end, base$gene_end[1]) - pmax(start, base$gene_start[1])) / gene_len else 0
        if (cov + 1e-9 >= min_cov) {
          emit(base, idx, start, end)
          break
        }
      }
    }
  }

  emit_short_merges(sorted)

  redundant_same_start <- vapply(seq_len(nrow(sorted)), function(i) {
    any(sorted$candidate_start == sorted$candidate_start[i] & sorted$candidate_end < sorted$candidate_end[i])
  }, logical(1))
  redundant_same_end <- vapply(seq_len(nrow(sorted)), function(i) {
    any(sorted$candidate_end == sorted$candidate_end[i] & sorted$candidate_start > sorted$candidate_start[i])
  }, logical(1))
  reduced <- sorted[!(redundant_same_start | redundant_same_end), , drop = FALSE]
  emit_short_merges(reduced)

  bind_rows(rows) %>%
    distinct(gene_id, cluster_start, cluster_end, source_candidate_ids, .keep_all = TRUE) %>%
    mutate(cluster_start_distance = abs(cluster_start - gene_start)) %>%
    arrange(
      cluster_length,
      desc(cluster_coverage),
      desc(mean_identity),
      cluster_start_distance,
      source_candidate_count,
      cluster_start,
      cluster_end
    ) %>%
    mutate(cluster_rank = row_number())
}

call_duplicated_cores <- function(genes, shadows, min_cov, candidate_merge_dist) {
  candidates <- bind_rows(lapply(seq_len(nrow(genes)), function(i) {
    gene <- genes[i, ]
    x <- shadows %>% filter(haplotype == gene$haplotype)
    if (nrow(x) == 0) return(NULL)

    gene_len <- gene$gene_end - gene$gene_start
    x %>%
      mutate(
        gene_id = gene$gene_id,
        gene_name = gene$gene_name,
        gene_index = gene$gene_index,
        gene_start = gene$gene_start,
        gene_end = gene$gene_end,
        prev_gene_end = gene$prev_gene_end,
        next_gene_start = gene$next_gene_start,
        gene_overlap = pmax(0, pmin(candidate_end, gene$gene_end) - pmax(candidate_start, gene$gene_start)),
        gene_coverage = if (gene_len > 0) gene_overlap / gene_len else 0,
        start_distance = abs(candidate_start - gene$gene_start)
      ) %>%
      filter(gene_overlap > 0) %>%
      arrange(candidate_length, desc(gene_coverage), desc(mean_identity), start_distance, candidate_start, candidate_end) %>%
      mutate(
        candidate_rank = row_number(),
        primitive_candidate_id = paste0("PC", candidate_rank)
      )
  }))

  if (nrow(candidates) == 0) return(tibble())

  clusters <- candidates %>%
    group_split(gene_id, .keep = TRUE) %>%
    lapply(make_cluster_rows, merge_dist = candidate_merge_dist, min_cov = min_cov) %>%
    bind_rows()

  clusters %>%
    filter(cluster_coverage + 1e-9 >= min_cov, !cluster_overlaps_adjacent_gene) %>%
    group_by(haplotype, gene_id, gene_name, gene_index, gene_start, gene_end) %>%
    slice_min(order_by = cluster_rank, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(haplotype, gene_start, cluster_start, cluster_end) %>%
    group_by(haplotype) %>%
    mutate(core_id = paste0("CORE", row_number())) %>%
    ungroup() %>%
    transmute(
      haplotype,
      gene_id,
      gene_name,
      gene_index,
      gene_start,
      gene_end,
      core_id,
      core_start = cluster_start,
      core_end = cluster_end,
      core_length = cluster_length,
      gene_coverage = cluster_coverage,
      mean_identity,
      source_candidate_count,
      source_candidate_ids,
      x_block_id,
      y_block_id,
      x_paf_rows,
      y_paf_rows
    )
}

build_duplicons <- function(cores) {
  if (nrow(cores) == 0) return(cores)
  cores %>%
    arrange(haplotype, core_start, core_end) %>%
    group_by(haplotype) %>%
    mutate(
      next_core_id = lead(core_id),
      next_core_start = lead(core_start),
      has_internal_spacer = !is.na(next_core_start) & next_core_start > core_end,
      duplicon_start = core_start,
      duplicon_end = if_else(has_internal_spacer, next_core_start, core_end),
      duplicon_length = duplicon_end - duplicon_start,
      spacer_start = if_else(has_internal_spacer, core_end, NA_real_),
      spacer_end = if_else(has_internal_spacer, next_core_start, NA_real_),
      spacer_length = if_else(has_internal_spacer, next_core_start - core_end, 0)
    ) %>%
    ungroup() %>%
    select(-next_core_start, -has_internal_spacer)
}

endpoint_repeats <- function(duplicons, rm_all) {
  if (nrow(duplicons) == 0 || nrow(rm_all) == 0) return(tibble())

  bind_rows(lapply(seq_len(nrow(duplicons)), function(i) {
    d <- duplicons[i, ]
    rm <- rm_all %>% filter(haplotype == d$haplotype)
    if (nrow(rm) == 0) return(NULL)

    endpoints <- tibble(
      endpoint_type = c("duplicon_start", "duplicon_end"),
      endpoint_pos = c(d$duplicon_start, d$duplicon_end)
    )

    bind_rows(lapply(seq_len(nrow(endpoints)), function(j) {
      e <- endpoints[j, ]
      rm %>%
        filter(te_end >= e$endpoint_pos, te_start <= e$endpoint_pos) %>%
        mutate(
          gene_id = d$gene_id,
          gene_name = d$gene_name,
          core_id = d$core_id,
          duplicon_start = d$duplicon_start,
          duplicon_end = d$duplicon_end,
          endpoint_type = e$endpoint_type,
          endpoint_pos = e$endpoint_pos
        )
    }))
  })) %>%
    arrange(haplotype, duplicon_start, endpoint_pos, te_start, te_end)
}

manifest <- read_tsv(opt$manifest, show_col_types = FALSE, col_types = cols(.default = col_character()))
required_manifest <- c("haplotype", "paf", "repeatmasker")
if (!all(required_manifest %in% colnames(manifest))) {
  stop("Manifest must contain: ", paste(required_manifest, collapse = ", "), call. = FALSE)
}
manifest <- manifest %>% mutate(haplotype = normalize_haplotype(haplotype))

for (f in c(opt$gene_track, manifest$paf, manifest$repeatmasker)) {
  if (!file.exists(f)) stop("Input file not found: ", f, call. = FALSE)
}

genes <- prepare_genes(opt$gene_track) %>% filter(haplotype %in% manifest$haplotype)

shadows <- bind_rows(lapply(seq_len(nrow(manifest)), function(i) {
  row <- manifest[i, ]
  paf <- read_paf_file(row$paf, row$haplotype)
  build_shadow_candidates(
    paf,
    row$haplotype,
    min_len = opt$min_alignment_length,
    diagonal_buffer = opt$diagonal_buffer,
    merge_dist = opt$merge_distance
  )
}))

cores <- call_duplicated_cores(
  genes,
  shadows,
  min_cov = opt$min_gene_coverage,
  candidate_merge_dist = opt$candidate_merge_distance
)
duplicons <- build_duplicons(cores)

rm_all <- bind_rows(lapply(seq_len(nrow(manifest)), function(i) {
  read_repeatmasker(manifest$repeatmasker[i], manifest$haplotype[i])
}))
te_hits <- endpoint_repeats(duplicons, rm_all)

core_out <- paste0(opt$output_prefix, ".duplicons.tsv")
te_out <- paste0(opt$output_prefix, ".duplicon_endpoint_repeats.tsv")
write_tsv(duplicons, core_out, na = ".")
write_tsv(te_hits, te_out, na = ".")

message("Wrote: ", core_out)
message("Wrote: ", te_out)
