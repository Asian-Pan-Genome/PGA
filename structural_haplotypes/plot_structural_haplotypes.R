#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggtree)
  library(gggenes)
  library(aplot)
  library(RColorBrewer)
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
  stop(
    "Usage: Rscript plot_structural_haplotypes.R ",
    "<representative_bundles.bed> <gene_track.bed> <duplicons.tsv> ",
    "<duplicon_endpoint_repeats.tsv> <tree.newick> <output.pdf>"
  )
}

bundle_file <- args[1]
gene_file <- args[2]
duplicon_file <- args[3]
te_file <- args[4]
tree_file <- args[5]
output_pdf <- args[6]

normalize_haplotype <- function(x) {
  x <- as.character(x)
  x <- str_split_fixed(x, "::", 2)[, 1]
  if_else(str_starts(x, "apr"), str_replace(x, "\\.hap([12])$", ".\\1"), x)
}

hap_data <- read.table(
  bundle_file,
  header = FALSE,
  sep = "\t",
  col.names = c(
    "long_contig", "bstart", "bend", "bundle_info",
    "contig", "bundle_id", "bundle_type", "bundle_path"
  ),
  comment.char = "",
  stringsAsFactors = FALSE
) %>%
  mutate(contig = normalize_haplotype(contig))

genes_data <- read.table(
  gene_file,
  header = FALSE,
  sep = "\t",
  col.names = c("long_contig", "tstart", "tend", "gene_name"),
  comment.char = "",
  stringsAsFactors = FALSE
) %>%
  mutate(contig = normalize_haplotype(long_contig), strand = "+")

duplicons <- read_tsv(duplicon_file, show_col_types = FALSE) %>%
  mutate(haplotype = normalize_haplotype(haplotype))

te_hits <- read_tsv(te_file, show_col_types = FALSE) %>%
  mutate(haplotype = normalize_haplotype(haplotype))

tree <- read.tree(tree_file)
tree$tip.label <- normalize_haplotype(tree$tip.label)
haps_to_plot <- unique(hap_data$contig)
tree <- keep.tip(tree, tip = haps_to_plot)

p_tree <- ggtree(tree, size = 0.5) +
  theme(plot.margin = margin(0, 0, 0, 0))

# Use the actual ggtree tip order rather than the raw Newick label order.
plot_order <- p_tree$data %>%
  filter(isTip) %>%
  arrange(desc(y)) %>%
  pull(label)

y_levels <- rev(plot_order)
y_lookup <- setNames(seq_along(y_levels), y_levels)

hap_plot_data <- hap_data %>%
  filter(contig %in% plot_order) %>%
  mutate(
    rel_start = bstart,
    rel_end = bend,
    y = y_lookup[contig]
  )

genes_plot_data <- genes_data %>%
  filter(contig %in% plot_order) %>%
  mutate(y = y_lookup[contig] + 0.45)

duplicon_plot_data <- duplicons %>%
  filter(haplotype %in% plot_order) %>%
  mutate(y = y_lookup[haplotype] + 0.22)

te_plot_data <- te_hits %>%
  filter(haplotype %in% plot_order) %>%
  mutate(
    y = y_lookup[haplotype] + 0.22,
    te_mid = (te_start + te_end) / 2
  )

bundle_ids <- sort(unique(hap_plot_data$bundle_id))
bundle_colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(bundle_ids))
names(bundle_colors) <- bundle_ids

gene_colors <- c(
  "PGA34A" = "#4dbbd5",
  "PGA34B" = "#ff9900",
  "PGA5" = "#e64b35"
)

te_colors <- c(
  "SINE/Alu" = "#e41a1c",
  "LTR/ERV" = "#377eb8",
  "LINE/L1" = "#4daf4a",
  "Retroposon/SVA" = "#984ea3"
)

te_plot_data <- te_plot_data %>%
  mutate(
    repeat_group = case_when(
      str_detect(class_family, "SINE/Alu") ~ "SINE/Alu",
      str_detect(class_family, "LTR/ERV") ~ "LTR/ERV",
      str_detect(class_family, "LINE/L1") ~ "LINE/L1",
      str_detect(class_family, "Retroposon/SVA") ~ "Retroposon/SVA",
      TRUE ~ "Other"
    ),
    te_color = if_else(
      repeat_group %in% names(te_colors),
      unname(te_colors[repeat_group]),
      "#999999"
    )
  )

contig_labels <- hap_plot_data %>%
  group_by(contig, y) %>%
  summarise(label_x = max(rel_end, na.rm = TRUE), .groups = "drop")

p_structure <- ggplot() +
  geom_gene_arrow(
    data = hap_plot_data,
    aes(
      xmin = rel_start,
      xmax = rel_end,
      y = y,
      fill = as.factor(bundle_id),
      forward = TRUE
    ),
    arrowhead_height = unit(3, "mm"),
    arrow_body_height = unit(2, "mm"),
    arrowhead_width = unit(3, "mm"),
    size = 0.4
  ) +
  geom_rect(
    data = duplicon_plot_data,
    aes(
      xmin = duplicon_start,
      xmax = duplicon_end,
      ymin = y - 0.07,
      ymax = y + 0.07
    ),
    fill = "grey65",
    color = NA,
    alpha = 0.65,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = te_plot_data,
    aes(
      x = te_mid,
      xend = te_mid,
      y = y - 0.10,
      yend = y + 0.10,
      color = I(te_color)
    ),
    linewidth = 0.45,
    inherit.aes = FALSE
  ) +
  geom_gene_arrow(
    data = genes_plot_data,
    aes(
      xmin = tstart,
      xmax = tend,
      y = y,
      color = gene_name,
      forward = (strand == "+")
    ),
    arrowhead_height = unit(1.5, "mm"),
    arrow_body_height = unit(1, "mm"),
    arrowhead_width = unit(3, "mm"),
    fill = "white",
    size = 0.5
  ) +
  geom_text(
    data = genes_plot_data,
    aes(x = tstart, y = y, label = gene_name, color = gene_name),
    hjust = 1,
    size = 2.5,
    position = position_nudge(x = -650)
  ) +
  geom_text(
    data = contig_labels,
    aes(x = label_x, y = y, label = contig),
    hjust = 0,
    size = 4.5,
    color = "grey40",
    position = position_nudge(x = 1000)
  ) +
  scale_fill_manual(values = bundle_colors) +
  scale_color_manual(values = gene_colors) +
  scale_x_continuous(
    "Size (kb)",
    expand = c(0.01, 0.01),
    breaks = seq(0, max(hap_plot_data$rel_end, na.rm = TRUE) + 5000, by = 50000),
    labels = function(x) x / 1000
  ) +
  scale_y_continuous(
    limits = c(0.5, length(y_levels) + 0.8),
    expand = c(0, 0)
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.4),
    axis.ticks.x = element_line(color = "black", linewidth = 0.4),
    plot.margin = margin(5.5, 50, 5.5, 0, "pt")
  ) +
  coord_cartesian(clip = "off")

final_plot <- p_structure %>% insert_left(p_tree, width = 0.1)

ggsave(output_pdf, final_plot, width = 12, height = 10)
