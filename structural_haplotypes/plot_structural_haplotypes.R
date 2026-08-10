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
if (length(args) != 4) {
    stop(
        "Usage: Rscript plot_structural_haplotypes.R ",
        "<representative_bundles.bed> <gene_track.bed> <tree.newick> <output.pdf>"
    )
}

bundle_file <- args[1]
gene_file <- args[2]
tree_file <- args[3]
output_pdf <- args[4]

hap_data <- read.table(
    bundle_file,
    header = FALSE,
    sep = "\t",
    col.names = c(
        "long_contig", "bstart", "bend", "bundle_info",
        "contig", "bundle_id", "bundle_type", "bundle_path"
    ),
    comment.char = ""
)

genes_data <- read.table(
    gene_file,
    header = FALSE,
    sep = "\t",
    col.names = c("long_contig", "tstart", "tend", "gene_name"),
    comment.char = ""
) %>%
    mutate(
        contig = str_split_fixed(long_contig, "::", 2)[, 1],
        contig = if_else(
            str_starts(contig, "apr"),
            str_replace(contig, "\\.hap([12])$", ".\\1"),
            contig
        ),
        strand = "+"
    )

tree <- read.tree(tree_file)
tree$tip.label <- ifelse(
    str_starts(tree$tip.label, "apr"),
    str_replace(tree$tip.label, "\\.hap([12])$", ".\\1"),
    tree$tip.label
)

haps_to_plot <- unique(hap_data$contig)
tree <- keep.tip(tree, tip = haps_to_plot)
plot_order <- tree$tip.label

hap_plot_data <- hap_data %>%
    mutate(
        rel_start = bstart,
        rel_end = bend,
        orientation = FALSE,
        contig = factor(contig, levels = rev(plot_order))
    ) %>%
    filter(!is.na(contig)) %>%
    arrange(contig)

genes_plot_data <- genes_data %>%
    mutate(contig = factor(contig, levels = rev(plot_order))) %>%
    filter(!is.na(contig)) %>%
    arrange(contig)

bundle_ids <- sort(unique(hap_plot_data$bundle_id))
bundle_colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(bundle_ids))
names(bundle_colors) <- bundle_ids

gene_colors <- c(
    "PGA34A" = "#4dbbd5",
    "PGA34B" = "#ff9900",
    "PGA5" = "#e64b35"
)

contig_labels <- hap_plot_data %>%
    group_by(contig) %>%
    summarise(label_x = max(rel_end, na.rm = TRUE), .groups = "drop")

p_tree <- ggtree(tree, size = 0.5) +
    theme(plot.margin = margin(0, 0, 0, 0))

p_structure <- ggplot(hap_plot_data, aes(y = contig)) +
    geom_gene_arrow(
        aes(
            xmin = rel_start,
            xmax = rel_end,
            fill = as.factor(bundle_id),
            forward = !orientation
        ),
        arrowhead_height = unit(3, "mm"),
        arrow_body_height = unit(2, "mm"),
        arrowhead_width = unit(3, "mm"),
        size = 0.6
    ) +
    geom_gene_arrow(
        data = genes_plot_data,
        aes(
            xmin = tstart,
            xmax = tend,
            color = gene_name,
            forward = (strand == "+")
        ),
        arrowhead_height = unit(1.5, "mm"),
        arrow_body_height = unit(1, "mm"),
        arrowhead_width = unit(3, "mm"),
        fill = "white",
        position = position_nudge(y = 0.45),
        size = 0.5
    ) +
    geom_text(
        data = genes_plot_data,
        aes(x = tstart, label = gene_name, color = gene_name),
        hjust = 1,
        size = 2.5,
        position = position_nudge(y = 0.45, x = -650)
    ) +
    geom_text(
        data = contig_labels,
        aes(x = label_x, label = contig),
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
        breaks = seq(
            0,
            max(hap_plot_data$rel_end, na.rm = TRUE) + 5000,
            by = 50000
        ),
        labels = function(x) x / 1000
    ) +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line.x = element_line(color = "black", size = 0.4),
        axis.ticks.x = element_line(color = "black", size = 0.4),
        plot.margin = margin(5.5, 50, 5.5, 0, "pt")
    )

final_plot <- p_structure %>% insert_left(p_tree, width = 0.1)

ggsave(
    output_pdf,
    final_plot,
    width = 12,
    height = 10
)
