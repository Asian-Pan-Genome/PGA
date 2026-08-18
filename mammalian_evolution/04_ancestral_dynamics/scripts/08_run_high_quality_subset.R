#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 7) {
  stop("
Usage:
Rscript 08_run_high_quality_subset.R tree.nwk PGA_CN.tsv diet.tsv quality.tsv out_prefix min_contig_N50 require_gap_free

Example:
Rscript 08_run_high_quality_subset.R \\
  mammal.tree.nwk \\
  PGA_CN.tsv \\
  diet.tsv \\
  quality.tsv \\
  301.PGA.analysis5.N50_10Mb \\
  10000000 \\
  TRUE
")
}

tree_file <- args[1]
cn_file <- args[2]
diet_file <- args[3]
quality_file <- args[4]
out_prefix <- args[5]
min_n50 <- as.numeric(args[6])
require_gap_free <- as.logical(args[7])

message("[1] Reading files...")

tr <- read.tree(tree_file)

cn <- read.table(
  cn_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

diet <- read.table(
  diet_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

qual <- read.table(
  quality_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

if (!all(c("species", "PGA_CN") %in% colnames(cn))) {
  stop("PGA_CN.tsv must contain columns: species, PGA_CN")
}

if (!all(c("species", "diet_state") %in% colnames(diet))) {
  stop("diet.tsv must contain columns: species, diet_state")
}

if (!"species" %in% colnames(qual)) {
  stop("quality.tsv must contain column: species")
}

if (!"contig_N50" %in% colnames(qual)) {
  stop("quality.tsv must contain column: contig_N50")
}

cn$species <- gsub(" ", "_", cn$species)
diet$species <- gsub(" ", "_", diet$species)
qual$species <- gsub(" ", "_", qual$species)

message("[2] Filtering high-quality species...")

qual$pass_N50 <- qual$contig_N50 >= min_n50

if ("gap_free" %in% colnames(qual) && require_gap_free) {
  qual$gap_free_norm <- tolower(as.character(qual$gap_free))
  qual$pass_gap <- qual$gap_free_norm %in% c("true", "t", "yes", "y", "1", "pass", "gap_free")
} else {
  qual$pass_gap <- TRUE
}

if ("quality_pass" %in% colnames(qual)) {
  qual$quality_pass_norm <- tolower(as.character(qual$quality_pass))
  qual$pass_quality_flag <- qual$quality_pass_norm %in% c("true", "t", "yes", "y", "1", "pass")
} else {
  qual$pass_quality_flag <- TRUE
}

qual$pass_high_quality <- qual$pass_N50 & qual$pass_gap & qual$pass_quality_flag

hq_species <- qual$species[qual$pass_high_quality]

message("High-quality species: ", length(hq_species))

dat <- merge(cn, diet, by = "species")
dat <- merge(dat, qual, by = "species", all.x = TRUE)

dat <- dat[dat$species %in% hq_species, ]
dat <- dat[!is.na(dat$PGA_CN) & !is.na(dat$diet_state), ]

common_tips <- intersect(tr$tip.label, dat$species)

if (length(common_tips) < 20) {
  stop("Too few high-quality matched species: ", length(common_tips))
}

message("High-quality matched tips in tree: ", length(common_tips))

tr <- drop.tip(tr, setdiff(tr$tip.label, common_tips))
dat <- dat[match(tr$tip.label, dat$species), ]

if (any(dat$species != tr$tip.label)) {
  stop("Species order mismatch after pruning.")
}

if (is.null(tr$edge.length)) {
  stop("Tree must have branch lengths.")
}

min_positive_bl <- min(tr$edge.length[tr$edge.length > 0], na.rm = TRUE)
tr$edge.length[tr$edge.length <= 0] <- min_positive_bl / 10

Ntip <- length(tr$tip.label)
Nnode <- tr$Nnode
internal_nodes <- (Ntip + 1):(Ntip + Nnode)

message("[3] Reconstructing ancestral PGA CN...")

cn_vec <- dat$PGA_CN
names(cn_vec) <- dat$species

cn_ace <- ace(
  x = cn_vec,
  phy = tr,
  type = "continuous",
  method = "REML"
)

node_CN <- rep(NA_real_, Ntip + Nnode)
names(node_CN) <- as.character(1:(Ntip + Nnode))

node_CN[as.character(1:Ntip)] <- cn_vec[tr$tip.label]

if (is.null(names(cn_ace$ace))) {
  names(cn_ace$ace) <- internal_nodes
}

node_CN[as.character(names(cn_ace$ace))] <- cn_ace$ace

node_CN_round <- round(node_CN)
node_CN_round[node_CN_round < 0] <- 0

message("[4] Reconstructing ancestral diet state...")

diet_vec <- as.factor(dat$diet_state)
names(diet_vec) <- dat$species

diet_ace <- ace(
  x = diet_vec,
  phy = tr,
  type = "discrete",
  model = "ER"
)

node_diet <- rep(NA_character_, Ntip + Nnode)
names(node_diet) <- as.character(1:(Ntip + Nnode))

node_diet[as.character(1:Ntip)] <- as.character(diet_vec[tr$tip.label])

lik <- diet_ace$lik.anc

if (is.null(rownames(lik))) {
  rownames(lik) <- internal_nodes
}

internal_diet <- apply(lik, 1, function(x) {
  colnames(lik)[which.max(x)]
})

node_diet[as.character(rownames(lik))] <- internal_diet

message("[5] Building high-quality branch event table...")

edge <- as.data.frame(tr$edge)
colnames(edge) <- c("parent_node", "child_node")

branch <- data.frame(
  branch_id = paste0("branch_", seq_len(nrow(edge))),
  parent_node = edge$parent_node,
  child_node = edge$child_node,
  branch_length = tr$edge.length,
  stringsAsFactors = FALSE
)

branch$parent_CN <- node_CN_round[as.character(branch$parent_node)]
branch$child_CN  <- node_CN_round[as.character(branch$child_node)]

branch$delta_CN <- branch$child_CN - branch$parent_CN
branch$gain_count <- ifelse(branch$delta_CN > 0, branch$delta_CN, 0)
branch$loss_count <- ifelse(branch$delta_CN < 0, abs(branch$delta_CN), 0)
branch$gain_event <- ifelse(branch$gain_count > 0, 1, 0)
branch$loss_event <- ifelse(branch$loss_count > 0, 1, 0)

branch$starting_CN <- branch$parent_CN
branch$ending_CN <- branch$child_CN

branch$parent_diet <- node_diet[as.character(branch$parent_node)]
branch$child_diet  <- node_diet[as.character(branch$child_node)]
branch$diet_state <- branch$child_diet

branch$diet_binary <- ifelse(
  branch$diet_state == "Plant_dominant",
  "Plant_dominant",
  "Non_plant"
)

branch$stable_diet_branch <- ifelse(branch$parent_diet == branch$child_diet, 1, 0)

branch$child_label <- ifelse(
  branch$child_node <= Ntip,
  tr$tip.label[branch$child_node],
  paste0("Node", branch$child_node)
)

branch$parent_label <- ifelse(
  branch$parent_node <= Ntip,
  tr$tip.label[branch$parent_node],
  paste0("Node", branch$parent_node)
)

write.table(
  branch,
  file = paste0(out_prefix, ".branch_events.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

calc_stats <- function(df, tag) {
  
  plant <- df[df$diet_binary == "Plant_dominant", ]
  non <- df[df$diet_binary == "Non_plant", ]
  
  plant_time <- sum(plant$branch_length, na.rm = TRUE)
  non_time <- sum(non$branch_length, na.rm = TRUE)
  
  plant_gain <- sum(plant$gain_count, na.rm = TRUE)
  non_gain <- sum(non$gain_count, na.rm = TRUE)
  
  plant_loss <- sum(plant$loss_count, na.rm = TRUE)
  non_loss <- sum(non$loss_count, na.rm = TRUE)
  
  plant_gain_events <- sum(plant$gain_event > 0, na.rm = TRUE)
  non_gain_events <- sum(non$gain_event > 0, na.rm = TRUE)
  
  gain_rr_corrected <- ((plant_gain + 0.5) / plant_time) /
    ((non_gain + 0.5) / non_time)
  
  loss_rr_corrected <- ((plant_loss + 0.5) / plant_time) /
    ((non_loss + 0.5) / non_time)
  
  gain_pt <- poisson.test(
    x = c(plant_gain, non_gain),
    T = c(plant_time, non_time)
  )
  
  loss_pt <- poisson.test(
    x = c(plant_loss, non_loss),
    T = c(plant_time, non_time)
  )
  
  gain_tab <- matrix(
    c(
      nrow(non) - non_gain_events, non_gain_events,
      nrow(plant) - plant_gain_events, plant_gain_events
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  rownames(gain_tab) <- c("Non_plant", "Plant_dominant")
  colnames(gain_tab) <- c("no_event", "event")
  
  gain_ft <- fisher.test(gain_tab)
  
  data.frame(
    analysis = tag,
    n_branches = nrow(df),
    plant_branches = nrow(plant),
    nonplant_branches = nrow(non),
    plant_time = plant_time,
    nonplant_time = non_time,
    plant_gain = plant_gain,
    nonplant_gain = non_gain,
    plant_gain_rate = plant_gain / plant_time,
    nonplant_gain_rate = non_gain / non_time,
    gain_rate_ratio = (plant_gain / plant_time) / (non_gain / non_time),
    gain_rate_ratio_corrected = gain_rr_corrected,
    gain_poisson_p = gain_pt$p.value,
    gain_poisson_low95 = gain_pt$conf.int[1],
    gain_poisson_high95 = gain_pt$conf.int[2],
    plant_gain_events = plant_gain_events,
    nonplant_gain_events = non_gain_events,
    gain_fisher_or = unname(gain_ft$estimate),
    gain_fisher_p = gain_ft$p.value,
    plant_loss = plant_loss,
    nonplant_loss = non_loss,
    plant_loss_rate = plant_loss / plant_time,
    nonplant_loss_rate = non_loss / non_time,
    loss_rate_ratio = (plant_loss / plant_time) / (non_loss / non_time),
    loss_rate_ratio_corrected = loss_rr_corrected,
    loss_poisson_p = loss_pt$p.value,
    loss_poisson_low95 = loss_pt$conf.int[1],
    loss_poisson_high95 = loss_pt$conf.int[2],
    stringsAsFactors = FALSE
  )
}

res_all <- calc_stats(branch, "all_high_quality_branches")
res_stable <- calc_stats(
  branch[branch$stable_diet_branch == 1, ],
  "stable_high_quality_branches"
)

res <- rbind(res_all, res_stable)

write.table(
  res,
  file = paste0(out_prefix, ".binary_rate_tests.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  dat,
  file = paste0(out_prefix, ".high_quality_species_used.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("[Done]")
message("Outputs:")
message("  ", out_prefix, ".branch_events.tsv")
message("  ", out_prefix, ".binary_rate_tests.tsv")
message("  ", out_prefix, ".high_quality_species_used.tsv")
