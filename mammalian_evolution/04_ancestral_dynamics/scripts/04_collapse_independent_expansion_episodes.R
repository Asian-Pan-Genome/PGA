#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("
Usage:
Rscript 04_collapse_independent_expansion_episodes.R tree.nwk branch_events.tsv species_order.tsv out_prefix

Example:
Rscript 04_collapse_independent_expansion_episodes.R \\
  mammal.tree.nwk \\
  301.PGA_branch_rate.branch_events.tsv \\
  species_order.tsv \\
  301.PGA.analysis2
")
}

tree_file   <- args[1]
branch_file <- args[2]
order_file  <- args[3]
out_prefix  <- args[4]

message("[1] Reading files...")

tr <- read.tree(tree_file)

branch <- read.table(
  branch_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

ord <- read.table(
  order_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

if (!all(c("species", "order") %in% colnames(ord))) {
  stop("species_order.tsv must contain columns: species, order")
}

ord$species <- gsub(" ", "_", ord$species)
order_map <- setNames(ord$order, ord$species)

missing_tips <- setdiff(tr$tip.label, names(order_map))
if (length(missing_tips) > 0) {
  warning("Some tree tips do not have order annotation. They will be marked Unknown.")
}

Ntip <- length(tr$tip.label)
Nnode <- tr$Nnode

children <- split(tr$edge[, 2], tr$edge[, 1])

desc_cache <- new.env(parent = emptyenv())

get_desc_tips <- function(node) {
  key <- as.character(node)
  if (exists(key, envir = desc_cache)) {
    return(get(key, envir = desc_cache))
  }
  if (node <= Ntip) {
    res <- node
  } else {
    ch <- children[[as.character(node)]]
    res <- unlist(lapply(ch, get_desc_tips))
  }
  assign(key, res, envir = desc_cache)
  return(res)
}

message("[2] Inferring node-level order annotation...")

all_nodes <- 1:(Ntip + Nnode)
node_order <- rep(NA_character_, length(all_nodes))
names(node_order) <- as.character(all_nodes)

node_order_set <- rep(NA_character_, length(all_nodes))
names(node_order_set) <- as.character(all_nodes)

for (node in all_nodes) {
  tips <- get_desc_tips(node)
  spp <- tr$tip.label[tips]
  orders <- unique(order_map[spp])
  orders <- orders[!is.na(orders)]
  
  if (length(orders) == 0) {
    node_order[as.character(node)] <- "Unknown"
    node_order_set[as.character(node)] <- "Unknown"
  } else if (length(orders) == 1) {
    node_order[as.character(node)] <- orders
    node_order_set[as.character(node)] <- orders
  } else {
    node_order[as.character(node)] <- "Mixed"
    node_order_set[as.character(node)] <- paste(sort(orders), collapse = ",")
  }
}

branch$parent_order <- node_order[as.character(branch$parent_node)]
branch$child_order  <- node_order[as.character(branch$child_node)]
branch$child_order_set <- node_order_set[as.character(branch$child_node)]

if (!"diet_binary" %in% colnames(branch)) {
  branch$diet_binary <- ifelse(
    branch$diet_state == "Plant_dominant",
    "Plant_dominant",
    "Non_plant"
  )
}

if (!"stable_diet_branch" %in% colnames(branch)) {
  branch$stable_diet_branch <- ifelse(branch$parent_diet == branch$child_diet, 1, 0)
}

message("[3] Finding gain branches...")

gain_branch <- branch[branch$gain_count > 0, ]

if (nrow(gain_branch) == 0) {
  stop("No gain branches found.")
}

gain_by_child <- setNames(gain_branch$branch_id, as.character(gain_branch$child_node))

parent_of_node <- setNames(branch$parent_node, as.character(branch$child_node))

find_nearest_ancestral_gain <- function(parent_node) {
  current <- parent_node
  while (!is.na(current)) {
    key <- as.character(current)
    if (key %in% names(gain_by_child)) {
      return(gain_by_child[[key]])
    }
    if (!(key %in% names(parent_of_node))) {
      return(NA_character_)
    }
    current <- parent_of_node[[key]]
  }
  return(NA_character_)
}

message("[4] Collapsing nested gain branches into independent episodes...")

gain_branch$nearest_ancestral_gain <- NA_character_
gain_branch$episode_root_branch_id <- NA_character_

for (i in seq_len(nrow(gain_branch))) {
  anc_gain <- find_nearest_ancestral_gain(gain_branch$parent_node[i])
  gain_branch$nearest_ancestral_gain[i] <- anc_gain
  
  if (is.na(anc_gain)) {
    gain_branch$episode_root_branch_id[i] <- gain_branch$branch_id[i]
  } else {
    gain_branch$episode_root_branch_id[i] <- anc_gain
  }
}

episode_list <- split(gain_branch, gain_branch$episode_root_branch_id)

episode_summary <- do.call(
  rbind,
  lapply(names(episode_list), function(ep) {
    x <- episode_list[[ep]]
    root <- branch[branch$branch_id == ep, ]
    if (nrow(root) != 1) {
      root <- x[1, ]
    }
    
    desc_tips <- get_desc_tips(root$child_node[1])
    desc_species <- tr$tip.label[desc_tips]
    
    data.frame(
      episode_id = ep,
      root_parent_node = root$parent_node[1],
      root_child_node = root$child_node[1],
      root_parent_label = root$parent_label[1],
      root_child_label = root$child_label[1],
      root_child_order = root$child_order[1],
      root_child_order_set = root$child_order_set[1],
      root_diet_state = root$diet_state[1],
      root_diet_binary = root$diet_binary[1],
      root_stable_diet_branch = root$stable_diet_branch[1],
      n_gain_branches_in_episode = nrow(x),
      total_gain_count_in_episode = sum(x$gain_count, na.rm = TRUE),
      max_single_branch_gain = max(x$gain_count, na.rm = TRUE),
      min_parent_CN = min(x$parent_CN, na.rm = TRUE),
      max_child_CN = max(x$child_CN, na.rm = TRUE),
      descendant_tip_count = length(desc_species),
      representative_descendant_tips = paste(head(desc_species, 10), collapse = ","),
      member_branch_ids = paste(x$branch_id, collapse = ","),
      member_child_labels = paste(x$child_label, collapse = ","),
      member_orders = paste(sort(unique(x$child_order)), collapse = ","),
      member_diets = paste(sort(unique(x$diet_binary)), collapse = ","),
      stringsAsFactors = FALSE
    )
  })
)

episode_summary <- episode_summary[order(
  -episode_summary$total_gain_count_in_episode,
  episode_summary$root_child_order
), ]

write.table(
  branch,
  file = paste0(out_prefix, ".branch_events.with_order.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gain_branch,
  file = paste0(out_prefix, ".gain_branches.with_episode.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  episode_summary,
  file = paste0(out_prefix, ".independent_expansion_episodes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

episode_order_summary <- aggregate(
  total_gain_count_in_episode ~ root_child_order + root_diet_binary,
  data = episode_summary,
  FUN = sum
)

episode_number_summary <- aggregate(
  episode_id ~ root_child_order + root_diet_binary,
  data = episode_summary,
  FUN = length
)

colnames(episode_number_summary)[3] <- "n_episodes"

episode_order_summary <- merge(
  episode_number_summary,
  episode_order_summary,
  by = c("root_child_order", "root_diet_binary"),
  all = TRUE
)

write.table(
  episode_order_summary,
  file = paste0(out_prefix, ".episode_order_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("[Done]")
message("Outputs:")
message("  ", out_prefix, ".branch_events.with_order.tsv")
message("  ", out_prefix, ".gain_branches.with_episode.tsv")
message("  ", out_prefix, ".independent_expansion_episodes.tsv")
message("  ", out_prefix, ".episode_order_summary.tsv")
