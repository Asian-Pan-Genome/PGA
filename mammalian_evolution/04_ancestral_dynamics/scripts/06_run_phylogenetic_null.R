#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("
Usage:
Rscript 06_run_phylogenetic_null.R tree.nwk branch_events.tsv diet.tsv out_prefix n_sim

Example:
Rscript 06_run_phylogenetic_null.R \\
  mammal.tree.nwk \\
  301.PGA_branch_rate.branch_events.tsv \\
  diet.tsv \\
  301.PGA.analysis4.phylo_null \\
  10000
")
}

tree_file   <- args[1]
branch_file <- args[2]
diet_file   <- args[3]
out_prefix  <- args[4]
n_sim       <- as.integer(args[5])

set.seed(12345)

message("[1] Reading files...")

tr <- read.tree(tree_file)
tr <- reorder.phylo(tr, "cladewise")

branch <- read.table(
  branch_file,
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

if (!all(c("species", "diet_state") %in% colnames(diet))) {
  stop("diet.tsv must contain columns: species, diet_state")
}

diet$species <- gsub(" ", "_", diet$species)
diet$diet_binary <- ifelse(
  diet$diet_state == "Plant_dominant",
  "Plant_dominant",
  "Non_plant"
)

missing <- setdiff(tr$tip.label, diet$species)
if (length(missing) > 0) {
  stop("Missing diet annotation for tree tips: ", paste(head(missing, 20), collapse = ", "))
}

diet <- diet[match(tr$tip.label, diet$species), ]

tip_state <- factor(
  diet$diet_binary,
  levels = c("Non_plant", "Plant_dominant")
)
names(tip_state) <- diet$species

message("[2] Fitting binary diet transition model using ace...")

fit <- ace(
  x = tip_state,
  phy = tr,
  type = "discrete",
  model = "ARD"
)

states <- c("Non_plant", "Plant_dominant")
Q <- matrix(0, nrow = 2, ncol = 2, dimnames = list(states, states))

idx <- fit$index.matrix
rates <- fit$rates

for (i in seq_len(nrow(idx))) {
  for (j in seq_len(ncol(idx))) {
    if (i != j && !is.na(idx[i, j])) {
      Q[i, j] <- rates[idx[i, j]]
    }
  }
}

Q[Q < 1e-12] <- 1e-12
diag(Q) <- -rowSums(Q)

q01 <- Q["Non_plant", "Plant_dominant"]
q10 <- Q["Plant_dominant", "Non_plant"]

message("Estimated transition rates:")
message("  Non_plant -> Plant_dominant: ", signif(q01, 4))
message("  Plant_dominant -> Non_plant: ", signif(q10, 4))

Ntip <- length(tr$tip.label)
Nnode <- tr$Nnode
root_node <- Ntip + 1

children <- split(tr$edge[, 2], tr$edge[, 1])
edge_length <- setNames(tr$edge.length, paste(tr$edge[,1], tr$edge[,2], sep = "_"))

simulate_states <- function() {
  
  node_state <- rep(NA_character_, Ntip + Nnode)
  names(node_state) <- as.character(1:(Ntip + Nnode))
  
  s <- q01 + q10
  
  # stationary root prior
  pi_non <- q10 / s
  pi_plant <- q01 / s
  
  node_state[as.character(root_node)] <- sample(
    states,
    size = 1,
    prob = c(pi_non, pi_plant)
  )
  
  for (e in seq_len(nrow(tr$edge))) {
    parent <- tr$edge[e, 1]
    child  <- tr$edge[e, 2]
    t <- tr$edge.length[e]
    
    parent_state <- node_state[as.character(parent)]
    
    if (parent_state == "Non_plant") {
      p_plant <- q01 / s * (1 - exp(-s * t))
      child_state <- sample(
        states,
        size = 1,
        prob = c(1 - p_plant, p_plant)
      )
    } else {
      p_non <- q10 / s * (1 - exp(-s * t))
      child_state <- sample(
        states,
        size = 1,
        prob = c(p_non, 1 - p_non)
      )
    }
    
    node_state[as.character(child)] <- child_state
  }
  
  return(node_state)
}

calc_stat <- function(branch_df, node_state = NULL, use_observed = FALSE, stable_only = FALSE) {
  
  df <- branch_df
  
  if (use_observed) {
    df$diet_binary_tmp <- ifelse(
      df$diet_state == "Plant_dominant",
      "Plant_dominant",
      "Non_plant"
    )
    df$parent_binary_tmp <- ifelse(
      df$parent_diet == "Plant_dominant",
      "Plant_dominant",
      "Non_plant"
    )
    df$child_binary_tmp <- ifelse(
      df$child_diet == "Plant_dominant",
      "Plant_dominant",
      "Non_plant"
    )
  } else {
    df$parent_binary_tmp <- node_state[as.character(df$parent_node)]
    df$child_binary_tmp  <- node_state[as.character(df$child_node)]
    df$diet_binary_tmp   <- df$child_binary_tmp
  }
  
  if (stable_only) {
    df <- df[df$parent_binary_tmp == df$child_binary_tmp, ]
  }
  
  plant <- df[df$diet_binary_tmp == "Plant_dominant", ]
  non   <- df[df$diet_binary_tmp == "Non_plant", ]
  
  plant_time <- sum(plant$branch_length, na.rm = TRUE)
  non_time   <- sum(non$branch_length, na.rm = TRUE)
  
  plant_gain <- sum(plant$gain_count, na.rm = TRUE)
  non_gain   <- sum(non$gain_count, na.rm = TRUE)
  
  plant_event <- sum(plant$gain_event > 0, na.rm = TRUE)
  non_event   <- sum(non$gain_event > 0, na.rm = TRUE)
  
  total_gain <- plant_gain + non_gain
  total_event <- plant_event + non_event
  
  corrected_rr <- ((plant_gain + 0.5) / plant_time) /
    ((non_gain + 0.5) / non_time)
  
  data.frame(
    n_branches = nrow(df),
    plant_time = plant_time,
    nonplant_time = non_time,
    plant_gain = plant_gain,
    nonplant_gain = non_gain,
    plant_gain_fraction = ifelse(total_gain > 0, plant_gain / total_gain, NA_real_),
    corrected_gain_rate_ratio = corrected_rr,
    log_corrected_gain_rate_ratio = log(corrected_rr),
    plant_gain_events = plant_event,
    nonplant_gain_events = non_event,
    plant_gain_event_fraction = ifelse(total_event > 0, plant_event / total_event, NA_real_)
  )
}

message("[3] Calculating observed statistics...")

obs_all <- calc_stat(branch, use_observed = TRUE, stable_only = FALSE)
obs_stable <- calc_stat(branch, use_observed = TRUE, stable_only = TRUE)

message("[4] Running phylogenetic null simulations: ", n_sim)

null_all <- vector("list", n_sim)
null_stable <- vector("list", n_sim)

for (i in seq_len(n_sim)) {
  if (i %% 1000 == 0) {
    message("  simulation ", i, " / ", n_sim)
  }
  
  sim_state <- simulate_states()
  
  null_all[[i]] <- calc_stat(
    branch,
    node_state = sim_state,
    use_observed = FALSE,
    stable_only = FALSE
  )
  
  null_stable[[i]] <- calc_stat(
    branch,
    node_state = sim_state,
    use_observed = FALSE,
    stable_only = TRUE
  )
}

null_all <- do.call(rbind, null_all)
null_stable <- do.call(rbind, null_stable)

emp_p <- function(obs, null) {
  (sum(null >= obs, na.rm = TRUE) + 1) / (sum(!is.na(null)) + 1)
}

summary_one <- function(obs, null, tag) {
  data.frame(
    analysis = tag,
    observed_plant_gain = obs$plant_gain,
    observed_nonplant_gain = obs$nonplant_gain,
    observed_plant_gain_fraction = obs$plant_gain_fraction,
    observed_corrected_gain_RR = obs$corrected_gain_rate_ratio,
    observed_log_corrected_gain_RR = obs$log_corrected_gain_rate_ratio,
    null_median_plant_gain_fraction = median(null$plant_gain_fraction, na.rm = TRUE),
    null_low95_plant_gain_fraction = quantile(null$plant_gain_fraction, 0.025, na.rm = TRUE),
    null_high95_plant_gain_fraction = quantile(null$plant_gain_fraction, 0.975, na.rm = TRUE),
    empirical_p_plant_gain_fraction = emp_p(
      obs$plant_gain_fraction,
      null$plant_gain_fraction
    ),
    null_median_log_RR = median(null$log_corrected_gain_rate_ratio, na.rm = TRUE),
    null_low95_log_RR = quantile(null$log_corrected_gain_rate_ratio, 0.025, na.rm = TRUE),
    null_high95_log_RR = quantile(null$log_corrected_gain_rate_ratio, 0.975, na.rm = TRUE),
    empirical_p_log_RR = emp_p(
      obs$log_corrected_gain_rate_ratio,
      null$log_corrected_gain_rate_ratio
    ),
    observed_plant_gain_event_fraction = obs$plant_gain_event_fraction,
    null_median_plant_gain_event_fraction = median(null$plant_gain_event_fraction, na.rm = TRUE),
    empirical_p_plant_gain_event_fraction = emp_p(
      obs$plant_gain_event_fraction,
      null$plant_gain_event_fraction
    ),
    stringsAsFactors = FALSE
  )
}

out_summary <- rbind(
  summary_one(obs_all, null_all, "all_branches"),
  summary_one(obs_stable, null_stable, "stable_branches")
)

write.table(
  out_summary,
  file = paste0(out_prefix, ".phylogenetic_null.summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  null_all,
  file = paste0(out_prefix, ".null_all_branches.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  null_stable,
  file = paste0(out_prefix, ".null_stable_branches.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("[Done]")
message("Outputs:")
message("  ", out_prefix, ".phylogenetic_null.summary.tsv")
message("  ", out_prefix, ".null_all_branches.tsv")
message("  ", out_prefix, ".null_stable_branches.tsv")
