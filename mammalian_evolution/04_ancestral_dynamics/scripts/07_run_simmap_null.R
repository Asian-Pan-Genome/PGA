#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(phytools)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("
Usage:
Rscript 07_run_simmap_null.R tree.nwk branch_events.tsv diet.tsv out_prefix n_sim [model]

Required:
1. tree.nwk
2. branch_events.tsv
3. diet.tsv, columns = species, diet_state
4. out_prefix
5. n_sim

Optional:
6. model = ER or ARD, default = ARD

Example:
Rscript 07_run_simmap_null.R \\
  mammal.tree.nwk \\
  301.PGA_branch_rate.branch_events.tsv \\
  diet.tsv \\
  301.PGA.analysis4.simmap_null \\
  1000 \\
  ARD
")
}

tree_file   <- args[1]
branch_file <- args[2]
diet_file   <- args[3]
out_prefix  <- args[4]
n_sim       <- as.integer(args[5])
model_type  <- ifelse(length(args) >= 6, args[6], "ARD")

set.seed(12345)

message("[1] Reading input files...")

tr <- read.tree(tree_file)

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

required_branch_cols <- c(
  "parent_node",
  "child_node",
  "branch_length",
  "gain_count",
  "gain_event",
  "diet_state",
  "parent_diet",
  "child_diet"
)

missing_cols <- setdiff(required_branch_cols, colnames(branch))
if (length(missing_cols) > 0) {
  stop("Missing required columns in branch_events.tsv: ",
       paste(missing_cols, collapse = ", "))
}

if (!all(c("species", "diet_state") %in% colnames(diet))) {
  stop("diet.tsv must contain columns: species, diet_state")
}

diet$species <- gsub(" ", "_", diet$species)

# branch table 中使用的是原始 tree 的 node number；
# 因此这里不要 drop tree，也不要改变 tip order，除非 tree 和 branch table 不是同一棵树。
if (!all(tr$tip.label %in% diet$species)) {
  missing_diet <- setdiff(tr$tip.label, diet$species)
  stop("Missing diet annotation for tree tips: ",
       paste(head(missing_diet, 20), collapse = ", "))
}

diet <- diet[match(tr$tip.label, diet$species), ]

diet$diet_binary <- ifelse(
  diet$diet_state == "Plant_dominant",
  "Plant_dominant",
  "Non_plant"
)

tip_state <- diet$diet_binary
names(tip_state) <- diet$species

tip_state <- factor(
  tip_state,
  levels = c("Non_plant", "Plant_dominant")
)

# 确认 branch table 和 tree edge 能够匹配
tree_edge_key <- paste(tr$edge[, 1], tr$edge[, 2], sep = "_")
branch_edge_key <- paste(branch$parent_node, branch$child_node, sep = "_")

if (!all(branch_edge_key %in% tree_edge_key)) {
  bad <- setdiff(branch_edge_key, tree_edge_key)
  stop(
    "Some branch_events edges do not match the input tree. ",
    "Make sure this is the same tree used to generate branch_events.tsv. Examples: ",
    paste(head(bad, 10), collapse = ", ")
  )
}

# 对齐 branch table 到 tree edge 顺序，方便和 simmap mapped.edge 对应
branch$edge_key <- branch_edge_key
branch_in_tree_order <- branch[match(tree_edge_key, branch$edge_key), ]

if (any(is.na(branch_in_tree_order$edge_key))) {
  stop("Failed to reorder branch table according to tree edges.")
}

message("[2] Preparing observed deterministic statistics...")

branch_in_tree_order$diet_binary <- ifelse(
  branch_in_tree_order$diet_state == "Plant_dominant",
  "Plant_dominant",
  "Non_plant"
)

branch_in_tree_order$parent_binary <- ifelse(
  branch_in_tree_order$parent_diet == "Plant_dominant",
  "Plant_dominant",
  "Non_plant"
)

branch_in_tree_order$child_binary <- ifelse(
  branch_in_tree_order$child_diet == "Plant_dominant",
  "Plant_dominant",
  "Non_plant"
)

branch_in_tree_order$stable_binary_branch <- ifelse(
  branch_in_tree_order$parent_binary == branch_in_tree_order$child_binary,
  1,
  0
)

calc_observed_deterministic <- function(df, tag) {
  
  plant <- df[df$diet_binary == "Plant_dominant", ]
  non   <- df[df$diet_binary == "Non_plant", ]
  
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
    analysis = tag,
    n_branches = nrow(df),
    plant_time = plant_time,
    nonplant_time = non_time,
    plant_gain = plant_gain,
    nonplant_gain = non_gain,
    plant_gain_fraction = plant_gain / total_gain,
    plant_gain_events = plant_event,
    nonplant_gain_events = non_event,
    plant_gain_event_fraction = plant_event / total_event,
    corrected_gain_RR = corrected_rr,
    log_corrected_gain_RR = log(corrected_rr),
    stringsAsFactors = FALSE
  )
}

obs_all <- calc_observed_deterministic(
  branch_in_tree_order,
  "observed_all_branches"
)

obs_stable <- calc_observed_deterministic(
  branch_in_tree_order[branch_in_tree_order$stable_binary_branch == 1, ],
  "observed_stable_branches"
)

message("[3] Running conditional stochastic character mapping with phytools::make.simmap...")
message("Model: ", model_type)
message("Number of simulations: ", n_sim)

# make.simmap 会固定 tip states，并在内部节点和 branches 上抽样 diet history
simmap_list <- make.simmap(
  tree = tr,
  x = tip_state,
  model = model_type,
  nsim = n_sim,
  pi = "estimated",
  message = FALSE
)

if (n_sim == 1) {
  simmap_list <- list(simmap_list)
}

extract_state_time <- function(sim_tree) {
  
  mapped <- sim_tree$mapped.edge
  
  if (is.null(mapped)) {
    stop("Simmap object does not contain mapped.edge.")
  }
  
  # 确保两个状态列都存在
  if (!"Plant_dominant" %in% colnames(mapped)) {
    mapped <- cbind(mapped, Plant_dominant = 0)
  }
  
  if (!"Non_plant" %in% colnames(mapped)) {
    mapped <- cbind(mapped, Non_plant = 0)
  }
  
  mapped <- mapped[, c("Non_plant", "Plant_dominant"), drop = FALSE]
  
  # make.simmap 后 edge 顺序通常与 sim_tree$edge 对应；
  # 这里再次用 edge_key 对齐到原始 tree edge / branch table。
  sim_edge_key <- paste(sim_tree$edge[, 1], sim_tree$edge[, 2], sep = "_")
  
  out <- data.frame(
    edge_key = sim_edge_key,
    nonplant_time_on_branch = mapped[, "Non_plant"],
    plant_time_on_branch = mapped[, "Plant_dominant"],
    stringsAsFactors = FALSE
  )
  
  out <- out[match(tree_edge_key, out$edge_key), ]
  
  if (any(is.na(out$edge_key))) {
    stop("Failed to match simmap edge order to original tree.")
  }
  
  return(out)
}

calc_simmap_stat <- function(state_time, branch_df, tag) {
  
  df <- cbind(branch_df, state_time[, c("nonplant_time_on_branch", "plant_time_on_branch")])
  
  # 每条 branch 的 gain_count 按 branch 上 Plant / Non-plant 时间比例分配
  # 这样不需要强行把整个 branch 归到某一个状态
  total_mapped_time <- df$plant_time_on_branch + df$nonplant_time_on_branch
  
  # 避免极小数值误差
  bad_time <- which(total_mapped_time <= 0 | is.na(total_mapped_time))
  if (length(bad_time) > 0) {
    df <- df[-bad_time, ]
    total_mapped_time <- total_mapped_time[-bad_time]
  }
  
  plant_fraction_on_branch <- df$plant_time_on_branch / total_mapped_time
  nonplant_fraction_on_branch <- df$nonplant_time_on_branch / total_mapped_time
  
  plant_time <- sum(df$plant_time_on_branch, na.rm = TRUE)
  nonplant_time <- sum(df$nonplant_time_on_branch, na.rm = TRUE)
  
  plant_gain <- sum(df$gain_count * plant_fraction_on_branch, na.rm = TRUE)
  nonplant_gain <- sum(df$gain_count * nonplant_fraction_on_branch, na.rm = TRUE)
  
  plant_event <- sum(df$gain_event * plant_fraction_on_branch, na.rm = TRUE)
  nonplant_event <- sum(df$gain_event * nonplant_fraction_on_branch, na.rm = TRUE)
  
  total_gain <- plant_gain + nonplant_gain
  total_event <- plant_event + nonplant_event
  
  corrected_rr <- NA_real_
  log_corrected_rr <- NA_real_
  
  if (plant_time > 0 && nonplant_time > 0) {
    corrected_rr <- ((plant_gain + 0.5) / plant_time) /
      ((nonplant_gain + 0.5) / nonplant_time)
    log_corrected_rr <- log(corrected_rr)
  }
  
  data.frame(
    analysis = tag,
    n_branches = nrow(df),
    plant_time = plant_time,
    nonplant_time = nonplant_time,
    plant_gain = plant_gain,
    nonplant_gain = nonplant_gain,
    plant_gain_fraction = ifelse(total_gain > 0, plant_gain / total_gain, NA_real_),
    plant_gain_events = plant_event,
    nonplant_gain_events = nonplant_event,
    plant_gain_event_fraction = ifelse(total_event > 0, plant_event / total_event, NA_real_),
    corrected_gain_RR = corrected_rr,
    log_corrected_gain_RR = log_corrected_rr,
    stringsAsFactors = FALSE
  )
}

message("[4] Calculating null statistics...")

null_all <- vector("list", n_sim)
null_stable <- vector("list", n_sim)

for (i in seq_len(n_sim)) {
  
  if (i %% 100 == 0) {
    message("  processed ", i, " / ", n_sim)
  }
  
  st <- extract_state_time(simmap_list[[i]])
  
  null_all[[i]] <- calc_simmap_stat(
    state_time = st,
    branch_df = branch_in_tree_order,
    tag = "null_all_branches"
  )
  
  # stable 版本：不是使用 simulated stable branch，而是沿用 observed stable branches。
  # 这样和你之前 stable analysis 的定义保持一致。
  stable_idx <- branch_in_tree_order$stable_binary_branch == 1
  
  null_stable[[i]] <- calc_simmap_stat(
    state_time = st[stable_idx, ],
    branch_df = branch_in_tree_order[stable_idx, ],
    tag = "null_stable_branches"
  )
}

null_all <- do.call(rbind, null_all)
null_stable <- do.call(rbind, null_stable)

emp_p_upper <- function(obs, null) {
  null <- null[is.finite(null) & !is.na(null)]
  if (length(null) == 0) return(NA_real_)
  (sum(null >= obs, na.rm = TRUE) + 1) / (length(null) + 1)
}

summarize_null <- function(obs, null, tag) {
  
  data.frame(
    analysis = tag,
    
    observed_plant_gain = obs$plant_gain,
    observed_nonplant_gain = obs$nonplant_gain,
    observed_plant_gain_fraction = obs$plant_gain_fraction,
    null_median_plant_gain_fraction = median(null$plant_gain_fraction, na.rm = TRUE),
    null_low95_plant_gain_fraction = quantile(null$plant_gain_fraction, 0.025, na.rm = TRUE),
    null_high95_plant_gain_fraction = quantile(null$plant_gain_fraction, 0.975, na.rm = TRUE),
    empirical_p_plant_gain_fraction = emp_p_upper(
      obs$plant_gain_fraction,
      null$plant_gain_fraction
    ),
    
    observed_plant_gain_event_fraction = obs$plant_gain_event_fraction,
    null_median_plant_gain_event_fraction = median(null$plant_gain_event_fraction, na.rm = TRUE),
    null_low95_plant_gain_event_fraction = quantile(null$plant_gain_event_fraction, 0.025, na.rm = TRUE),
    null_high95_plant_gain_event_fraction = quantile(null$plant_gain_event_fraction, 0.975, na.rm = TRUE),
    empirical_p_plant_gain_event_fraction = emp_p_upper(
      obs$plant_gain_event_fraction,
      null$plant_gain_event_fraction
    ),
    
    observed_corrected_gain_RR = obs$corrected_gain_RR,
    observed_log_corrected_gain_RR = obs$log_corrected_gain_RR,
    null_median_log_RR = median(null$log_corrected_gain_RR, na.rm = TRUE),
    null_low95_log_RR = quantile(null$log_corrected_gain_RR, 0.025, na.rm = TRUE),
    null_high95_log_RR = quantile(null$log_corrected_gain_RR, 0.975, na.rm = TRUE),
    empirical_p_log_RR = emp_p_upper(
      obs$log_corrected_gain_RR,
      null$log_corrected_gain_RR
    ),
    
    null_valid_log_RR_n = sum(is.finite(null$log_corrected_gain_RR)),
    null_total_n = nrow(null),
    
    observed_plant_time = obs$plant_time,
    observed_nonplant_time = obs$nonplant_time,
    null_median_plant_time = median(null$plant_time, na.rm = TRUE),
    null_median_nonplant_time = median(null$nonplant_time, na.rm = TRUE),
    
    stringsAsFactors = FALSE
  )
}

summary_all <- summarize_null(
  obs = obs_all,
  null = null_all,
  tag = "all_branches"
)

summary_stable <- summarize_null(
  obs = obs_stable,
  null = null_stable,
  tag = "stable_branches"
)

out_summary <- rbind(summary_all, summary_stable)

message("[5] Writing outputs...")

write.table(
  out_summary,
  file = paste0(out_prefix, ".phylogenetic_null.summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  obs_all,
  file = paste0(out_prefix, ".observed_all_branches.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  obs_stable,
  file = paste0(out_prefix, ".observed_stable_branches.tsv"),
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
message("Main output:")
message("  ", out_prefix, ".phylogenetic_null.summary.tsv")
message("")
message("Check these columns:")
message("  empirical_p_plant_gain_fraction")
message("  empirical_p_plant_gain_event_fraction")
message("  empirical_p_log_RR")
message("  null_valid_log_RR_n")
message("  null_median_plant_time")
message("  null_median_nonplant_time")
