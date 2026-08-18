#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(MASS)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("
Usage:
Rscript 02_reconstruct_branch_cn_and_test_rates.R tree.nwk PGA_CN.tsv diet.tsv out_prefix

Input:
1. tree.nwk: Newick tree with branch lengths
2. PGA_CN.tsv: columns = species, PGA_CN
3. diet.tsv: columns = species, diet_state
4. out_prefix: output prefix

Example:
Rscript 02_reconstruct_branch_cn_and_test_rates.R mammal.tree.nwk PGA_CN.tsv diet.tsv PGA_branch_rate
")
}

tree_file <- args[1]
cn_file   <- args[2]
diet_file <- args[3]
out_pref  <- args[4]

message("[1] Reading input files...")

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

if (!all(c("species", "PGA_CN") %in% colnames(cn))) {
  stop("PGA_CN.tsv must contain columns: species, PGA_CN")
}

if (!all(c("species", "diet_state") %in% colnames(diet))) {
  stop("diet.tsv must contain columns: species, diet_state")
}

cn$species <- gsub(" ", "_", cn$species)
diet$species <- gsub(" ", "_", diet$species)

dat <- merge(cn, diet, by = "species")
dat <- dat[!is.na(dat$PGA_CN) & !is.na(dat$diet_state), ]

common_tips <- intersect(tr$tip.label, dat$species)

if (length(common_tips) < 10) {
  stop("Too few matched species between tree, CN, and diet files.")
}

message("[2] Matched species: ", length(common_tips))

tr <- drop.tip(tr, setdiff(tr$tip.label, common_tips))
dat <- dat[match(tr$tip.label, dat$species), ]

if (any(dat$species != tr$tip.label)) {
  stop("Species order mismatch after matching.")
}

if (is.null(tr$edge.length)) {
  stop("Tree must have branch lengths.")
}

if (any(is.na(tr$edge.length))) {
  stop("Tree contains NA branch lengths.")
}

min_positive_bl <- min(tr$edge.length[tr$edge.length > 0], na.rm = TRUE)

if (!is.finite(min_positive_bl)) {
  stop("All branch lengths are zero or invalid.")
}

tr$edge.length[tr$edge.length <= 0] <- min_positive_bl / 10

Ntip <- length(tr$tip.label)
Nnode <- tr$Nnode
all_nodes <- 1:(Ntip + Nnode)
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
names(node_CN) <- all_nodes

node_CN[1:Ntip] <- cn_vec[tr$tip.label]

if (is.null(names(cn_ace$ace))) {
  names(cn_ace$ace) <- internal_nodes
}

node_CN[as.integer(names(cn_ace$ace))] <- cn_ace$ace

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
names(node_diet) <- all_nodes

node_diet[1:Ntip] <- as.character(diet_vec[tr$tip.label])

lik <- diet_ace$lik.anc

if (is.null(rownames(lik))) {
  rownames(lik) <- internal_nodes
}

internal_diet <- apply(lik, 1, function(x) {
  colnames(lik)[which.max(x)]
})

node_diet[as.integer(rownames(lik))] <- internal_diet

message("[5] Building branch-level event table...")

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

# 最简策略：把 branch 的 diet 定义为 child node 的 diet
branch$diet_state <- branch$child_diet

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

branch$diet_transition <- ifelse(branch$parent_diet == branch$child_diet, 0, 1)

branch <- branch[!is.na(branch$parent_CN) &
                   !is.na(branch$child_CN) &
                   !is.na(branch$diet_state) &
                   !is.na(branch$branch_length), ]

write.table(
  branch,
  file = paste0(out_pref, ".branch_events.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("[6] Branch event table written: ", out_pref, ".branch_events.tsv")

fit_count_model <- function(df, response, out_pref) {
  
  message("[7] Fitting model for: ", response)
  
  if (sum(df[[response]], na.rm = TRUE) == 0) {
    warning("No events for ", response, ". Skip model.")
    return(NULL)
  }
  
  model_df <- df
  model_df$diet_state <- as.factor(model_df$diet_state)
  model_df$log_branch_length <- log(model_df$branch_length)
  
  if (sd(model_df$starting_CN, na.rm = TRUE) == 0) {
    model_df$starting_CN_z <- 0
  } else {
    model_df$starting_CN_z <- as.numeric(scale(model_df$starting_CN))
  }
  
  formula_null <- as.formula(
    paste0(response, " ~ starting_CN_z + offset(log_branch_length)")
  )
  
  formula_diet <- as.formula(
    paste0(response, " ~ diet_state + starting_CN_z + offset(log_branch_length)")
  )
  
  pois_null <- glm(
    formula_null,
    family = poisson(),
    data = model_df
  )
  
  pois_diet <- glm(
    formula_diet,
    family = poisson(),
    data = model_df
  )
  
  overdisp <- sum(residuals(pois_diet, type = "pearson")^2) / df.residual(pois_diet)
  
  sink(paste0(out_pref, ".", response, ".model.txt"))
  
  cat("Response:", response, "\n")
  cat("Total events:", sum(model_df[[response]], na.rm = TRUE), "\n")
  cat("Number of branches:", nrow(model_df), "\n")
  cat("Poisson overdispersion:", overdisp, "\n\n")
  
  cat("=== Poisson null model ===\n")
  print(summary(pois_null))
  
  cat("\n=== Poisson diet model ===\n")
  print(summary(pois_diet))
  
  cat("\n=== Poisson model comparison: null vs diet ===\n")
  print(anova(pois_null, pois_diet, test = "Chisq"))
  
  best_model <- pois_diet
  best_model_name <- "Poisson"
  
  if (overdisp > 1.5) {
    cat("\nOverdispersion > 1.5, trying negative binomial model...\n")
    
    nb_null <- try(glm.nb(formula_null, data = model_df), silent = TRUE)
    nb_diet <- try(glm.nb(formula_diet, data = model_df), silent = TRUE)
    
    if (!inherits(nb_null, "try-error") && !inherits(nb_diet, "try-error")) {
      cat("\n=== Negative binomial null model ===\n")
      print(summary(nb_null))
      
      cat("\n=== Negative binomial diet model ===\n")
      print(summary(nb_diet))
      
      cat("\n=== Negative binomial AIC comparison ===\n")
      print(AIC(nb_null, nb_diet))
      
      best_model <- nb_diet
      best_model_name <- "Negative_binomial"
    } else {
      cat("\nNegative binomial model failed. Keep Poisson result.\n")
    }
  }
  
  cat("\n=== Best model used for coefficient table ===\n")
  cat(best_model_name, "\n")
  
  sink()
  
  coef_mat <- summary(best_model)$coefficients
  coef_df <- data.frame(
    term = rownames(coef_mat),
    estimate = coef_mat[, "Estimate"],
    std_error = coef_mat[, "Std. Error"],
    statistic = coef_mat[, ifelse("z value" %in% colnames(coef_mat), "z value", "t value")],
    p_value = coef_mat[, grep("Pr\\(", colnames(coef_mat), value = TRUE)[1]],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  coef_df$rate_ratio <- exp(coef_df$estimate)
  coef_df$rate_ratio_low95 <- exp(coef_df$estimate - 1.96 * coef_df$std_error)
  coef_df$rate_ratio_high95 <- exp(coef_df$estimate + 1.96 * coef_df$std_error)
  coef_df$model <- best_model_name
  coef_df$response <- response
  
  write.table(
    coef_df,
    file = paste0(out_pref, ".", response, ".coefficients.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  return(best_model)
}

gain_model <- fit_count_model(branch, "gain_count", out_pref)
loss_model <- fit_count_model(branch, "loss_count", out_pref)

message("[8] Done.")
message("Main outputs:")
message("  ", out_pref, ".branch_events.tsv")
message("  ", out_pref, ".gain_count.model.txt")
message("  ", out_pref, ".gain_count.coefficients.tsv")
message("  ", out_pref, ".loss_count.model.txt")
message("  ", out_pref, ".loss_count.coefficients.tsv")
