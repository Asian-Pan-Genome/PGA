#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(MASS)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("
Usage:
Rscript 03_test_binary_branch_rates.R branch_events.tsv out_prefix

Example:
Rscript 03_test_binary_branch_rates.R \\
  301.PGA_branch_rate.branch_events.tsv \\
  301.PGA_branch_rate.binary
")
}

branch_file <- args[1]
out_prefix  <- args[2]

message("[1] Reading branch event table...")

branch <- read.table(
  branch_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

required_cols <- c(
  "branch_length",
  "gain_count",
  "loss_count",
  "gain_event",
  "loss_event",
  "starting_CN",
  "diet_state",
  "parent_diet",
  "child_diet"
)

missing_cols <- setdiff(required_cols, colnames(branch))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

branch <- branch[!is.na(branch$branch_length) &
                   branch$branch_length > 0 &
                   !is.na(branch$diet_state) &
                   !is.na(branch$starting_CN), ]

message("[2] Branches retained: ", nrow(branch))

# 二分类：Plant_dominant vs Non_plant
branch$diet_binary <- ifelse(
  branch$diet_state == "Plant_dominant",
  "Plant_dominant",
  "Non_plant"
)

branch$diet_binary <- factor(
  branch$diet_binary,
  levels = c("Non_plant", "Plant_dominant")
)

# stable branch: parent_diet == child_diet
branch$stable_diet_branch <- ifelse(
  branch$parent_diet == branch$child_diet,
  1,
  0
)

# 输出重新编码后的 branch table
write.table(
  branch,
  file = paste0(out_prefix, ".branch_events.binary_recode.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summarize_rates <- function(df, tag, out_prefix) {
  
  message("[3] Summarizing rates for: ", tag)
  
  out <- do.call(
    rbind,
    lapply(split(df, df$diet_binary), function(x) {
      data.frame(
        diet_binary = unique(x$diet_binary),
        n_branches = nrow(x),
        total_branch_length = sum(x$branch_length, na.rm = TRUE),
        gain_branches = sum(x$gain_event > 0, na.rm = TRUE),
        gain_count = sum(x$gain_count, na.rm = TRUE),
        loss_branches = sum(x$loss_event > 0, na.rm = TRUE),
        loss_count = sum(x$loss_count, na.rm = TRUE),
        gain_rate_per_Myr = sum(x$gain_count, na.rm = TRUE) / sum(x$branch_length, na.rm = TRUE),
        loss_rate_per_Myr = sum(x$loss_count, na.rm = TRUE) / sum(x$branch_length, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(out) <- NULL
  
  # Haldane correction，避免 Non_plant 为 0 时 rate ratio 变成 Inf
  plant <- out[out$diet_binary == "Plant_dominant", ]
  non   <- out[out$diet_binary == "Non_plant", ]
  
  if (nrow(plant) == 1 && nrow(non) == 1) {
    
    gain_rr_corrected <- ((plant$gain_count + 0.5) / plant$total_branch_length) /
      ((non$gain_count + 0.5) / non$total_branch_length)
    
    loss_rr_corrected <- ((plant$loss_count + 0.5) / plant$total_branch_length) /
      ((non$loss_count + 0.5) / non$total_branch_length)
    
    rr <- data.frame(
      comparison = c("gain", "loss"),
      plant_count = c(plant$gain_count, plant$loss_count),
      nonplant_count = c(non$gain_count, non$loss_count),
      plant_time = c(plant$total_branch_length, plant$total_branch_length),
      nonplant_time = c(non$total_branch_length, non$total_branch_length),
      plant_rate = c(plant$gain_rate_per_Myr, plant$loss_rate_per_Myr),
      nonplant_rate = c(non$gain_rate_per_Myr, non$loss_rate_per_Myr),
      raw_rate_ratio = c(
        plant$gain_rate_per_Myr / non$gain_rate_per_Myr,
        plant$loss_rate_per_Myr / non$loss_rate_per_Myr
      ),
      corrected_rate_ratio = c(gain_rr_corrected, loss_rr_corrected),
      stringsAsFactors = FALSE
    )
    
    write.table(
      rr,
      file = paste0(out_prefix, ".", tag, ".rate_ratio.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  
  write.table(
    out,
    file = paste0(out_prefix, ".", tag, ".rate_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  return(out)
}

run_poisson_exact_test <- function(df, response_count, tag, out_prefix) {
  
  message("[4] Poisson exact test for ", response_count, " in ", tag)
  
  tmp <- aggregate(
    cbind(count = df[[response_count]], time = df$branch_length),
    by = list(diet_binary = df$diet_binary),
    FUN = sum
  )
  
  if (!all(c("Non_plant", "Plant_dominant") %in% tmp$diet_binary)) {
    warning("Both groups are not present for ", tag, " ", response_count)
    return(NULL)
  }
  
  plant <- tmp[tmp$diet_binary == "Plant_dominant", ]
  non   <- tmp[tmp$diet_binary == "Non_plant", ]
  
  pt <- poisson.test(
    x = c(plant$count, non$count),
    T = c(plant$time, non$time),
    alternative = "two.sided"
  )
  
  out <- data.frame(
    analysis = tag,
    response = response_count,
    plant_count = plant$count,
    nonplant_count = non$count,
    plant_time = plant$time,
    nonplant_time = non$time,
    plant_rate = plant$count / plant$time,
    nonplant_rate = non$count / non$time,
    rate_ratio = unname(pt$estimate),
    conf_low95 = pt$conf.int[1],
    conf_high95 = pt$conf.int[2],
    p_value = pt$p.value,
    method = pt$method,
    stringsAsFactors = FALSE
  )
  
  write.table(
    out,
    file = paste0(out_prefix, ".", tag, ".", response_count, ".poisson_exact.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  return(out)
}

run_fisher_event_test <- function(df, event_col, tag, out_prefix) {
  
  message("[5] Fisher event enrichment test for ", event_col, " in ", tag)
  
  event_binary <- ifelse(df[[event_col]] > 0, "event", "no_event")
  
  tab <- table(
    diet_binary = df$diet_binary,
    event = event_binary
  )
  
  # 确保 2x2 表完整
  all_rows <- c("Non_plant", "Plant_dominant")
  all_cols <- c("no_event", "event")
  
  full_tab <- matrix(
    0,
    nrow = 2,
    ncol = 2,
    dimnames = list(diet_binary = all_rows, event = all_cols)
  )
  
  full_tab[rownames(tab), colnames(tab)] <- tab
  
  ft <- fisher.test(full_tab)
  
  out <- data.frame(
    analysis = tag,
    event_col = event_col,
    nonplant_no_event = full_tab["Non_plant", "no_event"],
    nonplant_event = full_tab["Non_plant", "event"],
    plant_no_event = full_tab["Plant_dominant", "no_event"],
    plant_event = full_tab["Plant_dominant", "event"],
    odds_ratio = unname(ft$estimate),
    conf_low95 = ft$conf.int[1],
    conf_high95 = ft$conf.int[2],
    p_value = ft$p.value,
    stringsAsFactors = FALSE
  )
  
  write.table(
    out,
    file = paste0(out_prefix, ".", tag, ".", event_col, ".fisher.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  return(out)
}

fit_count_model <- function(df, response, tag, out_prefix) {
  
  message("[6] Fitting count model: ", tag, " / ", response)
  
  if (sum(df[[response]], na.rm = TRUE) == 0) {
    warning("No events for ", response, " in ", tag)
    return(NULL)
  }
  
  model_df <- df
  model_df$diet_binary <- factor(
    model_df$diet_binary,
    levels = c("Non_plant", "Plant_dominant")
  )
  
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
    paste0(response, " ~ diet_binary + starting_CN_z + offset(log_branch_length)")
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
  
  best_model <- pois_diet
  best_model_name <- "Poisson"
  
  nb_null <- NULL
  nb_diet <- NULL
  
  if (overdisp > 1.5) {
    nb_null <- try(glm.nb(formula_null, data = model_df), silent = TRUE)
    nb_diet <- try(glm.nb(formula_diet, data = model_df), silent = TRUE)
    
    if (!inherits(nb_null, "try-error") && !inherits(nb_diet, "try-error")) {
      best_model <- nb_diet
      best_model_name <- "Negative_binomial"
    }
  }
  
  sink(paste0(out_prefix, ".", tag, ".", response, ".model.txt"))
  
  cat("Analysis:", tag, "\n")
  cat("Response:", response, "\n")
  cat("Number of branches:", nrow(model_df), "\n")
  cat("Total count:", sum(model_df[[response]], na.rm = TRUE), "\n")
  cat("Poisson overdispersion:", overdisp, "\n\n")
  
  cat("=== Event distribution by diet_binary ===\n")
  print(aggregate(
    model_df[[response]],
    by = list(diet_binary = model_df$diet_binary),
    FUN = sum
  ))
  
  cat("\n=== Poisson null model ===\n")
  print(summary(pois_null))
  
  cat("\n=== Poisson diet model ===\n")
  print(summary(pois_diet))
  
  cat("\n=== Poisson model comparison: null vs diet ===\n")
  print(anova(pois_null, pois_diet, test = "Chisq"))
  
  if (!is.null(nb_null) && !is.null(nb_diet) &&
      !inherits(nb_null, "try-error") && !inherits(nb_diet, "try-error")) {
    
    cat("\n=== Negative binomial null model ===\n")
    print(summary(nb_null))
    
    cat("\n=== Negative binomial diet model ===\n")
    print(summary(nb_diet))
    
    cat("\n=== Negative binomial AIC comparison ===\n")
    print(AIC(nb_null, nb_diet))
  }
  
  cat("\n=== Best model used for coefficient table ===\n")
  cat(best_model_name, "\n")
  
  sink()
  
  coef_mat <- summary(best_model)$coefficients
  
  p_col <- grep("Pr\\(", colnames(coef_mat), value = TRUE)[1]
  stat_col <- ifelse("z value" %in% colnames(coef_mat), "z value", "t value")
  
  coef_df <- data.frame(
    analysis = tag,
    response = response,
    model = best_model_name,
    term = rownames(coef_mat),
    estimate = coef_mat[, "Estimate"],
    std_error = coef_mat[, "Std. Error"],
    statistic = coef_mat[, stat_col],
    p_value = coef_mat[, p_col],
    rate_ratio = exp(coef_mat[, "Estimate"]),
    rate_ratio_low95 = exp(coef_mat[, "Estimate"] - 1.96 * coef_mat[, "Std. Error"]),
    rate_ratio_high95 = exp(coef_mat[, "Estimate"] + 1.96 * coef_mat[, "Std. Error"]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  write.table(
    coef_df,
    file = paste0(out_prefix, ".", tag, ".", response, ".coefficients.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  return(best_model)
}

run_one_analysis <- function(df, tag, out_prefix) {
  
  message("========================================")
  message("Running analysis: ", tag)
  message("Branches: ", nrow(df))
  message("========================================")
  
  if (nrow(df) < 10) {
    warning("Too few branches for ", tag)
    return(NULL)
  }
  
  df <- df[!is.na(df$diet_binary) &
             !is.na(df$branch_length) &
             df$branch_length > 0, ]
  
  summarize_rates(df, tag, out_prefix)
  
  run_poisson_exact_test(df, "gain_count", tag, out_prefix)
  run_poisson_exact_test(df, "loss_count", tag, out_prefix)
  
  run_fisher_event_test(df, "gain_event", tag, out_prefix)
  run_fisher_event_test(df, "loss_event", tag, out_prefix)
  
  fit_count_model(df, "gain_count", tag, out_prefix)
  fit_count_model(df, "loss_count", tag, out_prefix)
}

# 1. all branches
run_one_analysis(
  df = branch,
  tag = "all_branches",
  out_prefix = out_prefix
)

# 2. stable diet branches only
stable_branch <- branch[branch$stable_diet_branch == 1, ]

run_one_analysis(
  df = stable_branch,
  tag = "stable_diet_branches",
  out_prefix = out_prefix
)

message("[Done]")
message("Main outputs:")
message("  ", out_prefix, ".all_branches.rate_summary.tsv")
message("  ", out_prefix, ".all_branches.rate_ratio.tsv")
message("  ", out_prefix, ".all_branches.gain_count.poisson_exact.tsv")
message("  ", out_prefix, ".all_branches.gain_event.fisher.tsv")
message("  ", out_prefix, ".all_branches.gain_count.coefficients.tsv")
message("  ", out_prefix, ".stable_diet_branches.rate_summary.tsv")
message("  ", out_prefix, ".stable_diet_branches.rate_ratio.tsv")
message("  ", out_prefix, ".stable_diet_branches.gain_count.poisson_exact.tsv")
message("  ", out_prefix, ".stable_diet_branches.gain_event.fisher.tsv")
message("  ", out_prefix, ".stable_diet_branches.gain_count.coefficients.tsv")
