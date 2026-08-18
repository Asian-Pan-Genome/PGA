#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("
Usage:
Rscript 05_leave_one_order_out.R branch_events.with_order.tsv out_prefix

Example:
Rscript 05_leave_one_order_out.R \\
  301.PGA.analysis2.branch_events.with_order.tsv \\
  301.PGA.analysis3.LOO
")
}

branch_file <- args[1]
out_prefix  <- args[2]

branch <- read.table(
  branch_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

required <- c(
  "branch_length", "gain_count", "loss_count",
  "gain_event", "loss_event",
  "diet_state", "parent_diet", "child_diet",
  "diet_binary", "child_order"
)

missing <- setdiff(required, colnames(branch))
if (length(missing) > 0) {
  stop("Missing columns: ", paste(missing, collapse = ", "))
}

branch <- branch[!is.na(branch$branch_length) & branch$branch_length > 0, ]

branch$diet_binary <- ifelse(
  branch$diet_state == "Plant_dominant",
  "Plant_dominant",
  "Non_plant"
)

branch$stable_diet_branch <- ifelse(branch$parent_diet == branch$child_diet, 1, 0)

calc_stats <- function(df, tag) {
  
  df$diet_binary <- ifelse(df$diet_state == "Plant_dominant", "Plant_dominant", "Non_plant")
  
  plant <- df[df$diet_binary == "Plant_dominant", ]
  non   <- df[df$diet_binary == "Non_plant", ]
  
  plant_time <- sum(plant$branch_length, na.rm = TRUE)
  non_time   <- sum(non$branch_length, na.rm = TRUE)
  
  plant_gain <- sum(plant$gain_count, na.rm = TRUE)
  non_gain   <- sum(non$gain_count, na.rm = TRUE)
  
  plant_loss <- sum(plant$loss_count, na.rm = TRUE)
  non_loss   <- sum(non$loss_count, na.rm = TRUE)
  
  plant_gain_events <- sum(plant$gain_event > 0, na.rm = TRUE)
  non_gain_events   <- sum(non$gain_event > 0, na.rm = TRUE)
  
  plant_loss_events <- sum(plant$loss_event > 0, na.rm = TRUE)
  non_loss_events   <- sum(non$loss_event > 0, na.rm = TRUE)
  
  gain_rate_plant <- plant_gain / plant_time
  gain_rate_non   <- non_gain / non_time
  
  loss_rate_plant <- plant_loss / plant_time
  loss_rate_non   <- non_loss / non_time
  
  gain_rr <- gain_rate_plant / gain_rate_non
  loss_rr <- loss_rate_plant / loss_rate_non
  
  gain_rr_corrected <- ((plant_gain + 0.5) / plant_time) /
    ((non_gain + 0.5) / non_time)
  
  loss_rr_corrected <- ((plant_loss + 0.5) / plant_time) /
    ((non_loss + 0.5) / non_time)
  
  gain_pois_p <- NA_real_
  loss_pois_p <- NA_real_
  gain_pois_low <- NA_real_
  gain_pois_high <- NA_real_
  loss_pois_low <- NA_real_
  loss_pois_high <- NA_real_
  
  if (plant_time > 0 && non_time > 0) {
    gain_pt <- poisson.test(
      x = c(plant_gain, non_gain),
      T = c(plant_time, non_time)
    )
    loss_pt <- poisson.test(
      x = c(plant_loss, non_loss),
      T = c(plant_time, non_time)
    )
    
    gain_pois_p <- gain_pt$p.value
    gain_pois_low <- gain_pt$conf.int[1]
    gain_pois_high <- gain_pt$conf.int[2]
    
    loss_pois_p <- loss_pt$p.value
    loss_pois_low <- loss_pt$conf.int[1]
    loss_pois_high <- loss_pt$conf.int[2]
  }
  
  gain_fisher_p <- NA_real_
  gain_fisher_or <- NA_real_
  loss_fisher_p <- NA_real_
  loss_fisher_or <- NA_real_
  
  if (nrow(plant) > 0 && nrow(non) > 0) {
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
    
    loss_tab <- matrix(
      c(
        nrow(non) - non_loss_events, non_loss_events,
        nrow(plant) - plant_loss_events, plant_loss_events
      ),
      nrow = 2,
      byrow = TRUE
    )
    rownames(loss_tab) <- c("Non_plant", "Plant_dominant")
    colnames(loss_tab) <- c("no_event", "event")
    
    gf <- fisher.test(gain_tab)
    lf <- fisher.test(loss_tab)
    
    gain_fisher_p <- gf$p.value
    gain_fisher_or <- unname(gf$estimate)
    
    loss_fisher_p <- lf$p.value
    loss_fisher_or <- unname(lf$estimate)
  }
  
  data.frame(
    analysis = tag,
    n_branches = nrow(df),
    plant_branches = nrow(plant),
    nonplant_branches = nrow(non),
    plant_time = plant_time,
    nonplant_time = non_time,
    plant_gain = plant_gain,
    nonplant_gain = non_gain,
    plant_loss = plant_loss,
    nonplant_loss = non_loss,
    plant_gain_events = plant_gain_events,
    nonplant_gain_events = non_gain_events,
    plant_loss_events = plant_loss_events,
    nonplant_loss_events = non_loss_events,
    plant_gain_rate = gain_rate_plant,
    nonplant_gain_rate = gain_rate_non,
    gain_rate_ratio = gain_rr,
    gain_rate_ratio_corrected = gain_rr_corrected,
    gain_poisson_p = gain_pois_p,
    gain_poisson_low95 = gain_pois_low,
    gain_poisson_high95 = gain_pois_high,
    gain_fisher_or = gain_fisher_or,
    gain_fisher_p = gain_fisher_p,
    plant_loss_rate = loss_rate_plant,
    nonplant_loss_rate = loss_rate_non,
    loss_rate_ratio = loss_rr,
    loss_rate_ratio_corrected = loss_rr_corrected,
    loss_poisson_p = loss_pois_p,
    loss_poisson_low95 = loss_pois_low,
    loss_poisson_high95 = loss_pois_high,
    loss_fisher_or = loss_fisher_or,
    loss_fisher_p = loss_fisher_p,
    stringsAsFactors = FALSE
  )
}

orders <- sort(unique(branch$child_order))
orders <- setdiff(orders, c("Mixed", "Unknown", NA))

res <- list()

res[["all_orders_all_branches"]] <- calc_stats(branch, "all_orders_all_branches")

stable <- branch[branch$stable_diet_branch == 1, ]
res[["all_orders_stable_branches"]] <- calc_stats(stable, "all_orders_stable_branches")

for (ord in orders) {
  message("[LOO] Removing order: ", ord)
  
  df1 <- branch[branch$child_order != ord, ]
  df2 <- stable[stable$child_order != ord, ]
  
  res[[paste0("remove_", ord, "_all_branches")]] <- calc_stats(
    df1,
    paste0("remove_", ord, "_all_branches")
  )
  
  res[[paste0("remove_", ord, "_stable_branches")]] <- calc_stats(
    df2,
    paste0("remove_", ord, "_stable_branches")
  )
}

out <- do.call(rbind, res)

write.table(
  out,
  file = paste0(out_prefix, ".leave_one_order_out.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("[Done]")
message("Output:")
message("  ", out_prefix, ".leave_one_order_out.tsv")
