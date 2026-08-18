#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(nlme)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("
Usage:
Rscript 04_run_pgls_diet_models.R <input.tsv> <tree.nwk> <out_prefix>

Example:
Rscript 04_run_pgls_diet_models.R \\
  301.PGA_diet_bodymass.high_order_ratio.tsv \\
  mammal_time_tree.nwk \\
  PGA_high_order_diet_PGLS
")
}

infile <- args[1]
treefile <- args[2]
out_prefix <- args[3]

dat <- read.table(
  infile,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE
)

tree <- read.tree(treefile)

# -----------------------------------------
# basic cleanup
# -----------------------------------------

dat$tip_label <- gsub(" ", "_", dat$tip_label)

dat$CNV_log1p <- log1p(dat$CNV)
dat$log10_body_mass <- log10(dat$`BodyMass-Value`)

dat$plant_dominant_prop <- as.numeric(dat$plant_dominant_prop)
dat$carnivore_prop <- as.numeric(dat$carnivore_prop)
dat$insectivore_prop <- as.numeric(dat$insectivore_prop)
dat$nectarivore_prop <- as.numeric(dat$nectarivore_prop)

dat$High_Order_Diet <- as.factor(dat$High_Order_Diet)

# Reference group
# Carnivore is a useful reference because Plant_dominant coefficient
# then means Plant_dominant vs Carnivore after controlling body mass.
if ("Carnivore" %in% levels(dat$High_Order_Diet)) {
  dat$High_Order_Diet <- relevel(dat$High_Order_Diet, ref = "Carnivore")
}

dat$Use_High_Order_Diet <- as.character(dat$Use_High_Order_Diet)
dat$Use_High_Order_Diet <- dat$Use_High_Order_Diet %in% c("TRUE", "True", "true", "1")
# -----------------------------------------
# match tree and data
# -----------------------------------------

common <- intersect(tree$tip.label, dat$tip_label)

cat("Species in data:", length(unique(dat$tip_label)), "\n")
cat("Tips in tree:", length(tree$tip.label), "\n")
cat("Matched species:", length(common), "\n")

tree <- drop.tip(tree, setdiff(tree$tip.label, common))
dat <- dat[dat$tip_label %in% common, ]

# avoid zero branch length problems
if (!is.null(tree$edge.length)) {
  positive_edges <- tree$edge.length[tree$edge.length > 0]
  if (length(positive_edges) > 0) {
    min_positive <- min(positive_edges)
    tree$edge.length[tree$edge.length <= 0] <- min_positive / 2
  }
}

# -----------------------------------------
# PGLS function
# -----------------------------------------

fit_pgls <- function(formula, data, tree, model_name) {

  vars <- all.vars(formula)
  needed <- unique(c("tip_label", vars))

  d <- data[complete.cases(data[, needed]), ]

  tr <- drop.tip(tree, setdiff(tree$tip.label, d$tip_label))
  d <- d[match(tr$tip.label, d$tip_label), ]

  d$tip_label <- factor(d$tip_label, levels = tr$tip.label)

  fit <- gls(
    formula,
    data = d,
    correlation = corPagel(
      value = 0.5,
      phy = tr,
      fixed = FALSE,
      form = ~ tip_label
    ),
    method = "ML",
    na.action = na.omit
  )

  sm <- summary(fit)

  coef_tab <- as.data.frame(sm$tTable)
  coef_tab$term <- rownames(coef_tab)
  rownames(coef_tab) <- NULL

  coef_tab$model <- model_name
  coef_tab$n_species <- nrow(d)
  coef_tab$AIC <- AIC(fit)

  lambda <- tryCatch(
    coef(fit$modelStruct$corStruct, unconstrained = FALSE),
    error = function(e) NA
  )

  coef_tab$lambda <- as.numeric(lambda)
  coef_tab$CI_low <- coef_tab$Value - 1.96 * coef_tab$Std.Error
  coef_tab$CI_high <- coef_tab$Value + 1.96 * coef_tab$Std.Error

  list(
    fit = fit,
    coef = coef_tab,
    data = d,
    tree = tr
  )
}

# -----------------------------------------
# Model 1:
# CNV ~ high-order diet + body mass
# -----------------------------------------

dat_class <- dat[
  dat$Use_High_Order_Diet == TRUE &
    dat$High_Order_Diet %in% c(
      "Plant_dominant",
      "Omnivore",
      "Carnivore",
      "Insectivore"
    ),
]

dat_class$High_Order_Diet <- droplevels(dat_class$High_Order_Diet)

if ("Carnivore" %in% levels(dat_class$High_Order_Diet)) {
  dat_class$High_Order_Diet <- relevel(dat_class$High_Order_Diet, ref = "Carnivore")
}

cat("Rows in dat_class:", nrow(dat_class), "\n")
print(table(dat_class$High_Order_Diet, useNA = "ifany"))
cat("Matched species:", length(intersect(tree$tip.label, dat_class$tip_label)), "\n")

m_class <- fit_pgls(
  CNV_log1p ~ High_Order_Diet + log10_body_mass,
  dat_class,
  tree,
  "CNV_log1p ~ High_Order_Diet + log10_body_mass"
)

# -----------------------------------------
# Model 2:
# CNV ~ corresponding food proportions + body mass
# -----------------------------------------
#
# Important:
# Diet proportions are compositional.
# Do NOT include plant + carnivore + insectivore + intercept together
# if they sum to 1, because this can cause collinearity.
#
# Here we use carnivore_prop as the implicit reference component.
# plant_dominant_prop tests plant food proportion relative to carnivory.
# insectivore_prop controls insect-based diet proportion.
#
# If nectarivore_prop is common in your data, you can include it too.
# By default, we include it only if any species has nectarivore_prop > 0.
# -----------------------------------------

if (any(dat$nectarivore_prop > 0, na.rm = TRUE)) {

  m_prop <- fit_pgls(
    CNV_log1p ~ plant_dominant_prop + insectivore_prop + nectarivore_prop + log10_body_mass,
    dat,
    tree,
    "CNV_log1p ~ plant_dominant_prop + insectivore_prop + nectarivore_prop + log10_body_mass"
  )

} else {

  m_prop <- fit_pgls(
    CNV_log1p ~ plant_dominant_prop + insectivore_prop + log10_body_mass,
    dat,
    tree,
    "CNV_log1p ~ plant_dominant_prop + insectivore_prop + log10_body_mass"
  )

}

# Also run a simple plant-ratio model, useful for main-text interpretation
m_prop_simple <- fit_pgls(
  CNV_log1p ~ plant_dominant_prop + log10_body_mass,
  dat,
  tree,
  "CNV_log1p ~ plant_dominant_prop + log10_body_mass"
)

# -----------------------------------------
# output
# -----------------------------------------

sink(paste0(out_prefix, ".model_summary.txt"))

cat("\n==============================\n")
cat("Model 1: High-order diet + body mass\n")
cat("==============================\n")
print(summary(m_class$fit))

cat("\n==============================\n")
cat("Model 2: Food proportions + body mass\n")
cat("==============================\n")
print(summary(m_prop$fit))

cat("\n==============================\n")
cat("Model 2 simple: Plant-dominant proportion + body mass\n")
cat("==============================\n")
print(summary(m_prop_simple$fit))

cat("\n==============================\n")
cat("AIC\n")
cat("==============================\n")
print(AIC(m_class$fit, m_prop$fit, m_prop_simple$fit))

sink()

coef_all <- rbind(
  m_class$coef,
  m_prop$coef,
  m_prop_simple$coef
)

write.table(
  coef_all,
  paste0(out_prefix, ".coefficients.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

aic_tab <- data.frame(
  model = c(
    "High_Order_Diet",
    "Food_Proportion",
    "Plant_Proportion_Simple"
  ),
  AIC = c(
    AIC(m_class$fit),
    AIC(m_prop$fit),
    AIC(m_prop_simple$fit)
  ),
  n_species = c(
    nrow(m_class$data),
    nrow(m_prop$data),
    nrow(m_prop_simple$data)
  )
)

write.table(
  aic_tab,
  paste0(out_prefix, ".AIC.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  m_class$data,
  paste0(out_prefix, ".used_species.high_order_diet.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  m_prop$data,
  paste0(out_prefix, ".used_species.food_proportion.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

diet_count <- as.data.frame(table(dat_class$High_Order_Diet))
colnames(diet_count) <- c("High_Order_Diet", "n_species")

write.table(
  diet_count,
  paste0(out_prefix, ".high_order_diet.counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Done.\n")
cat("Summary:", paste0(out_prefix, ".model_summary.txt"), "\n")
cat("Coefficients:", paste0(out_prefix, ".coefficients.tsv"), "\n")
cat("AIC:", paste0(out_prefix, ".AIC.tsv"), "\n")
