#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(nlme)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("
Usage:
Rscript 03_run_pgls_body_mass.R <pgls_input.tsv> <tree.nwk> <out_prefix>

Example:
Rscript 03_run_pgls_body_mass.R \\
  301.PGA_diet_bodymass.best_assembly.tsv \\
  mammal_time_tree.nwk \\
  PGA_diet_bodymass_PGLS
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

# Make sure species names match tree tip labels
dat$tip_label <- gsub(" ", "_", dat$tip_label)

# Basic variables
dat$PGA_log1p <- log1p(dat$PGA_CN_primary)
dat$log10_body_mass <- log10(dat$`BodyMass-Value`)
dat$log10_contig_N50 <- log10(dat$`contig N50 (bp)`)

dat$plant_total_prop <- dat$plant_total_pct / 100
dat$plant_structural_prop <- dat$plant_structural_prop

dat$Diet_Class <- as.factor(dat$Diet_Class)
dat$Diet_Group <- as.factor(dat$Diet_Group)

# Set reference level
if ("Animal_dominant" %in% levels(dat$Diet_Group)) {
  dat$Diet_Group <- relevel(dat$Diet_Group, ref = "Animal_dominant")
}

# Keep species shared by tree and data
common <- intersect(tree$tip.label, dat$tip_label)

cat("Species in data:", length(unique(dat$tip_label)), "\n")
cat("Tips in tree:", length(tree$tip.label), "\n")
cat("Matched species:", length(common), "\n")

if (length(common) < 50) {
  warning("Matched species < 50. Please check tree tip labels and data species names.")
}

tree <- drop.tip(tree, setdiff(tree$tip.label, common))
dat <- dat[dat$tip_label %in% common, ]

# Function to run PGLS with Pagel's lambda
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
    method = "ML"
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

# Main models
m1 <- fit_pgls(
  PGA_log1p ~ plant_total_prop + log10_body_mass,
  dat,
  tree,
  "PGA_log1p ~ plant_total_prop + log10_body_mass"
)

m2 <- fit_pgls(
  PGA_log1p ~ plant_structural_prop + log10_body_mass,
  dat,
  tree,
  "PGA_log1p ~ plant_structural_prop + log10_body_mass"
)

m3 <- fit_pgls(
  PGA_log1p ~ Diet_Group + log10_body_mass,
  dat,
  tree,
  "PGA_log1p ~ Diet_Group + log10_body_mass"
)

m4 <- fit_pgls(
  PGA_log1p ~ plant_total_prop + log10_body_mass + log10_contig_N50,
  dat,
  tree,
  "PGA_log1p ~ plant_total_prop + log10_body_mass + log10_contig_N50"
)

# Save model summaries
sink(paste0(out_prefix, ".model_summary.txt"))

cat("\n==============================\n")
cat("Model 1: plant_total_prop + body mass\n")
cat("==============================\n")
print(summary(m1$fit))

cat("\n==============================\n")
cat("Model 2: plant_structural_prop + body mass\n")
cat("==============================\n")
print(summary(m2$fit))

cat("\n==============================\n")
cat("Model 3: Diet_Group + body mass\n")
cat("==============================\n")
print(summary(m3$fit))

cat("\n==============================\n")
cat("Model 4: plant_total_prop + body mass + contig N50\n")
cat("==============================\n")
print(summary(m4$fit))

cat("\n==============================\n")
cat("AIC\n")
cat("==============================\n")
print(AIC(m1$fit, m2$fit, m3$fit, m4$fit))

sink()

# Save coefficient table
coef_all <- rbind(
  m1$coef,
  m2$coef,
  m3$coef,
  m4$coef
)

write.table(
  coef_all,
  paste0(out_prefix, ".coefficients.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Save matched data used in primary model
write.table(
  m1$data,
  paste0(out_prefix, ".primary_model.used_species.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Done.\n")
cat("Summary:", paste0(out_prefix, ".model_summary.txt"), "\n")
cat("Coefficients:", paste0(out_prefix, ".coefficients.tsv"), "\n")
cat("Used species:", paste0(out_prefix, ".primary_model.used_species.tsv"), "\n")
