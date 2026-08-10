suppressPackageStartupMessages({
  library(PheWAS)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3 || length(args) > 4) {
  stop(
    "Usage: Rscript PheWAS.R <cohort.tsv> <icd10.tsv> <output.csv> [cores]"
  )
}

cohort_file <- args[1]
icd10_file <- args[2]
output_file <- args[3]
cores <- if (length(args) == 4) as.integer(args[4]) else 1

covariates <- read.delim(cohort_file, check.names = FALSE)
icd10_data <- read.delim(icd10_file, check.names = FALSE)

required_covariates <- c(
  "eid", "Pred_PGA34A", "Age", "Sex", paste0("PC", 1:10)
)
missing_covariates <- setdiff(required_covariates, names(covariates))
if (length(missing_covariates) > 0) {
  stop(
    "Missing cohort columns: ",
    paste(missing_covariates, collapse = ", ")
  )
}

required_icd10 <- c("id", "vocabulary_id", "code", "count")
missing_icd10 <- setdiff(required_icd10, names(icd10_data))
if (length(missing_icd10) > 0) {
  stop(
    "Missing ICD-10 columns: ",
    paste(missing_icd10, collapse = ", ")
  )
}

covariates <- covariates[, required_covariates]
names(covariates)[names(covariates) == "eid"] <- "id"
icd10_data <- icd10_data[, required_icd10]

normalize_ukb_icd10 <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("\\.", "", x)
  ifelse(
    nchar(x) > 3,
    paste0(substr(x, 1, 3), ".", substr(x, 4, nchar(x))),
    x
  )
}

icd10_data$code <- normalize_ukb_icd10(icd10_data$code)

id_sex <- data.frame(
  id = covariates$id,
  sex = ifelse(covariates$Sex == 1, "M", "F")
)

phenotypes <- createPhenotypes(
  id.vocab.code.index = icd10_data,
  min.code.count = 1,
  add.phecode.exclusions = TRUE,
  translate = TRUE,
  vocabulary.map = phecode_map_icd10,
  id.sex = id_sex,
  full.population.ids = covariates$id
)

if (!setequal(phenotypes$id, covariates$id)) {
  stop("PheCode matrix does not contain the complete eligible cohort.")
}

data_for_phewas <- inner_join(phenotypes, covariates, by = "id")

results <- phewas(
  phenotypes = names(phenotypes)[-1],
  genotypes = "Pred_PGA34A",
  data = data_for_phewas,
  covariates = c("Age", "Sex", paste0("PC", 1:10)),
  significance.threshold = "bonferroni",
  additive.genotypes = FALSE,
  min.records = 20,
  cores = cores
)

n_tests <- attr(results, "n.tests")
results$p_bonferroni <- p.adjust(
  results$p,
  method = "bonferroni",
  n = n_tests
)
results <- addPhecodeInfo(results)

write.csv(results, output_file, row.names = FALSE)
