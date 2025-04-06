#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Check arguments
if (length(args) < 2) {
    stop("Usage: Rscript script.R <input.csv> <output.csv>", call. = FALSE)
}

input_file <- args[1]
name <- args[2]
cat(sprintf("Processed '%s'\n", input_file))

pacman::p_load(GEOquery, openxlsx)

gse <- getGEO(filename = input_file, getGPL = FALSE)

expression_matrix <- exprs(gse)
phenotype_data <- pData(gse)

dim(expression_matrix) # Genes x Samples
dim(phenotype_data) # Samples x Phenotype variables

dir_path <- paste("data\\", name, "\\", sep = "")
expr_path <- file.path(dir_path, paste(name, "_expr.csv", sep = ""))
phen_path <- file.path(dir_path, paste(name, "_phen.csv", sep = ""))

if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE, showWarnings = TRUE)
}

write.csv(expression_matrix, expr_path)
write.csv(phenotype_data, phen_path)
# Print confirmation
cat(sprintf("Processed '%s'\n", input_file))
