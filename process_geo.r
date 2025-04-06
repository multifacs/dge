#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

# Check arguments
if (length(args) < 2) {
  stop("Usage: Rscript script.R <input_file> <output_name>", call. = FALSE)
}

input_file <- args[1]
name <- args[2]
cat(sprintf("Processing '%s'\n", input_file))

pacman::p_load(GEOquery, dplyr, tidyr, stringr)

# Load GEO data with platform annotations (GPL)
gse <- getGEO(filename = input_file, getGPL = TRUE)

# Extract expression matrix and feature data (probe annotations)
expression_matrix <- exprs(gse)
feature_data <- fData(gse)  # Probe-to-gene mappings from GPL

# Find the gene symbol column (handle naming variations)
symbol_col <- grep("gene.*symbol", colnames(feature_data), ignore.case = TRUE, value = TRUE)[1]
if (is.na(symbol_col)) {
  stop("No 'Gene Symbol' column found in the platform annotations.")
}

# Create probe-to-gene mapping, filter ambiguous/empty symbols, and drop duplicates
gene_mapping <- data.frame(
  ProbeID = rownames(feature_data),
  Gene_Symbol = feature_data[[symbol_col]],
  stringsAsFactors = FALSE
) %>%
  filter(
    !str_detect(Gene_Symbol, "///"),  # Remove ambiguous probes
    Gene_Symbol != ""                 # Remove empty symbols
  ) %>%
  distinct(Gene_Symbol, .keep_all = TRUE)  # Keep only first occurrence of each gene

# Replace ProbeIDs with Gene Symbols and filter expression matrix
expression_df <- as.data.frame(expression_matrix) %>%
  tibble::rownames_to_column("ProbeID") %>%
  inner_join(gene_mapping, by = "ProbeID") %>%  # Keep only filtered probes
  select(-ProbeID) %>%
  tibble::column_to_rownames("Gene_Symbol")

# Save cleaned data
dir_path <- file.path("data", name)
if (!dir.exists(dir_path)) {
  dir.create(dir_path, recursive = TRUE)
}

expr_path <- file.path(dir_path, paste0(name, "_expr.csv"))
phen_path <- file.path(dir_path, paste0(name, "_phen.csv"))

write.csv(expression_df, expr_path)
write.csv(pData(gse), phen_path)

cat(sprintf("Cleaned expression matrix saved to: %s\n", expr_path))