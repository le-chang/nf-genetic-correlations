# Load packages -----------------------------------------------------------
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)

# Arguments -----------------------------------------------------------
args <- commandArgs(TRUE)
metadata_file <- as.character(args[1])
munged_dir <- as.character(args[2])
run_id <- as.character(args[3]) # run identifier for output file prefix

# Set output file prefix based on run_id
file_prefix <- ifelse(run_id == "" || is.na(run_id), "", paste0(run_id, "_"))

# Read metadata and derive suffixes -------------------------------------------

metadata <- read.table(metadata_file, sep = "\t", header = TRUE)

# Derive suffix from filename (same logic as in main pipeline)
suffixes <- metadata$filename %>%
    stringr::str_remove(".gz$") %>%
    stringr::str_remove(".[:alpha:]+$")

# Get all pairwise combinations of suffixes
pairwise_combos <- combn(suffixes, 2)

# Convert to data frame with file paths and suffixes
df_out <- data.frame(
    suffix1 = pairwise_combos[1, ],
    suffix2 = pairwise_combos[2, ]
) %>%
    mutate(
        file1 = stringr::str_c(munged_dir, "/", suffix1, ".sumstats.gz"),
        file2 = stringr::str_c(munged_dir, "/", suffix2, ".sumstats.gz")
    ) %>%
    dplyr::select(file1, suffix1, file2, suffix2)

write.table(df_out, paste0(file_prefix, "pairs_to_test.tsv"),
            sep = "\t", row.names = F, quote = F, col.names = T)