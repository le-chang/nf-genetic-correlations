# Load packages -----------------------------------------------------------
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)

# Arguments -----------------------------------------------------------
args <- commandArgs(TRUE)
munged_path <- as.character(args[1])
munged_dir <- as.character(args[2])

# Read munged files -----------------------------------------------------------

munged_list <- list.files(stringr::str_c(munged_path), pattern = "*.sumstats.gz", full.names = TRUE)

# Get all pairwise combinations
pairwise_combos <- combn(munged_list, 2)

# Convert to data frame with one trait per column
pairwise_df <- data.frame(
    file1v1 = pairwise_combos[1, ],
    file2v1 = pairwise_combos[2, ]
)

df_out <- pairwise_df %>%
    mutate(file1 = stringr::str_c(munged_dir, "/", file1v1),
           suffix1 = stringr::str_remove(basename(file1), ".sumstats.gz"),
           file2 = stringr::str_c(munged_dir, "/", file2v1),
           suffix2 = stringr::str_remove(basename(file2), ".sumstats.gz")) %>%
    dplyr::select(file1, suffix1, file2, suffix2)

write.table(df_out, "pairs_to_test.tsv", sep = "\t", row.names = F, quote = F, col.names = T)