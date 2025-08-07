# Description: 
# 1. Obtain the sample overlaps between pairs of traits for LAVA.
# 2. Obtain the info file for LAVA from the metadata file.

# Load packages -----------------------------------------------------------

library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(readr)

# Set arguments -----------------------------------------------------------

args <- commandArgs(TRUE)
metadata_file <- as.character(args[1]) # metadata path
formatted_dir <- as.character(args[2]) # directory for formatted summary statistics
ldsc_rg <- as.character(args[3]) # results dir from LDSC rg

# -----------------------------------------------------------
# -----------------------------------------------------------
# Format metadata file into the info file for LAVA:

metadata <- read.table(metadata_file, sep = "\t", header = TRUE) 

info <- metadata %>%
  mutate(suffix = stringr::str_remove(filename, ".gz$") %>% stringr::str_remove(., ".[:alpha:]+$"),
        formatted_path = stringr::str_c(formatted_dir, "/formatted_", suffix, ".tsv")) %>%
  dplyr::select(phenotype = suffix, cases, controls, filename = formatted_path)
write.table(info, "info_file.txt", sep = "\t", row.names = F, quote = F)

# -----------------------------------------------------------
# -----------------------------------------------------------

###### Extracting LDSC results ######

# Extract all outputs into single file
file_paths <-
  list.files(
    ldsc_rg,
    pattern = "*.rg_results",
    full.names = TRUE
  )

files <-
  vector(mode = "list",
         length = length(file_paths))

for(i in 1:length(file_paths)){
  
  files[[i]] <-
    read_table(
      file = file_paths[i]
    )
  
}

# extract rg values in matrix
all_rg <-
  files %>%
  rbindlist(., fill = TRUE) %>%
  dplyr::mutate(
    p1 = basename(p1) %>%
      stringr::str_remove(".sumstats.gz"),
    p2 = basename(p2) %>%
      stringr::str_remove(".sumstats.gz")
  ) %>%
  select(
    p1,
    p2,
    rg,
    se,
    z,
    p,
    h2_obs,
    h2_obs_se,
    h2_int,
    h2_int_se,
    gcov_int,
    gcov_int_se
  ) %>% distinct(., p1, p2, .keep_all = TRUE)

phenotypes <- c(all_rg$p1, all_rg$p2) %>% unique()

###### Creating sample overlap matrix by extracting the intercept from LDSC results ######

n <- length(phenotypes)
covar_matrix <- matrix(NA, n, n)

rownames(covar_matrix) <- colnames(covar_matrix) <- phenotypes

# FIX: Remove the unnecessary outer k loop
for(i in phenotypes) {
  for(j in phenotypes) {
    
    cat("Getting genetic covariance intercept for", i, "and", j, ".\n")
    
    if (i == j) {
      gcov_int <- 1
    } else {
      subset <- dplyr::filter(all_rg, p1 == i, p2 == j)
      if (nrow(subset) > 0) {
        gcov_int <- subset[["gcov_int"]]
      } else {
        subset <- dplyr::filter(all_rg, p1 == j, p2 == i)
        if (nrow(subset) > 0) {
          gcov_int <- subset[["gcov_int"]]
        } else {
          cat("Warning: No LDSC output found for traits", i, "and", j, ". Setting to NA.\n")
          gcov_int <- NA
          next
        }
      }
    }
    
    covar_matrix[i,j] <- gcov_int   
    cat("Done.\n")
  }
}

# Standardise the matrix (moved outside the loop)
# Check for NA values before standardization
if(any(is.na(covar_matrix))) {
  warning("Covariance matrix contains NA values. These will be handled in standardization.")
  # Option 1: Replace NA with 0 for missing correlations
  covar_matrix[is.na(covar_matrix)] <- 0
}

covar_matrix <-
  covar_matrix %>%
  cov2cor() %>%
  round(digits = 5)

# Save data ---------------------------------------------------------------

write.table(
  covar_matrix,
  file = "sample_overlap.txt",
  quote = F,
  row.names = TRUE,  # Changed to TRUE to preserve row names
  col.names = NA,     # This ensures proper format with row names
  sep = "\t"
)

write.table(
  all_rg,
  file = "all_rg_results.tsv",
  quote = F,
  row.names = F,
  sep = "\t"
)