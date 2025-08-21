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
    pattern = "*.log",  # Changed to .log files
    full.names = TRUE
  )

# Function to extract rg results from LDSC log files
extract_rg_from_log <- function(log_file) {
  lines <- readLines(log_file)
  
  # Find the line with genetic correlation results
  rg_line_idx <- grep("^Summary of Genetic Correlation Results", lines)
  
  if(length(rg_line_idx) == 0) {
    # Try to find the results table
    rg_line_idx <- grep("^p1\\s+p2\\s+rg\\s+se\\s+z\\s+p\\s+h2_obs", lines)
  }
  
  if(length(rg_line_idx) > 0) {
    # Read the results table (usually 2 lines after the header)
    result_line <- lines[rg_line_idx + 2]
    
    # Parse the result line
    result_parts <- strsplit(trimws(result_line), "\\s+")[[1]]
    
    if(length(result_parts) >= 12) {
      return(data.frame(
        p1 = result_parts[1],
        p2 = result_parts[2],
        rg = as.numeric(result_parts[3]),
        se = as.numeric(result_parts[4]),
        z = as.numeric(result_parts[5]),
        p = as.numeric(result_parts[6]),
        h2_obs = as.numeric(result_parts[7]),
        h2_obs_se = as.numeric(result_parts[8]),
        h2_int = as.numeric(result_parts[9]),
        h2_int_se = as.numeric(result_parts[10]),
        gcov_int = as.numeric(result_parts[11]),
        gcov_int_se = as.numeric(result_parts[12]),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(NULL)
}

# Extract results from all log files
files <- lapply(file_paths, extract_rg_from_log)
files <- files[!sapply(files, is.null)]

if(length(files) == 0) {
  stop("No valid LDSC rg results found in log files")
}

# Combine all results
all_rg <- rbindlist(files, fill = TRUE) %>%
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
  ) %>% 
  distinct(., p1, p2, .keep_all = TRUE)

phenotypes <- c(all_rg$p1, all_rg$p2) %>% unique()

###### Creating sample overlap matrix by extracting the intercept from LDSC results ######

n <- length(phenotypes)
covar_matrix <- matrix(NA, n, n)

rownames(covar_matrix) <- colnames(covar_matrix) <- phenotypes

# Fill the covariance matrix
for(i in 1:n) {
  for(j in 1:n) {
    
    pheno_i <- phenotypes[i]
    pheno_j <- phenotypes[j]
    
    cat("Getting genetic covariance intercept for", pheno_i, "and", pheno_j, ".\n")
    
    if (pheno_i == pheno_j) {
      gcov_int <- 1
    } else {
      subset <- dplyr::filter(all_rg, p1 == pheno_i & p2 == pheno_j)
      if (nrow(subset) > 0) {
        gcov_int <- subset[["gcov_int"]][1]  # Take first value
      } else {
        subset <- dplyr::filter(all_rg, p1 == pheno_j & p2 == pheno_i)
        if (nrow(subset) > 0) {
          gcov_int <- subset[["gcov_int"]][1]  # Take first value
        } else {
          cat("Warning: No LDSC output found for traits", pheno_i, "and", pheno_j, ". Setting to 0.\n")
          gcov_int <- 0  # Use 0 instead of NA for missing correlations
        }
      }
    }
    
    covar_matrix[pheno_i, pheno_j] <- gcov_int   
    cat("Done.\n")
  }
}

# Standardise the matrix
# Check for NA values before standardization
if(any(is.na(covar_matrix))) {
  warning("Covariance matrix contains NA values. Replacing with 0.")
  covar_matrix[is.na(covar_matrix)] <- 0
}

# Ensure matrix is symmetric
covar_matrix <- (covar_matrix + t(covar_matrix)) / 2

# Standardize
covar_matrix <-
  covar_matrix %>%
  cov2cor() %>%
  round(digits = 5)

# Save data ---------------------------------------------------------------

write.table(
  covar_matrix,
  file = "sample_overlap.txt",
  quote = F,
  row.names = TRUE,
  col.names = NA,
  sep = "\t"
)

write.table(
  all_rg,
  file = "all_rg_results.tsv",
  quote = F,
  row.names = F,
  sep = "\t"
)