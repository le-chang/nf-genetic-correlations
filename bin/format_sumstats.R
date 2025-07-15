# Load packages -----------------------------------------------------------
library(dplyr)
library(tidyr)
library(data.table)

# Arguments -----------------------------------------------------------
args <- commandArgs(TRUE)
input <- as.character(args[1])
N_samples <- as.character(args[2])
output <- as.character(args[3])

# Read harmonized summary stats -----------------------------------------------------------
df <- fread(input)

# Check for required columns -----------------------------------------------------------

cols_expected = c("variant_id", "effect_allele", "other_allele", "beta", "standard_error", "p_value")

check_colnames <- cols_expected %in% colnames(df)

if ("FALSE" %in% check_colnames == TRUE) {
  stop("The summary statistics do not have one of the expected column names. Please check that the input has the following column names (in no specific order):\n
      variant_id, effect_allele, other_allele, beta, standard_error, p_value.\n")
} else {
  df <- df %>%
    dplyr::select(
      variant_id,
      effect_allele,
      other_allele,
      beta,
      standard_error,
      p_value
    )
}

# if (ncol(df) == 6) {
#   colnames(df) <- c("variant_id", "effect_allele", "other_allele", "beta", "standard_error", "p_value")

# } else if (ncol(df) < 6) {
#   stop("The summary statistics file has less than 6 columns.\n Please check source data and make sure to have exactly 6 columns corresponding to:\n
#       variant ID, effect allele, other allele, effect size (beta), standard error of the effect size, p-value.\n
#       Order of columns needs to be respected, column names do not need to match, and separator type does not matter.\n")

# } else if (ncol(df) > 6) {
#   stop("The summary statistics file has more than 6 columns.\n Please check source data and make sure to have exactly 6 columns corresponding to:\n
#       variant ID, effect allele, other allele, effect size (beta), standard error of the effect size, p-value.\n
#       Order of columns needs to be respected, column names do not need to match, and separator type does not matter.\n")
# }

# Check if rsids are indeed rsids -----------------------------------------------------------
df <- df %>%
  filter(grepl("^rs", variant_id))

if (nrow(df) == 0) {
  warning("No rows with rsid starting with 'rs' found. No output will be generated.")
}

# Include sample size in sumstats -----------------------------------------------------------
df <- df %>% mutate(N = N_samples)

# Calculate Z-score and select required columns -----------------------------------------------------------
df_formatted <- df %>%
  mutate(Z = beta / standard_error) %>%
  transmute(
    SNP = variant_id,
    N = N,
    Z = Z,
    A1 = effect_allele,
    A2 = other_allele,
    P = p_value
  )

# Write to output -----------------------------------------------------------
if (nrow(df_formatted) > 0) {
  write.table(df_formatted, output, sep = "\t", row.names = F, quote = F)
}