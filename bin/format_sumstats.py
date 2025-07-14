import pandas as pd
import argparse
import sys

# Parse arguments
parser = argparse.ArgumentParser(description="Format harmonized summary stats for LDSC/LAVA")
parser.add_argument("input", type=str, help="Input file (harmonized summary stats)")
parser.add_argument("N_samples", type=str, help="Sample size to include")
parser.add_argument("output", type=str, help="Output file path")
args = parser.parse_args()

# Read harmonized summary stats
try:
    df = pd.read_csv(args.input, sep="\t", dtype=str)
except Exception as e:
    sys.exit(f"❌ Failed to read input file: {e}")

# Check for required columns
required_cols = ["variant_id", "effect_allele", "other_allele", "beta", "standard_error", "p_value"]
missing_cols = [col for col in required_cols if col not in df.columns]
if missing_cols:
    sys.exit(f"❌ Missing required columns: {', '.join(missing_cols)}")

# Convert numeric columns to float
df["beta"] = df["beta"].astype(float)
df["standard_error"] = df["standard_error"].astype(float)
df["p_value"] = df["p_value"].astype(float)

# Add N column
df["N"] = args.N_samples

# Calculate Z and select final columns
df["Z"] = df["beta"] / df["standard_error"]
df_formatted = df[["variant_id", "N", "Z", "effect_allele", "other_allele", "p_value"]]
df_formatted.columns = ["SNP", "N", "Z", "A1", "A2", "P"]

# Write to output
df_formatted.to_csv(args.output, sep="\t", index=False, quoting=3)  # quoting=3 means no quotes