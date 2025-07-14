# nf-genetic-correlations
Nextflow pipeline to perform genetic correlations with GWAS summary statistics (LDSC and LAVA)

This Nextflow pipeline receives as input GWAS summary statistics (at least two datasets - restricted to EUR genetic ancestry, for now) and processes them to compute global genetic correlations using LDSC (https://github.com/bulik/ldsc) and regional genetic correlations using LAVA (https://github.com/josefin-werme/LAVA).

-----------------------------------------------------------
DOCUMENTATION
-----------------------------------------------------------

**1. Install**
<git clone ape4fld/nf-genetic-correlations>

**2. Inputs needed from user:**
   
a) GWAS summary statistics files 
  - In text format (e.g. .tsv, .csv, .txt, etc.)
  - Required columns present (in any order, but with exact same column name): "variant_id", "effect_allele", "other_allele", "beta", "standard_error", "p_value"
    Note: variant_id are expected to be rsIDs; optimized for harmonized summary statistics downloaded from the GWAS catalog (https://www.ebi.ac.uk/gwas/).
  -  Other columns can be present and they will be ignored.
  -  Stored in /genetic_correlations/data/sumstats/

b) Metadata file
  - One file named 'metadata.txt' with the following column names (tab separated):
    "dataset", "filename", "N", "cases", "controls"
  - dataset = a short name for each dataset
  - filename = file name for each GWAS summary statistics file
  - N = total sample size; if it varies per variant one can use the max. N
  - Cases = number of cases (if continuous trait, set to NA)
  - Controls = number of controls (if continuous trait, set to NA)
  - Stored in /genetic_correlations/data/

