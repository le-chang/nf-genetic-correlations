🧬 nf-genetic-correlations
Nextflow pipeline for global and regional genetic correlations using GWAS summary statistics
Supports LDSC for genome-wide correlations and LAVA for local (regional) genetic correlations.

📖 Overview
This pipeline processes harmonized GWAS summary statistics (restricted to European ancestry for now) and computes:

Global genetic correlations using LDSC

Local genetic correlations using LAVA

🚀 Getting Started
1. Install
bash
Copy
Edit
git clone https://github.com/ape4fld/nf-genetic-correlations.git
cd nf-genetic-correlations
2. Inputs Required
📁 a) GWAS Summary Statistics
Accepted formats: .tsv, .csv, .txt, etc.

Required columns (names must match exactly, but order can vary):

matlab
Copy
Edit
variant_id, effect_allele, other_allele, beta, standard_error, p_value
⚠️ variant_id must be rsIDs. This pipeline is optimized for harmonized summary stats from the GWAS Catalog.

Store files under:

swift
Copy
Edit
/genetic_correlations/data/sumstats/
📝 b) Metadata File
A single file named metadata.txt, tab-separated, with the following columns:

Column	Description
dataset	Short name for each dataset
filename	File name of the GWAS summary statistics file
N	Total sample size (use max if per-variant varies)
cases	Number of cases (use NA for continuous traits)
controls	Number of controls (use NA for continuous traits)

Store this file at:

bash
Copy
Edit
/genetic_correlations/data/
📦 c) LD Reference Files
LD Scores (for LDSC)
Download and extract:

bash
Copy
Edit
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
tar -xvjf eur_w_ld_chr.tar.bz2
Place in:

swift
Copy
Edit
/genetic_correlations/data/ld_reference/eur_w_ld_chr/
1000 Genomes (for LAVA)
Download PLINK files for European population as described in the LAVA reference guide

Place in:

swift
Copy
Edit
/genetic_correlations/data/ld_reference/g1000_eur/
3. ⚙️ Nextflow Configuration
This pipeline is configured to run on Alliance Canada’s Béluga cluster.
You will need to customize the provided nextflow.config:

Update the paths to your data files and directories.

Set the Apptainer image path for LDSC (e.g., .sif file).

Modify mounted directory paths for Apptainer.

Replace --account=def-xxxxx with your Alliance Canada user ID.

📂 Directory Structure
bash
Copy
Edit
nf-genetic-correlations/
├── main.nf
├── nextflow.config
├── data/
│   ├── sumstats/          # Input GWAS summary statistics
│   ├── metadata.txt       # Metadata file
│   └── ld_reference/
│       ├── eur_w_ld_chr/  # LD scores for LDSC
│       └── g1000_eur/     # 1000 Genomes data for LAVA
├── results/               # Output directory (after run)
└── README.md
