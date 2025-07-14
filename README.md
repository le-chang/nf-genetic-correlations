# 🧬 nf-genetic-correlations

**Nextflow pipeline for global and regional genetic correlations using GWAS summary statistics**  
Supports **LDSC** for genome-wide correlations and **LAVA** for local (regional) genetic correlations.

---

## 📖 Overview

This pipeline processes **harmonized GWAS summary statistics** (restricted to **European ancestry** for now) and computes:

- **Global genetic correlations** using [LDSC](https://github.com/bulik/ldsc)  
- **Local genetic correlations** using [LAVA](https://github.com/josefin-werme/LAVA)

---

## 🚀 Getting Started

### 1. Install

```bash
git clone https://github.com/ape4fld/nf-genetic-correlations.git
cd nf-genetic-correlations

### 2. Inputs Required

---

#### 📁 a) GWAS Summary Statistics

- Accepted formats: `.tsv`, `.csv`, `.txt`, etc.
- Required columns (**names must match exactly**, order can vary):
variant_id, effect_allele, other_allele, beta, standard_error, p_value

> ⚠️ `variant_id` must be rsIDs. This pipeline is optimized for harmonized summary stats from the [GWAS Catalog](https://www.ebi.ac.uk/gwas/).

- Store your files in:
/genetic_correlations/data/sumstats/

#### 📝 b) Metadata File

Create a single file named `metadata.txt`, tab-separated, with the following columns:

| Column     | Description                                      |
|------------|--------------------------------------------------|
| `dataset`  | Short name for each dataset                      |
| `filename` | File name of the GWAS summary statistics file    |
| `N`        | Total sample size (use max if per-variant varies)|
| `cases`    | Number of cases (use `NA` for continuous traits) |
| `controls` | Number of controls (use `NA` for continuous traits)|

- Store the metadata file at:
 ```bash
/genetic_correlations/data/

#### 📦 c) LD Reference Files

1. **LD Scores for LDSC**  
 Download and extract:

 ```bash
 wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
 tar -xvjf eur_w_ld_chr.tar.bz2

Place contents in:
 ```bash
/genetic_correlations/data/ld_reference/eur_w_ld_chr/

2. **1000 Genomes Reference (for LAVA)**
Download European PLINK reference files as described in the LAVA reference guide

Place contents in:
 ```bash
/genetic_correlations/data/ld_reference/g1000_eur/

### 3. ⚙️ Nextflow Configuration

---

This pipeline is configured to run on **Alliance Canada’s Béluga cluster**, but can be adapted to other HPC environments.

You will need to **edit the provided `nextflow.config`** file to match your setup:

- 🔧 **Update parameters** to match paths to your input and reference files.
- 📦 **Set the path to the Apptainer image** (`.sif`) used to run LDSC.
- 📁 **Define the directories** that will be mounted into the Apptainer container.
- 🧑‍💻 **Replace your user account string** in the SLURM options:  
  Change  
  ```nextflow
  process.clusterOptions = '--account=def-xxxxx'
