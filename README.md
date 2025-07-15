# 🧬 nf-genetic-correlations

**Nextflow pipeline for global and regional genetic correlations using GWAS summary statistics**  
Supports **LDSC** for genome-wide correlations and **LAVA** for local (regional) genetic correlations.

---

## 📖 Overview

This pipeline processes **harmonized GWAS summary statistics** (restricted to **European ancestry** for now) and computes:

- **Global genetic correlations** using [LDSC](https://github.com/bulik/ldsc)  
- **Local genetic correlations** using [LAVA](https://github.com/josefin-werme/LAVA)

It will also compute the SNP-based heritability for each of the GWAS summary statistics. For LAVA, it will perform the univariate test for each trait across all loci, and will compute the local genetic correlations (i.e., bivariate test) for pairs of traits which univariate test is significant (0.05/number of loci tested).

There are several advantages of using the pipeline:
1) Given that it uses an LDSC .sif image, there is no need to load old python versions to run LDSC.
2) The pipeline formats and adapts the GWAS summary statistics for each tool.
3) The user does not need to prepare additional files to run LAVA.
4) It is reproducible and the user can easily re-run the analysis by adding/removing GWAS datasets.

---

## 🚀 Getting Started

### 1. Install

```bash
git clone https://github.com/ape4fld/nf-genetic-correlations.git
cd nf-genetic-correlations
```

### 2. Dependencies

Some R packages need to be pre-installed in R version 4.3.1:

- Tidyverse: dplyr, tidyr, stringr, readr
- Others: here, data.table
- LAVA (```R via remotes::install_github()```)

In Alliance Canada, you can follow these steps:

```bash
module load StdEnv/2023 r/4.3.1
mkdir -p ~/.local/R/$EBVERSIONR/
export R_LIBS=~/.local/R/$EBVERSIONR/
R -e 'install.packages(c("dplyr", "tidyr", "stringr", "readr", "here", "data.table"), repos="https://cloud.r-project.org/")'
R -e 'remotes::install_github("josefin-werme/LAVA")'
```

### 3. Inputs Required

---

#### 📁 a) GWAS Summary Statistics

- Accepted formats: `.tsv`, `.csv`, `.txt`, etc.
- Required columns (**names must match exactly**, order can vary):
variant_id, effect_allele, other_allele, beta, standard_error, p_value

> ⚠️ `variant_id` must be rsIDs. This pipeline is optimized for harmonized summary stats from the [GWAS Catalog](https://www.ebi.ac.uk/gwas/).

- Store your files in:
```bash
/genetic_correlations/data/sumstats/
```

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
```

#### 📦 c) LD Reference Files

1. **LD Scores for LDSC**  
 Download and extract:

 ```bash
 wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
 tar -xvjf eur_w_ld_chr.tar.bz2
```

Place contents in:
 ```bash
/genetic_correlations/data/ld_reference/eur_w_ld_chr/
```

2. **1000 Genomes Reference (for LAVA)**
Download European PLINK reference files as described in the LAVA reference guide

Place contents in:
 ```bash
/genetic_correlations/data/ld_reference/g1000_eur/
```

### 4. 📦 LDSC Apptainer/Singularity Image

---

Download the LDSC container image from Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15920751.svg)](https://doi.org/10.5281/zenodo.15920751)

```bash
# Download the image (1.2GB)
wget https://zenodo.org/records/15920751/files/ldsc_latest.sif
# Place it in the bin/ directory
mv ldsc_latest.sif bin/
```

### 5. ⚙️ Nextflow Configuration

---

This pipeline is configured to run on **Alliance Canada's Béluga cluster**, but can be adapted to other HPC environments.

You will need to **edit the provided `nextflow.config`** file to match your setup:

- 🔧 **Update parameters** to match paths to your input and reference files.
- 📦 **Set the path to the Apptainer image** (`.sif`) used to run LDSC.
- 📁 **Define the directories** that will be mounted into the Apptainer container.
- 🧑‍💻 **Replace your user account string** in the SLURM options:  
  Change  
  ```nextflow
  process.clusterOptions = '--account=def-xxxxx'
  ```

---

## 🚀 Running the Pipeline

Once you've completed the setup and configuration, you can run the pipeline:

### For Alliance Canada/Béluga Users:

1. **Edit the SLURM script** (`run_nextflow.sh`):
   - Replace `def-xxxxx` with your compute allocation
   - Update file paths to match your directory structure

2. **Submit the job**:
   ```bash
   sbatch run_nextflow.sh
   ```

### For Other HPC/Local Systems:

Run Nextflow directly:
```bash
nextflow run main_full.nf -profile <your_profile> -resume
```

The pipeline will:
- Process your GWAS summary statistics
- Calculate global genetic correlations using LDSC
- Calculate local genetic correlations using LAVA
- Output results to the `results/` directory

---

## 📊 Expected Outputs

The pipeline generates results in the following directory structure:

```
results/
├── formatted/                 # Formatted summary statistics
│   └── formatted_*.tsv        # One file per GWAS dataset
├── munged/                    # LDSC-ready files
│   └── *.sumstats.gz         # Munged summary statistics
├── ldsc_h2/                   # Heritability estimates
│   └── *.h2_results          # SNP-heritability for each trait
├── ldsc_rg/                   # Global genetic correlations
│   ├── *.rg_results          # Pairwise genetic correlations
│   └── all_rg_results.tsv    # Combined results table
└── LAVA/                      # Local genetic correlations
    ├── univ_*.rds            # Univariate test results per trait
    └── bivar_*.rds           # Bivariate test results for trait pairs

data/LAVA/                     # LAVA input files
├── info_file.txt             # Trait information
└── sample_overlap.txt        # Sample overlap matrix
```

### Key Output Files:

- **`all_rg_results.tsv`**: Summary table with all global genetic correlations (rg), standard errors, p-values, and heritability estimates
- **`univ_*.rds`**: Local heritability and association p-values for each genomic locus per trait
- **`bivar_*.rds`**: Local genetic correlations between trait pairs at specific loci where both traits show significant univariate signals
