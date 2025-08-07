# 🧬 nf-genetic-correlations

**Nextflow pipeline for global and regional genetic correlations using GWAS summary statistics**  
Supports **LDSC** for genome-wide correlations and **LAVA** for local (regional) genetic correlations.
(forked from [ape4fld/nf-genetic-correlations](https://github.com/ape4fld/nf-genetic-correlations))

---

## 📖 Overview

This pipeline processes **harmonized GWAS summary statistics** (restricted to **European ancestry** for now) and computes:

- **Global genetic correlations** using [LDSC](https://github.com/bulik/ldsc)  
- **Local genetic correlations** using [LAVA](https://github.com/josefin-werme/LAVA)

It will also compute the SNP-based heritability for each of the GWAS summary statistics. For LAVA, it will perform the univariate test for each trait across all loci, and will compute the local genetic correlations (i.e., bivariate test) for pairs of traits which univariate test is significant (0.05/number of loci tested).

---

## 🚀 Getting Started

### 1. Install

```bash
git clone https://github.com/le-chang/nf-genetic-correlations.git
cd nf-genetic-correlations
```

### 2. Dependencies

Some R packages need to be pre-installed in R version 4.3.1:

- Tidyverse: dplyr, tidyr, stringr, readr
- Others: data.table
- LAVA (via remotes::install_github())

In Alliance Canada, you can follow these steps:

```bash
module load StdEnv/2023 r/4.3.1
mkdir -p ~/.local/R/$EBVERSIONR/
export R_LIBS=~/.local/R/$EBVERSIONR/
R -e 'install.packages(c("dplyr", "tidyr", "stringr", "readr", "data.table"), repos="https://cloud.r-project.org/")'
R -e 'remotes::install_github("josefin-werme/LAVA")'
```

### 3. Directory Structure Setup

The pipeline uses a flexible directory structure. By default, it expects this layout **relative to where you run the pipeline**:

```
your-project-directory/
├── bin/                    # Pipeline scripts (from git clone)
│   ├── format_sumstats.R
│   ├── lava.R
│   ├── prep_for_LAVA.R
│   ├── prep_pairs_for_rg.R
│   └── ldsc_latest.sif    # LDSC container (you'll download this)
├── data/                   # Your input data
│   ├── metadata.txt        # Your metadata file
│   ├── sumstats/          # Your GWAS summary statistics
│   │   ├── trait1.tsv.gz
│   │   ├── trait2.tsv.gz
│   │   └── ...
│   └── ld_reference/      # Reference files (you'll download these)
│       ├── eur_w_ld_chr/
│       ├── g1000_eur/
│       └── w_hm3.snplist
├── results/               # Output directory (created automatically)
├── main_full.nf          # Main pipeline file
├── nextflow.config       # Configuration file
└── run_nextflow.sh       # SLURM submission script (if using HPC)
```

### 4. Flexible Path Configuration

You have **three options** for setting up paths:

#### Option 1: Use Default Relative Paths (Recommended)
Simply organize your files according to the structure above and run from the project directory.

#### Option 2: Override Paths at Runtime
You can specify custom paths when running the pipeline:

```bash
nextflow run main_full.nf \
  --data_dir /path/to/your/data \
  --output_dir /path/to/your/results \
  --bin_dir /path/to/scripts
```

#### Option 3: Modify the Config File
Edit `nextflow.config` to set your permanent custom paths:

```groovy
params {
    // Set your custom paths here
    base_dir = "/home/username/my-genetics-project"
    data_dir = "/scratch/username/gwas_data"
    output_dir = "/scratch/username/gwas_results"
    bin_dir = "${base_dir}/bin"
    
    // Reference files paths
    ref_ld_chr = "${data_dir}/references/g1000_eur/g1000_eur"
    w_ld_chr = "${data_dir}/references/eur_w_ld_chr"
    hapmap_ref = "${data_dir}/references/w_hm3.snplist"
    locus_file = "${data_dir}/references/blocks_s2500_EUR_chr1-23full.GRCh37_hg19.locfile"
}
```

---

## 📁 Input Requirements

### a) GWAS Summary Statistics

- **Location**: Place in `data/sumstats/` (or your custom path)
- **Format**: `.tsv`, `.csv`, `.txt`, etc.
- **Required columns** (names must match exactly, order can vary):
  - `variant_id` (must be rsIDs)
  - `effect_allele`
  - `other_allele`
  - `beta`
  - `standard_error`
  - `p_value`

> ⚠️ This pipeline is optimized for harmonized summary stats from the [GWAS Catalog](https://www.ebi.ac.uk/gwas/).

### b) Metadata File

Create `metadata.txt` in your `data/` directory with these columns (tab-separated):

| Column     | Description                                      |
|------------|--------------------------------------------------|
| `dataset`  | Short name for each dataset                      |
| `filename` | File name of the GWAS summary statistics file    |
| `N`        | Total sample size                                |
| `cases`    | Number of cases (use `NA` for continuous traits) |
| `controls` | Number of controls (use `NA` for continuous traits)|

**Example `metadata.txt`:**
```
dataset	filename	N	cases	controls
height	height_gwas.tsv.gz	500000	NA	NA
t2d	t2d_gwas.tsv.gz	100000	25000	75000
bmi	bmi_gwas.tsv.gz	450000	NA	NA
```

### c) LD Reference Files

Download and place in `data/ld_reference/` (or your custom path):

1. **LD Scores for LDSC**  
```bash
cd data/ld_reference
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2
tar -xvjf eur_w_ld_chr.tar.bz2
```

2. **HapMap3 SNP list**
```bash
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
bunzip2 w_hm3.snplist.bz2
```

3. **1000 Genomes Reference for LAVA**
Download European PLINK files from [LAVA reference](https://github.com/josefin-werme/LAVA)

4. **Locus file for LAVA**
Download from LAVA GitHub repository

### d) LDSC Container

Download the LDSC Singularity/Apptainer image:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15920751.svg)](https://doi.org/10.5281/zenodo.15920751)

```bash
# Download to bin/ directory
wget -O bin/ldsc_latest.sif https://zenodo.org/records/15920751/files/ldsc_latest.sif
```

---

## 🚀 Running the Pipeline

### For Alliance Canada HPC:

1. **Edit the SLURM script** (`run_nextflow.sh`):
   - Replace `def-xxxxx` with your compute allocation
   - Adjust time/memory if needed

2. **Submit the job**:
```bash
sbatch run_nextflow.sh
```

### For Other Systems:

```bash
# With default paths
nextflow run main_full.nf -resume

# With custom paths
nextflow run main_full.nf \
  --data_dir /my/data/path \
  --output_dir /my/results/path \
  -resume

# For local execution (no HPC)
nextflow run main_full.nf -profile local -resume
```

### Using Different Profiles:

```bash
# Béluga cluster
nextflow run main_full.nf -profile beluga -resume

# Narval cluster  
nextflow run main_full.nf -profile narval -resume

# Local machine (for testing)
nextflow run main_full.nf -profile local -resume
```

---

## 📊 Expected Outputs

Results will be created in your `results/` directory (or custom output path):

```
results/
├── formatted/              # Formatted summary statistics
├── munged/                # LDSC-ready files
├── ldsc_h2/              # Heritability estimates
├── ldsc_rg/              # Global genetic correlations
│   ├── *.log             # Individual pairwise results
│   └── all_rg_results.tsv # Combined results table
└── LAVA/                  # Local genetic correlations
    ├── *.univ.lava.rds   # Univariate test results
    └── *.bivar.lava.rds  # Bivariate test results
```

---

## ⏱️ Performance Considerations

### Time Requirements:
- **3 traits**: ~5 hours (3 pairwise comparisons)
- **5 traits**: ~10 hours (10 pairwise comparisons)
- **10 traits**: ~30 hours (45 pairwise comparisons)

Adjust in `nextflow.config`:
```groovy
withLabel: lava {
    time = 47.h  // Increase for more traits
}
```

### Memory Requirements:
- Minimum: 10GB per process
- For >10 traits: Consider increasing to 20GB

---

## 🔧 Troubleshooting

### Common Issues:

1. **"No sample size found for file"**
   - Check that filenames in `metadata.txt` match exactly with files in `sumstats/`

2. **LDSC container not found**
   - Ensure `ldsc_latest.sif` is in the `bin/` directory
   - Check path in `nextflow.config`

3. **R package errors**
   - Verify all R packages are installed
   - Check R version (4.3.1 recommended)

4. **Path issues**
   - Use absolute paths if relative paths aren't working
   - Check that all reference files are downloaded

### Debug Mode:

Run with more verbose output:
```bash
nextflow run main_full.nf -with-trace -with-report -with-timeline
```

---

## 📧 Support

For issues or questions:
- Check the [GitHub Issues](https://github.com/ape4fld/nf-genetic-correlations/issues)
- Contact: [your-email@example.com]

---

## 📄 Citation

If you use this pipeline, please cite:
- LDSC: [Bulik-Sullivan et al., 2015](https://www.nature.com/articles/ng.3211)
- LAVA: [Werme et al., 2022](https://www.nature.com/articles/s41588-022-01082-3)
- This pipeline: [DOI/citation pending]