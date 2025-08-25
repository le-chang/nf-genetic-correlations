#!/bin/bash

echo "=== Fixing broken symlinks for LDSC genetic correlation analysis ==="
echo ""

# Navigate to main project directory
cd /lustre06/project/6054517/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/

echo "1. Current directory: $(pwd)"
echo ""

echo "2. Checking what's in the root directory:"
ls -la *.sumstats.gz 2>/dev/null || echo "No .sumstats.gz files in root directory (expected)"
echo ""

echo "3. Checking what's in results/munged/:"
ls -la results/munged/*.sumstats.gz
echo ""

echo "4. Creating symlinks in root directory to point to correct files:"
ln -sf results/munged/AD_kunkle2019_formatted.sumstats.gz ./AD_kunkle2019_formatted.sumstats.gz
ln -sf results/munged/uacr_female_formatted.sumstats.gz ./uacr_female_formatted.sumstats.gz
ln -sf results/munged/uacr_male_formatted.sumstats.gz ./uacr_male_formatted.sumstats.gz
ln -sf results/munged/uacr_sex_combined_formatted.sumstats.gz ./uacr_sex_combined_formatted.sumstats.gz

echo "Created symlinks:"
ls -la *.sumstats.gz
echo ""

echo "5. Verifying symlinks work:"
for file in *.sumstats.gz; do
    if [ -f "$file" ]; then
        echo "✓ $file -> $(readlink $file) [$(ls -lh $file | awk '{print $5}')]"
    else
        echo "✗ $file (broken)"
    fi
done
echo ""

echo "6. You can now resume the pipeline with:"
echo "   sbatch run_nextflow.sh"
echo "   OR"
echo "   nextflow run main_full.nf -profile beluga -resume"
