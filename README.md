# nf-genetic-correlations
Nextflow pipeline to perform genetic correlations with GWAS summary statistics (LDSC and LAVA)

This Nextflow pipeline receives as input GWAS summary statistics (at least two datasets - restricted to EUR genetic ancestry, for now) and processes them to compute global genetic correlations using LDSC (https://github.com/bulik/ldsc) and regional genetic correlations using LAVA (https://github.com/josefin-werme/LAVA).


