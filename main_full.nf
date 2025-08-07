nextflow.enable.dsl=2

def metadata_file = file(params.data_dir).resolve('metadata.txt')
def sampleSizeMap = metadata_file
    .readLines()
    .drop(1)
    .collectEntries { line ->
        def row = line.split('\t')
        [row[1], row[2] as Integer]
    }

Channel
    .fromPath("${params.data_dir}/sumstats/*")
    .map { file ->
        def sampleN = sampleSizeMap[file.getName()]
        if (sampleN == null)
            throw new RuntimeException("No sample size found for file: ${file.getName()}")

        // def suffixName = file.getName().replaceAll(/\.tsv\.gz$/, '')
        def suffixName = file.getBaseName(file.name.endsWith('.gz')? 2: 1)

        tuple(file, sampleN, suffixName)
    }
    .set { raw_sumstats_with_N }

// Create channel of munged sumstats directory
Channel
    .value(file("${params.output_dir}/munged"))
    .set { munged_dir_ch }

// Create channel of LDSC rg results directory
Channel
    .value(file("${params.output_dir}/ldsc_rg"))
    .set { rg_dir_ch }

process FormatSumstats {
    input:
    tuple path(file), val(sampleN), val(suffixName)

    output:
    tuple path("formatted_${suffixName}.tsv"), val(sampleN), val(suffixName)

    publishDir "${params.output_dir}/formatted", mode: 'copy'

    script:
    """
    Rscript ${params.bin_dir}/format_sumstats.R \
        ${file} \
        ${sampleN} \
        formatted_${suffixName}.tsv
    """
}

process MungeSumstats {
    label 'ldsc'

    input:
    tuple path(file), val(sampleN), val(suffixName)

    output:
    tuple path("${suffixName}.sumstats.gz"), val(suffixName)

    publishDir "${params.output_dir}/munged", mode: 'copy'

    script:
    """
    munge_sumstats.py \
        --sumstats ${file} \
        --N ${sampleN} \
        --merge-alleles ${params.hapmap_ref} \
        --chunksize 500000 \
        --out ${suffixName}
    """
}

process RunLDSC_h2 {
    label 'ldsc'

    input:
    tuple path(file), val(suffixName)

    output:
    path("ldsc_h2_${suffixName}.log")

    publishDir "${params.output_dir}/ldsc_h2", mode: 'copy'

    script:
    """
    ldsc.py \
        --h2 ${file} \
        --ref-ld-chr ${params.w_ld_chr}/ \
        --w-ld-chr ${params.w_ld_chr}/ \
        --out ldsc_h2_${suffixName}
    """
}

process PrepRg {
    input:
    path munged

    output:
    path("pairs_to_test.tsv")

    publishDir "${params.output_dir}/ldsc_rg", mode: 'copy'

    script:
    """
    Rscript ${params.bin_dir}/prep_pairs_for_rg.R \
        ${munged} \
        ${munged}
    """
}

process RunLDSC_rg {
    label 'ldsc'

    input:
    tuple path(file1), val(suffix1), path(file2), val(suffix2)

    output:
    path("ldsc_rg_${suffix1}_${suffix2}.log")

    publishDir "${params.output_dir}/ldsc_rg", mode: 'copy'

    script:
    """
    ldsc.py \
        --rg ${file1},${file2} \
        --ref-ld-chr ${params.w_ld_chr}/ \
        --w-ld-chr ${params.w_ld_chr}/ \
        --out ldsc_rg_${suffix1}_${suffix2}
    """
}

process PrepLAVA {

    input:
    path rg_dir

    output:
    path("*.txt"), emit: data_files
    path("all_rg_results.tsv"), emit: rg_all_results

    publishDir "${params.data_dir}/LAVA", mode: 'copy', pattern: "*.txt"
    publishDir "${params.output_dir}/ldsc_rg", mode: 'copy', pattern: "all_rg_results*.tsv"

    script:
    """
    Rscript ${params.bin_dir}/prep_for_LAVA.R \
        ${metadata_file} \
        ${params.output_dir}/formatted \
        ${rg_dir}
    """
}

process RunLAVA {
    label 'lava'

    input:
    path lava_data_dir

    output:
    path("*.rds")

    publishDir "${params.output_dir}/LAVA", mode: 'copy'

    script:
    """
    # Copy LAVA data files to working directory
    cp ${lava_data_dir}/*.txt .
    
    # Run LAVA with current directory as data location
    Rscript ${params.bin_dir}/lava.R \
        ${params.ref_ld_chr} \
        ${params.locus_file} \
        .
    """
}

workflow {
    
    // Step 1: Format
    formatted_sumstats = raw_sumstats_with_N | FormatSumstats

    // Step 2: Munge
    munged_sumstats = formatted_sumstats | MungeSumstats

    // Step 3: Estimate heritability (LDSC)
    munged_sumstats | RunLDSC_h2

    // Step 4. Prepare pairs of traits for RunLDSC_rg
    all_munged_files_ch = munged_sumstats
        .map { it[0] }  // extract only the file paths
        .collect()

    pairs_test = PrepRg(all_munged_files_ch)

    // Step 5: Run LDSC rg
    pairs_ch = pairs_test
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            tuple(file(row.file1), row.suffix1, file(row.file2), row.suffix2)
        }

    output_rg = pairs_ch | RunLDSC_rg

    // Step 6: Get sample overlap and info files for LAVA
    all_rg_results_ch = output_rg.collect()

    prep_lava = PrepLAVA(all_rg_results_ch)

    // Step 7: Run LAVA
    lava_data_ch = prep_lava.data_files.collect()

    RunLAVA(lava_data_ch)
}