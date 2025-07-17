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
        ldsc munge_sumstats \
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
    path("ldsc_h2_${suffixName}.h2_results")

    publishDir "${params.output_dir}/ldsc_h2", mode: 'copy'

    script:
    """
        ldsc ldsc h2 \
            --h2 ${file} \
            --ref-ld-chr ${params.w_ld_chr} \
            --w-ld-chr ${params.w_ld_chr} \
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
        ${params.output_dir}
    """
}

process RunLDSC_rg {
    label 'ldsc'

    input:
    tuple path(file1), val(suffix1), path(file2), val(suffix2)

    output:
    path("ldsc_rg_${suffix1}_${suffix2}.rg_results")

    publishDir "${params.output_dir}/ldsc_rg", mode: 'copy'

    script:
    """
        ldsc ldsc rg \
            --rg ${file1} \
            --rg ${file2} \
            --ref-ld-chr ${params.w_ld_chr} \
            --w-ld-chr ${params.w_ld_chr} \
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
    // params.ref_ld_chr = LD reference plink files
    // loc = locus file for LAVA
    // info and sample overlap files
    """
        Rscript ${params.bin_dir}/lava.R \
            ${params.ref_ld_chr} \
            ${params.locus_file} \
            ${lava_data_dir}
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
        .collect()
        .map { it[0] }  // extract only the file paths

    pairs_test = all_munged_files_ch
        .combine(munged_dir_ch)
        .map { munged_files, munged_dir -> munged_dir }
        | PrepRg

    // Step 5: Run LDSC rg
    pairs_ch = pairs_test
        .map { file -> file.readLines().drop(1) }  // skip header
        .flatten()
        .map { line -> 
            def row = line.split('\t')
            tuple(row[0], row[1], row[2], row[3])
        }

    output_rg = pairs_ch | RunLDSC_rg

    // Step 6: Get sample overlap and info files for LAVA
    all_rg_results_ch = output_rg
    .collect()
    .map { it[0] }

    lava_input_ch = all_rg_results_ch
        .combine(rg_dir_ch)
        .map { results, rg_dir -> rg_dir }

    prep_lava = PrepLAVA(lava_input_ch)

    // Step 7: Run LAVA
    lava_data_ch = prep_lava.data_files
        .collect()
        .map { _ -> "${params.data_dir}/LAVA" }

    lava_data_ch | RunLAVA
}
