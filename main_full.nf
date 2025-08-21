nextflow.enable.dsl=2

// Check if metadata file exists
def metadata_file = file(params.data_dir).resolve('metadata.txt')
if (!metadata_file.exists()) {
    error "Metadata file not found: ${metadata_file}"
}

// Parse metadata to get sample sizes
def sampleSizeMap = metadata_file
    .readLines()
    .drop(1)  // Skip header
    .collectEntries { line ->
        def row = line.split('\t')
        if (row.size() < 3) {
            error "Invalid metadata line: ${line}"
        }
        // Handle both integer N and "NA" values
        def n_value = row[2].trim()
        def sample_n = n_value == "NA" ? null : n_value as Integer
        [row[1].trim(), sample_n]
    }

// Create input channel with proper error handling
Channel
    .fromPath("${params.data_dir}/sumstats/*")
    .map { file ->
        def sampleN = sampleSizeMap[file.getName()]
        if (sampleN == null) {
            error "No sample size found for file: ${file.getName()}. Check metadata.txt"
        }

        // Handle different file extensions properly
        def suffixName = file.getBaseName(file.name.endsWith('.gz') ? 2 : 1)
        
        log.info "Processing: ${file.getName()} with N=${sampleN}, suffix=${suffixName}"
        
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
    tag "${suffixName}"
    
    input:
    tuple path(file), val(sampleN), val(suffixName)

    output:
    tuple path("formatted_${suffixName}.tsv"), val(sampleN), val(suffixName)

    publishDir "${params.output_dir}/formatted", mode: 'copy'

    script:
    """
    echo "Formatting ${file} with N=${sampleN}"
    
    Rscript ${params.bin_dir}/format_sumstats.R \
        ${file} \
        ${sampleN} \
        formatted_${suffixName}.tsv
    
    # Check if output was created
    if [ ! -f formatted_${suffixName}.tsv ]; then
        echo "Error: formatted file was not created"
        exit 1
    fi
    
    echo "Created formatted_${suffixName}.tsv with \$(wc -l < formatted_${suffixName}.tsv) lines"
    """
}

process MungeSumstats {
    tag "${suffixName}"
    label 'ldsc'

    input:
    tuple path(file), val(sampleN), val(suffixName)

    output:
    tuple path("${suffixName}.sumstats.gz"), val(suffixName)

    publishDir "${params.output_dir}/munged", mode: 'copy'

    script:
    """
    echo "Munging ${file} with N=${sampleN}"
    
    munge_sumstats.py \
        --sumstats ${file} \
        --N ${sampleN} \
        --merge-alleles ${params.hapmap_ref} \
        --chunksize 500000 \
        --out ${suffixName}
    
    # Check output
    if [ ! -f ${suffixName}.sumstats.gz ]; then
        echo "Error: munged file was not created"
        exit 1
    fi
    """
}

process RunLDSC_h2 {
    tag "${suffixName}"
    label 'ldsc'

    input:
    tuple path(file), val(suffixName)

    output:
    path("ldsc_h2_${suffixName}.log")

    publishDir "${params.output_dir}/ldsc_h2", mode: 'copy'

    script:
    """
    echo "Running LDSC h2 for ${suffixName}"
    
    ldsc.py \
        --h2 ${file} \
        --ref-ld-chr ${params.w_ld_chr}/ \
        --w-ld-chr ${params.w_ld_chr}/ \
        --out ldsc_h2_${suffixName}
    
    # Check output
    if [ ! -f ldsc_h2_${suffixName}.log ]; then
        echo "Error: LDSC h2 log was not created"
        exit 1
    fi
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
    echo "Preparing pairs for genetic correlation"
    
    # Check if we have munged files
    n_files=\$(ls -1 *.sumstats.gz 2>/dev/null | wc -l)
    echo "Found \$n_files munged files"
    
    if [ "\$n_files" -lt 2 ]; then
        echo "Error: Need at least 2 munged files for correlation analysis"
        exit 1
    fi
    
    Rscript ${params.bin_dir}/prep_pairs_for_rg.R \
        . \
        .
    
    # Check output
    if [ ! -f pairs_to_test.tsv ]; then
        echo "Error: pairs file was not created"
        exit 1
    fi
    
    echo "Created pairs file with \$(wc -l < pairs_to_test.tsv) pairs"
    """
}

process RunLDSC_rg {
    tag "${suffix1}_${suffix2}"
    label 'ldsc'

    input:
    tuple path(file1), val(suffix1), path(file2), val(suffix2)

    output:
    path("ldsc_rg_${suffix1}_${suffix2}.log")

    publishDir "${params.output_dir}/ldsc_rg", mode: 'copy'

    script:
    """
    echo "Running LDSC rg for ${suffix1} vs ${suffix2}"
    
    ldsc.py \
        --rg ${file1},${file2} \
        --ref-ld-chr ${params.w_ld_chr}/ \
        --w-ld-chr ${params.w_ld_chr}/ \
        --out ldsc_rg_${suffix1}_${suffix2}
    
    # Check output and show summary
    if [ -f ldsc_rg_${suffix1}_${suffix2}.log ]; then
        echo "Genetic correlation computed successfully"
        grep "Genetic Correlation:" ldsc_rg_${suffix1}_${suffix2}.log || true
    else
        echo "Error: LDSC rg log was not created"
        exit 1
    fi
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
    echo "Preparing LAVA input files"
    
    # Check if we have rg results
    n_logs=\$(ls -1 ldsc_rg_*.log 2>/dev/null | wc -l)
    echo "Found \$n_logs LDSC rg log files"
    
    if [ "\$n_logs" -eq 0 ]; then
        echo "Warning: No LDSC rg results found"
    fi
    
    Rscript ${params.bin_dir}/prep_for_LAVA.R \
        ${metadata_file} \
        ${params.output_dir}/formatted \
        .
    
    # Check outputs
    if [ ! -f info_file.txt ]; then
        echo "Error: info_file.txt was not created"
        exit 1
    fi
    
    if [ ! -f sample_overlap.txt ]; then
        echo "Error: sample_overlap.txt was not created"
        exit 1
    fi
    
    echo "Created LAVA input files successfully"
    """
}

process RunLAVA {
    tag "LAVA_analysis"
    label 'lava'

    input:
    path lava_data_dir

    output:
    path("*.rds")

    publishDir "${params.output_dir}/LAVA", mode: 'copy'

    script:
    """
    echo "Running LAVA analysis"
    
    # Copy LAVA data files to working directory
    cp ${lava_data_dir}/*.txt .
    
    # Check that we have the required files
    if [ ! -f info_file.txt ]; then
        echo "Error: info_file.txt not found"
        exit 1
    fi
    
    if [ ! -f sample_overlap.txt ]; then
        echo "Error: sample_overlap.txt not found"
        exit 1
    fi
    
    # Check that reference files exist
    if [ ! -d ${params.ref_ld_chr} ]; then
        echo "Error: LAVA reference directory not found: ${params.ref_ld_chr}"
        echo "Please extract the LAVA reference files first"
        exit 1
    fi
    
    # Run LAVA with current directory as data location
    Rscript ${params.bin_dir}/lava.R \
        ${params.ref_ld_chr} \
        ${params.locus_file} \
        .
    
    # Check outputs
    n_rds=\$(ls -1 *.rds 2>/dev/null | wc -l)
    if [ "\$n_rds" -eq 0 ]; then
        echo "Warning: No RDS files created by LAVA"
    else
        echo "LAVA created \$n_rds output files"
    fi
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

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    println "Duration: $workflow.duration"
}