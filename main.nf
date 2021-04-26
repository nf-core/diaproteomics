#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/diaproteomics
========================================================================================
 nf-core/diaproteomics Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/diaproteomics
----------------------------------------------------------------------------------------
*/

log.info Headers.nf_core(workflow, params.monochrome_logs)

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////+
def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/diaproteomics --input sample_sheet.tsv --input_spectral_library library_sheet.tsv --irts irt_sheet.tsv"
    log.info NfcoreSchema.params_help(workflow, params, json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         VALIDATE PARAMETERS              -- */
////////////////////////////////////////////////////+
if (params.validate_params) {
    NfcoreSchema.validateParameters(params, json_schema, log)
}

////////////////////////////////////////////////////
/* --     Collect configuration parameters     -- */
////////////////////////////////////////////////////

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, 'Specify correct --awsqueue and --awsregion parameters on AWSBatch!'
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, 'Outdir not on S3 - specify S3 Bucket to run on AWSBatch!'
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, 'Specify a local tracedir or run without trace! S3 cannot be used for tracefiles.'
}

// Stage config files
projectDir=workflow.projectDir
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)
sample_sheet = file(params.input)
Channel
    .from( sample_sheet )
    .into { input_exp_design; input_exp_design_mztab}

params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()

// DIA MS input
Channel
    .from( sample_sheet )
    .splitCsv(header: true, sep:'\t')
    .map { col -> tuple("${col.Sample}", "${col.BatchID}", "${col.MSstats_Condition}", file("${col.Spectra_Filepath}", checkifExists: true))}
    .flatMap{it -> [tuple(it[0],it[1].toString(),it[2],it[3])]}
    .into { input_branch; check_dia }


// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

input_branch.branch {
    raw: hasExtension(it[3], 'raw')
    mzml: hasExtension(it[3], 'mzml')
    mzxml: hasExtension(it[3], 'mzxml')
    other: true
}.set{input_dia_ms_files}

input_dia_ms_files.other.subscribe { row -> log.warn("unknown format for entry " + row[3] + " in provided sample sheet. ignoring line."); exit 1 }
check_dia_n = check_dia.map{it[1]}.unique().toList().size().val

// Validate inputs
sample_sheet = file(params.input)

if( params.irts_from_outer_quantiles){
    quant_flag = '--quantiles True'
} else {
    quant_flag = ''
}

if( params.align_libraries) {
    align_flag = '--align True'
} else {
    align_flag = ''
}

/*
 * Create a channel for input spectral library
 */
if( params.generate_spectral_library) {

    // Spectral library input
    dda_sheet = file(params.input_sheet_dda)

    Channel.from( dda_sheet )
        .splitCsv(header: true, sep:'\t')
        .map { col -> tuple("${col.Sample}", "${col.BatchID}", file("${col.Spectra_Filepath}", checkifExists: true), file("${col.Id_Filepath}", checkifExists: true))}
        .flatMap{it -> [tuple(it[0],it[1],it[2],it[3])]}
        .into {input_dda;input_check;input_check_samples}

    check_n = input_check.toList().size().val
    check_n_sample = input_check_samples.map{it[1]}.unique().toList().size().val
    if ((check_n > 1) & (check_n != check_n_sample) & (!params.merge_libraries)) {
        print('You specified multiple DDA files to generate spectral libraries, but library merging is not set \n')
        print('Set --merge_libraries and possibly --align_libraries to align them in the same RT space \n')
        exit 1
    }

    if ((check_n_sample != check_dia_n)) {
        print('The number of batches in the sample input does not match the number of batches of the spectral library input \n')
        exit 1
    }

    input_dda.branch {
        raw: hasExtension(it[2], 'raw')
        mzml: hasExtension(it[2], 'mzML')
        mzxml: hasExtension(it[2], 'mzXML')
        other: true
    }.set{input_dda_ms_files}

    Channel
        .fromPath( params.unimod )
        .ifEmpty { exit 1, "params.unimod was empty - no unimod.xml supplied" }
        .set { input_unimod }

    input_lib = Channel.empty()
    input_lib_1 = Channel.empty()
    input_lib_nd = Channel.empty()
    input_lib_nd_1 = Channel.empty()
    params.input_spectral_library = "generate spectral library from DDA data"

} else if( params.skip_decoy_generation) {

    // Spectral library input
    library_sheet = file(params.input_spectral_library)

    Channel.from( library_sheet )
        .splitCsv(header: true, sep:'\t')
        .map { col -> tuple("${col.Sample}", "${col.BatchID}", file("${col.Library_Filepath}", checkifExists: true))}
        .flatMap{ it -> [tuple(it[0],it[1],it[2])] }
        .into { input_lib_nd; check_decoy; check_decoy_2 }

    check_decoy_n2 = check_decoy_2.toList().size().val
    check_decoy_n = check_decoy.map{it[1]}.unique().toList().size().val
    if ((check_decoy_n > 1) & (check_decoy_n2 != check_decoy_n)) {
        print('You specified multiple spectral libraries including decoys for one batch \n')
        print('This is not possible yet, merging and aligning would be possible without decoys\n')
        exit 1
    }

    if ((check_decoy_n != check_dia_n)) {
        print('The number of batches in the sample input do not match the number of batches of the spectral library input \n')
        exit 1
    }

    input_lib = Channel.empty()
    input_lib_1 = Channel.empty()
    input_dda = Channel.empty()
    input_dda.branch {
        raw: hasExtension(it[2], 'raw')
        mzml: hasExtension(it[2], 'mzML')
        mzxml: hasExtension(it[2], 'mzXML')
        other: true
    }.set { input_dda_ms_files }
    input_unimod = Channel.empty()

} else {

    // Spectral library input
    library_sheet = file(params.input_spectral_library)

    Channel.from( library_sheet )
        .splitCsv(header: true, sep:'\t')
        .map { col -> tuple("${col.Sample}", "${col.BatchID}", file("${col.Library_Filepath}", checkifExists: true))}
        .flatMap{it -> [tuple(it[0],it[1],it[2])]}
        .into {input_lib; input_lib_1; input_lib_n; input_lib_n2 }

    check_lib_n2 = input_lib_n2.toList().size().val
    check_lib_n = input_lib_n.map{it[1]}.unique().toList().size().val
    if ((check_lib_n > 1) & (check_lib_n2 != check_lib_n) & (!params.merge_libraries)) {
        print('You specified multiple spectral libraries for one batch \n')
        print('Set --merge_libraries and possibly --align_libraries to align them in the same RT space \n')
        exit 1
    }

    if ((check_lib_n != check_dia_n)) {
        print('The number of batches in the sample input does not match the number of batches of the spectral library input \n')
        exit 1
    }

    input_lib_nd = Channel.empty()
    input_lib_nd_1 = Channel.empty()
    input_dda = Channel.empty()
    input_dda.branch {
        raw: hasExtension(it[2], 'raw')
        mzml: hasExtension(it[2], 'mzML')
        mzxml: hasExtension(it[2], 'mzXML')
        other: true
    }.set { input_dda_ms_files }
    input_unimod = Channel.empty()
}


if( !params.generate_pseudo_irts){
    // iRT library input
    irt_sheet = file(params.irts)

    Channel.from( irt_sheet )
        .splitCsv(header: true, sep:'\t')
        .map { col -> tuple("${col.BatchID}", file("${col.irt_Filepath}", checkifExists: true))}
        .flatMap{it -> [tuple(it[0],it[1])]}
        .into {input_irts; input_irts_check; input_irts_check_2}

    check_irts_n = input_irts_check.toList().size().val
    check_irts_n2 = input_irts_check_2.map{it[0]}.unique().toList().size().val
    if ((check_irts_n > 1) & (check_irts_n != check_irts_n2)) {
        print('You specified multiple DDA files to generate spectral libraries, but library merging is not set \n')
        print('Set --merge_libraries and possibly --align_libraries to align them in the same RT space \n')
        exit 1
    }

    if ((check_irts_n2 != check_dia_n)) {
        print('The number of batches in the sample input does not match the number of batches of the irt input \n')
        exit 1
    }

} else {

    input_irts = Channel.empty()

}

// MS1 option
if (params.use_ms1){
    ms1_option = '-use_ms1_traces'
    ms1_scoring = '-Scoring:Scores:use_ms1_mi'
    ms1_mi = '-Scoring:Scores:use_mi_score'
} else {
    ms1_option = ''
    ms1_scoring = ''
    ms1_mi = ''
}

// Force option
if (params.force_option){
    force_option = '-force'
} else {
    force_option = ''
}

// DIAlignR multithreading
if (params.dialignr_parallelization){
    dialignr_parallel='parallel'
} else {
    dialignr_parallel=''
}

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////
log.info NfcoreSchema.params_summary_log(workflow, params, json_schema)

// Header log info
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Spectral Library']    = params.input_spectral_library
summary['Run Name']         = workflow.runName
summary['Input']            = params.input
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
if (workflow.containerEngine) summary['Container'] = "$workflow.containerEngine - $workflow.container"
summary['Output dir']       = params.outdir
summary['Launch dir']       = workflow.launchDir
summary['Working dir']      = workflow.workDir
summary['Script dir']       = workflow.projectDir
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
}
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Profile Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Profile Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config Profile URL']         = params.config_profile_url
summary['Config Files'] = workflow.configFiles.join(', ')
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
}

// Check the hostnames against configured profiles
checkHostname()

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-diaproteomics-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/diaproteomics Workflow Summary'
    section_href: 'https://github.com/nf-core/diaproteomics'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf('.csv') > 0) filename
                      else null
        }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file 'software_versions.csv'

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    FileInfo --help &> v_openms.txt
    pyprophet --version &> v_pyprophet.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}


/*
 * STEP 0 - Raw File Conversion
 */
process dda_raw_file_conversion {
    input:
    set val(id), val(sample), file(raw_file), file(dda_id_file) from input_dda_ms_files.raw

    output:
    set val(id), val(sample), file("${raw_file.baseName}.mzML"), file(dda_id_file) into converted_dda_input_mzmls

    when:
    params.generate_spectral_library

    script:
    """
    ThermoRawFileParser.sh -i=${raw_file} -f=2 -b=${raw_file.baseName}.mzML
    """
}


/*
 * STEP 1 - Convert IDs for Spectral Library Generation using EasyPQP
 */
process dda_id_format_conversion {

    input:
    set val(id), val(sample), file(dda_mzml), file(dda_id_file) from input_dda_ms_files.mzml.mix(input_dda_ms_files.mzxml).mix(converted_dda_input_mzmls)

    output:
    set val(id), val(sample), file(dda_mzml), file("${id}_${sample}_peptide_ids.idXML") into input_dda_converted

    when:
    params.generate_spectral_library

    script:
    """
    IDFileConverter -in ${dda_id_file} -out ${id}_${sample}_peptide_ids.idXML -threads ${task.cpus}
    """
}


/*
 * STEP 2 - Spectral Library Generation using EasyPQP
 */
process dda_library_generation {

    input:
    set val(id), val(sample), file(dda_mzml_file), file(idxml_file) from input_dda_converted
    file unimod_file from input_unimod.first()

    output:
    set val(id), val(sample), file("${id}_${sample}_library.tsv") into input_lib_dda_nd

    when:
    params.generate_spectral_library

    script:
    """
    easypqp convert \\
        --unimod ${unimod_file} \\
        --pepxml ${idxml_file} \\
        --spectra ${dda_mzml_file}

    easypqp library \\
        --out ${dda_mzml_file.baseName}_run_peaks.tsv \\
        --rt_psm_fdr_threshold ${params.library_rt_fdr} \\
        --nofdr \\
        ${dda_mzml_file.baseName}.psmpkl \\
        ${dda_mzml_file.baseName}.peakpkl

    mv ${dda_mzml_file.baseName}_run_peaks.tsv ${id}_${sample}_library.tsv
    """
}


/*
 * STEP 3 - Assay Generation for Spectral Library
 */
process assay_generation {

    input:
    set val(id), val(sample), file(lib_file_na) from input_lib.mix(input_lib_dda_nd)

    output:
    set val(id), val(sample), file("${id}_${sample}_assay.tsv") into (input_lib_assay, input_lib_assay_for_irt, input_lib_assay_for_merging)

    when:
    !params.skip_decoy_generation

    script:
    """
    TargetedFileConverter \\
        -in ${lib_file_na} \\
        -out ${lib_file_na.baseName}.tsv \\
        -threads ${task.cpus}

    OpenSwathAssayGenerator \\
        -in ${lib_file_na.baseName}.tsv \\
        -min_transitions ${params.min_transitions} \\
        -max_transitions ${params.max_transitions} \\
        -out ${id}_${sample}_assay.tsv \\
        -threads ${task.cpus}
    """
}


if(params.merge_libraries) {
    input_lib_assay = Channel.empty()
    input_lib_assay_for_irt = Channel.empty()
}


/*
 * STEP 4 - Merge and align spectral Libraries
 */
process library_merging_and_alignment {
    publishDir "${params.outdir}/spectral_library_files"

    input:
    set val(id), val(sample), file(lib_files_for_merging) from input_lib_assay_for_merging.groupTuple(by:1)

    output:
    set val(id), val(sample), file("${sample}_library_merged.tsv") into (input_lib_assay_merged, input_lib_assay_merged_for_irt)
    set val(id), val(sample), file("*.png") optional true

    when:
    params.merge_libraries

    script:
    """
    merge_and_align_libraries_from_easypqp.py \\
        --input_libraries ${lib_files_for_merging} \\
        --min_overlap ${params.min_overlap_for_merging} \\
        --rsq_threshold 0.75  \\
        --output ${sample}_library_merged.tsv \\
        ${align_flag}
    """
}


/*
 * STEP 5 - Pseudo iRT Library Generation
 */
process pseudo_irt_generation {
    publishDir "${params.outdir}/spectral_library_files"

    input:
    set val(id), val(sample), file(lib_file_assay_irt) from input_lib_assay_for_irt.mix(input_lib_assay_merged_for_irt)

    output:
    set val(sample), file("${lib_file_assay_irt.baseName}_pseudo_irts.pqp") into input_lib_assay_irt_2

    when:
    params.generate_pseudo_irts

    script:
    """
    select_pseudo_irts_from_lib.py \\
        --input_libraries ${lib_file_assay_irt} \\
        --min_rt 0 \\
        --n_irts ${params.n_irts} \\
        --max_rt 100 \\
        --output ${lib_file_assay_irt.baseName}_pseudo_irts.tsv \\
        ${quant_flag}

    TargetedFileConverter \\
        -in ${lib_file_assay_irt.baseName}_pseudo_irts.tsv \\
        -out ${lib_file_assay_irt.baseName}_pseudo_irts.pqp \\
        -threads ${task.cpus}
    """
}


/*
 * STEP 6 - Decoy Generation for Spectral Library
 */
process decoy_generation {
    publishDir "${params.outdir}/spectral_library_files"

    input:
    set val(id), val(sample), file(lib_file_nd) from input_lib_assay.mix(input_lib_assay_merged)

    output:
    set val(id), val(sample), file("${lib_file_nd.baseName}_decoy.pqp") into input_lib_decoy

    when:
    !params.skip_decoy_generation

    script:
    """
    TargetedFileConverter \\
        -in ${lib_file_nd} \\
        -out ${lib_file_nd.baseName}.pqp \\
        -threads ${task.cpus}

    OpenSwathDecoyGenerator \\
        -in ${lib_file_nd.baseName}.pqp \\
        -method ${params.decoy_method} \\
        -out ${lib_file_nd.baseName}_decoy.pqp \\
        -threads ${task.cpus}
    """
}


/*
 * STEP 7 - DIA Raw File Conversion
 */
process dia_raw_file_conversion {

    input:
    set val(id), val(sample), val(condition), file(raw_file) from input_dia_ms_files.raw

    output:
    set val(id), val(sample), val(condition), file("${raw_file.baseName}.mzML") into converted_dia_input_mzmls

    when:
    !params.skip_dia_processing

    script:
    """
    ThermoRawFileParser.sh -i=${raw_file} -f=2 -b=${raw_file.baseName}.mzML
    """
}


/*
 * STEP 8 - DIA library search with OpenSwathWorkFlow
 */
process dia_spectral_library_search {
    publishDir "${params.outdir}/openswathworkflow_output"

    label 'process_medium'

    input:
    set val(sample), val(id), val(condition), file(mzml_file), val(dummy_id), file(lib_file), file(irt_file) from converted_dia_input_mzmls.mix(input_dia_ms_files.mzml.mix(input_dia_ms_files.mzxml)).combine(input_lib_decoy.mix(input_lib_nd), by:1).combine(input_irts.mix(input_lib_assay_irt_2), by:0)

    output:
    set val(id), val(sample), val(condition), file("${mzml_file.baseName}_chrom.mzML") into chromatogram_files
    set val(id), val(sample), val(condition), file("${mzml_file.baseName}.osw") into osw_files
    set val(id), val(sample), file("${lib_file.baseName}.pqp") into (input_lib_used, input_lib_used_I, input_lib_used_I_mztab)

    when:
    !params.skip_dia_processing


    script:
    """
    mkdir tmp

    TargetedFileConverter \\
        -in ${lib_file} \\
        -out ${lib_file.baseName}.pqp \\
        -threads ${task.cpus}

    TargetedFileConverter \\
        -in ${irt_file} \\
        -out ${irt_file.baseName}.pqp \\
        -threads ${task.cpus}

    OpenSwathWorkflow \\
        -in ${mzml_file} \\
        -tr ${lib_file.baseName}.pqp \\
        -sort_swath_maps \\
        -tr_irt ${irt_file.baseName}.pqp \\
        -min_rsq ${params.irt_min_rsq} \\
        -out_osw ${mzml_file.baseName}.osw \\
        -out_chrom ${mzml_file.baseName}_chrom.mzML \\
        -mz_extraction_window ${params.mz_extraction_window} \\
        -mz_extraction_window_ms1 ${params.mz_extraction_window_ms1} \\
        -mz_extraction_window_unit ${params.mz_extraction_window_unit} \\
        -mz_extraction_window_ms1_unit ${params.mz_extraction_window_ms1_unit} \\
        -rt_extraction_window ${params.rt_extraction_window} \\
        -min_upper_edge_dist ${params.min_upper_edge_dist} \\
        -RTNormalization:alignmentMethod ${params.irt_alignment_method} \\
        -RTNormalization:estimateBestPeptides \\
        -RTNormalization:outlierMethod none \\
        -RTNormalization:NrRTBins ${params.irt_n_bins} \\
        -RTNormalization:MinBinsFilled ${params.irt_min_bins_covered} \\
        -mz_correction_function quadratic_regression_delta_ppm \\
        -Scoring:stop_report_after_feature 5 \\
        -Scoring:TransitionGroupPicker:compute_peak_quality false \\
        -Scoring:TransitionGroupPicker:peak_integration 'original' \\
        -Scoring:TransitionGroupPicker:background_subtraction 'none' \\
        -Scoring:TransitionGroupPicker:PeakPickerMRM:sgolay_frame_length 11 \\
        -Scoring:TransitionGroupPicker:PeakPickerMRM:sgolay_polynomial_order 3 \\
        -Scoring:TransitionGroupPicker:PeakPickerMRM:gauss_width 30 \\
        -Scoring:TransitionGroupPicker:PeakPickerMRM:use_gauss 'false' \\
        -Scoring:TransitionGroupPicker:PeakIntegrator:integration_type 'intensity_sum' \\
        -Scoring:TransitionGroupPicker:PeakIntegrator:baseline_type 'base_to_base' \\
        -Scoring:TransitionGroupPicker:PeakIntegrator:fit_EMG 'false' \\
        -batchSize 1000 \\
        -readOptions ${params.cache_option} \\
        -tempDirectory tmp \\
        -Scoring:DIAScoring:dia_nr_isotopes 3 \\
        -enable_uis_scoring \\
        -Scoring:uis_threshold_sn -1 \\
        -threads ${task.cpus} \\
        ${force_option} ${ms1_option} ${ms1_scoring} ${ms1_mi}
    """
}


/*
 * STEP 9 - Pyprophet merging of OpenSwath results
 */
process dia_search_output_merging {

    input:
    set val(sample), val(id), val(condition), file(all_osws), val(dummy_id), file(lib_file_template) from osw_files.groupTuple(by:1).join(input_lib_used, by:1)

    output:
    set val(id), val(sample), val(condition), file("${sample}_osw_file_merged.osw") into merged_osw_file_for_global

    when:
    !params.skip_dia_processing

    script:
    """
    pyprophet merge \\
        --template=${lib_file_template} \\
        --out=${sample}_osw_file_merged.osw \\
        --no-same_run \\
        ${all_osws}
    """
}


/*
 * STEP 10 - Pyprophet global FDR Scoring
 */
process global_false_discovery_rate_estimation {
    publishDir "${params.outdir}/pyprophet_output"

    label 'process_medium'

    input:
    set val(id), val(sample), val(condition), file(scored_osw) from merged_osw_file_for_global

    output:
    set val(id), val(sample), val(condition), file("${scored_osw.baseName}_global_merged.osw") into merged_osw_scored_global_for_pyprophet
    set val(id), val(sample), val(condition), file("*.pdf") into target_decoy_global_score_plots

    when:
    !params.skip_dia_processing

    script:
    if (params.pyprophet_classifier=='LDA'){
        """
        pyprophet score \\
            --in=${scored_osw} \\
            --level=${params.pyprophet_fdr_ms_level} \\
            --out=${scored_osw.baseName}_scored.osw \\
            --classifier=${params.pyprophet_classifier} \\
            --pi0_lambda ${params.pyprophet_pi0_start} ${params.pyprophet_pi0_end} ${params.pyprophet_pi0_steps} \\
            --threads=${task.cpus}

        pyprophet peptide \\
            --in=${scored_osw.baseName}_scored.osw \\
            --out=${scored_osw.baseName}_global_merged.osw \\
            --context=run-specific

        pyprophet peptide --in=${scored_osw.baseName}_global_merged.osw --context=experiment-wide

        pyprophet peptide --in=${scored_osw.baseName}_global_merged.osw --context=global

        pyprophet ${params.pyprophet_global_fdr_level} --in=${scored_osw.baseName}_global_merged.osw --context=run-specific

        pyprophet ${params.pyprophet_global_fdr_level} --in=${scored_osw.baseName}_global_merged.osw --context=experiment-wide

        pyprophet ${params.pyprophet_global_fdr_level} --in=${scored_osw.baseName}_global_merged.osw --context=global
        """
    } else {
        """
        pyprophet score \\
            --in=${scored_osw} \\
            --level=${params.pyprophet_fdr_ms_level} \\
            --out=${scored_osw.baseName}_scored.osw \\
            --classifier=${params.pyprophet_classifier} \\
            --threads=${task.cpus}

        pyprophet peptide \\
            --in=${scored_osw.baseName}_scored.osw \\
            --out=${scored_osw.baseName}_global_merged.osw \\
            --context=run-specific

        pyprophet peptide --in=${scored_osw.baseName}_global_merged.osw --context=experiment-wide

        pyprophet peptide --in=${scored_osw.baseName}_global_merged.osw --context=global

        pyprophet ${params.pyprophet_global_fdr_level} --in=${scored_osw.baseName}_global_merged.osw --context=run-specific

        pyprophet ${params.pyprophet_global_fdr_level} --in=${scored_osw.baseName}_global_merged.osw --context=experiment-wide

        pyprophet ${params.pyprophet_global_fdr_level} --in=${scored_osw.baseName}_global_merged.osw --context=global
        """
    }
}


/*
 * STEP 11 - Pyprophet Export
 */
process export_of_scoring_results {
    publishDir "${params.outdir}/pyprophet_output"

    input:
    set val(id), val(sample), val(condition), file(global_osw) from merged_osw_scored_global_for_pyprophet

    output:
    set val(id), val(sample), val(condition), file("*.tsv") into pyprophet_results
    set val(id), val(sample), val(condition), file(global_osw) into osw_for_dialignr

    when:
    !params.skip_dia_processing

    script:
    """
    pyprophet export \\
        --in=${global_osw} \\
        --max_rs_peakgroup_qvalue=${params.pyprophet_peakgroup_fdr} \\
        --max_global_peptide_qvalue=${params.pyprophet_peptide_fdr} \\
        --max_global_protein_qvalue=${params.pyprophet_protein_fdr} \\
        --out=legacy.tsv
    """
}


/*
 * STEP 12 - Index Chromatogram mzMLs
 */
process chromatogram_indexing {

    label 'process_high'

    input:
    set val(id), val(sample), val(condition), file(chrom_file_noindex) from chromatogram_files

    output:
    set val(id), val(sample), val(condition), file("${chrom_file_noindex.baseName.split('_chrom')[0]}.chrom.sqMass") into chromatogram_files_indexed

    when:
    !params.skip_dia_processing

    script:
    """
    FileConverter \\
        -in ${chrom_file_noindex} \\
        -process_lowmemory \\
        -out ${chrom_file_noindex.baseName.split('_chrom')[0]}.chrom.mzML

    OpenSwathMzMLFileCacher \\
        -in ${chrom_file_noindex.baseName.split('_chrom')[0]}.chrom.mzML \\
        -lossy_compression false \\
        -process_lowmemory \\
        -lowmem_batchsize 50000 \\
        -out ${chrom_file_noindex.baseName.split('_chrom')[0]}.chrom.sqMass
    """
}


// Combine channels of osw files and osw chromatograms
osw_for_dialignr
    .transpose()
    .join(chromatogram_files_indexed, by:1)
    .groupTuple(by:0)
    // Channel contains now the following elements:
    // ([id, [samples], [conditions], [osw_files], id_2, condition_2, [chromatogram_files]])
    .flatMap{it -> [tuple(it[0],it[1].unique()[0],it[2].unique()[0],it[3].unique()[0],it[4],it[5],it[6])]}
    .set{osw_and_chromatograms_combined_by_condition}

/*
 * STEP 13 - Align DIA Chromatograms using DIAlignR
 */
process chromatogram_alignment {
    publishDir "${params.outdir}/"

    label 'process_high_mem'

    input:
    set val(sample), val(id), val(condition), file(pyresults), val(id_dummy), val(condition_dummy), file(chrom_files_index) from osw_and_chromatograms_combined_by_condition

    output:
    set val(id), val(sample), val(condition), file("${sample}_peptide_quantities.csv") into (DIALignR_result, DIALignR_result_I, DIALignR_result_mztab)

    when:
    !params.skip_dia_processing

    script:
    """
    mkdir osw
    mv ${pyresults} osw/
    mkdir xics
    mv *.chrom.sqMass xics/

    DIAlignR.R \\
        ${params.dialignr_global_align_fdr} \\
        ${params.dialignr_analyte_fdr} \\
        ${params.dialignr_unalign_fdr} \\
        ${params.dialignr_align_fdr} \\
        ${params.dialignr_query_fdr} \\
        ${params.pyprophet_global_fdr_level} \\
        ${params.dialignr_xicfilter} \\
        ${dialignr_parallel} \\
        ${task.cpus}

    mv DIAlignR.tsv ${sample}_peptide_quantities.csv
    """
}


/*
 * STEP 14 - Reformat output for MSstats: Combine with experimental design and missing columns from input library
 */
process reformatting {
    publishDir "${params.outdir}/"

    input:
    set val(id), val(sample), val(condition), file(dialignr_file) from DIALignR_result
    file exp_design from input_exp_design.first()
    set val(id), val(sample_lib), file(lib_file) from input_lib_used_I.first()

    output:
    set val(id), val(sample), val(condition), file("${sample}_${condition}.csv") into msstats_file

    when:
    params.run_msstats

    script:

    if (params.pyprophet_global_fdr_level==''){

    """
    TargetedFileConverter \\
        -in ${lib_file} \\
        -out ${lib_file.baseName}.tsv

    reformat_output_for_msstats.py \\
        --input ${dialignr_file} \\
        --exp_design ${exp_design} \\
        --library ${lib_file.baseName}.tsv \\
        --fdr_level "none" \\
        --output "${sample}_${condition}.csv"
    """

    } else {

    """
    TargetedFileConverter -in ${lib_file} -out ${lib_file.baseName}.tsv

    reformat_output_for_msstats.py \\
        --input ${dialignr_file} \\
        --exp_design ${exp_design} \\
        --library ${lib_file.baseName}.tsv \\
        --fdr_level ${params.pyprophet_global_fdr_level} \\
        --output "${sample}_${condition}.csv"
    """
    }
}


/*
 * STEP 14.5 - export_mztab
 */
process mztab_export {
    publishDir "${params.outdir}/"

    input:
    set val(id), val(sample), val(condition), file(dialignr_file) from DIALignR_result_mztab
    file exp_design from input_exp_design_mztab.first()
    set val(id), val(sample_lib), file(lib_file) from input_lib_used_I_mztab.first()

    output:
    set val(id), val(sample), val(condition), file("${sample}_${condition}.mzTab") into mztab_file

    when:
    params.mztab_export

    script:

    """
    TargetedFileConverter -in ${lib_file} -out ${lib_file.baseName}.tsv

    mztab_output.py \\
        --input ${dialignr_file} \\
        --exp_design ${exp_design} \\
        --library ${lib_file.baseName}.tsv \\
        --fdr_level ${params.pyprophet_global_fdr_level} \\
        --fdr_threshold_pep ${params.pyprophet_peptide_fdr} \\
        --fdr_threshold_prot ${params.pyprophet_protein_fdr} \\
        --ms1_scoring ${params.use_ms1} \\
        --rt_extraction_window ${params.rt_extraction_window} \\
        --mz_extraction_window ${params.mz_extraction_window} \\
        --mz_extraction_window_ms1 ${params.mz_extraction_window_ms1} \\
        --mz_extraction_unit ${params.mz_extraction_window_unit} \\
        --mz_extraction_unit_ms1 ${params.mz_extraction_window_ms1_unit} \\
        --dialignr_global_align_fdr ${params.dialignr_global_align_fdr} \\
        --dialignr_analyte_fdr ${params.dialignr_analyte_fdr} \\
        --dialignr_unalign_fdr ${params.dialignr_unalign_fdr} \\
        --dialignr_align_fdr ${params.dialignr_align_fdr} \\
        --dialignr_query_fdr ${params.dialignr_query_fdr} \\
        --workflow_version $workflow.manifest.version \\
        --output "${sample}_${condition}.mzTab"
    """
}


/*
 * STEP 15 - Run MSstats
 */
process statistical_post_processing {
    publishDir "${params.outdir}/"

    label 'process_low'

    input:
    set val(id), val(sample), val(condition), file(csv) from msstats_file.groupTuple(by:1)

    output:
    file "*.pdf" optional true // Output plots: 1) Comparative plots across pairwise conditions, 2) VolcanoPlot
    file "*.csv" // Csv of normalized differential protein abundancies calculated by msstats
    file "*.log" // logfile of msstats run

    when:
    params.run_msstats

    script:
    """
    msstats.R > msstats.log || echo "Optional MSstats step failed. Please check logs and re-run or do a manual statistical analysis."
    """
}


/*
 * STEP 16 - Generate plots describing output:
 * 1) BarChartProtein/Peptide Counts
 * 2) Pie Chart: Peptide Charge distribution
 * 3) Density Scatter: Library vs run RT deviations for all identifications
 * 4) Heatmap: Peptide quantities across MS runs
 * 5) Pyprophet score plots
 */
process output_visualization {
    publishDir "${params.outdir}/"

    label 'process_low'

    input:
    set val(sample), val(id), val(condition), file(quantity_csv_file), val(dummy_id), val(dummy_condition), file(pyprophet_tsv_file) from DIALignR_result_I.transpose().join(pyprophet_results, by:1)

    output:
    file "*.pdf" into output_plots

    when:
    params.generate_plots

    script:
    """
    plot_quantities_and_counts.R ${sample}
    """
}

/*
 * Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

    output:
    file 'results_description.html'

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/diaproteomics] Successful: $workflow.runName"
    if (!workflow.success) {
        subject = "[nf-core/diaproteomics] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir" ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/diaproteomics] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            mail_cmd.execute() << email_html
            log.info "[nf-core/diaproteomics] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/diaproteomics]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/diaproteomics]${c_red} Pipeline completed with errors${c_reset}-"
    }

}

workflow.onError {
    // Print unexpected parameters - easiest is to just rerun validation
    NfcoreSchema.validateParameters(params, json_schema, log)
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = 'hostname'.execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error '====================================================\n' +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            '============================================================'
                }
            }
        }
    }
}
