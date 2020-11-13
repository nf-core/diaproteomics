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

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/diaproteomics --input 'sample_sheet.tsv' --input_spectral_library 'library_sheet.tsv' --irts 'irt_sheet.tsv' -profile standard,docker

    Mandatory arguments:
      --input                           Path to input DIA raw/mzML/mzXML data sheet (must be surrounded with quotes)
      -profile                          Configuration profile to use. Can use multiple (comma separated)
                                        Available: standard, conda, docker, singularity, awsbatch, test
    DIA Mass Spectrometry Search:
      --input_spectral_library          Path to spectral library input sheet
      --irts                            Path to internal retention time standards input_sheet
      --irt_min_rsq			Minimal rsq error for irt RT alignment (default=0.95)
      --irt_n_bins                      Number of RT bins for iRT alignment
      --irt_min_bins_covered            Minimal number of RT bins covered for iRT alignment
      --irt_alignment_method            Method for irt RT alignment ('linear','lowess')
      --generate_spectral_library       Set flag if spectral libraries should be generated from provided DDA data (pepXML and mzML)
      --merge_libraries                 Set flag if multiple input spectral libraries should be merged by SampleID Column
      --align_libraries                 Set flag if multiple input spectral libraries should be aligned to the same RT reference
      --min_overlap_for_merging         Minimal number of peptides overlapping between libraries for RT alignment when merging.
      --generate_pseudo_irts            Set flag if pseudo irts should be generated from provided DDA data (pepXML and mzML)
      --n_irts                          Number of pseudo irts to be selected from dda data (default 250)
      --input_sheet_dda                 Path to input sheet of mzML DDA MS raw data and mzid, idXML or other formats of DDA search results to use for spectral library generation
      --library_rt_fdr                  PSM fdr threshold to align peptide ids with reference run (default = 0.01)
      --unimod                          Path to unimod.xml file describing modifications (https://github.com/nf-core/test-datasets/tree/diaproteomics)
      --skip_decoy_generation           Use a spectral library that already includes decoy sequences
      --decoy_method                    Method for generating decoys ('shuffle','pseudo-reverse','reverse','shift')
      --min_transitions                 Minimum number of transitions for assay
      --max_transitions                 Maximum number of transitions for assay
      --mz_extraction_window            Mass tolerance for transition extraction (ppm)
      --mz_extraction_window_ms1        Mass tolerance for precursor transition extraction (ppm)
      --rt_extraction_window            RT window for transition extraction (seconds)
      --pyprophet_classifier            Classifier used for target / decoy separation ('LDA','XGBoost')
      --pyprophet_fdr_ms_level          MS Level of FDR calculation ('ms1', 'ms2', 'ms1ms2')
      --pyprophet_global_fdr_level      Level of FDR calculation ('peptide', 'protein')
      --pyprophet_peakgroup_fdr         Threshold for FDR filtering
      --pyprophet_peptide_fdr           Threshold for global Peptide FDR
      --pyprophet_protein_fdr           Threshold for global Protein FDR
      --pyprophet_pi0_start             Start for non-parametric pi0 estimation
      --pyprophet_pi0_end               End for non-parametric pi0 estimation
      --pyprophet_pi0_steps             Steps for non-parametric pi0 estimation
      --DIAlignR_global_align_FDR       DIAlignR global Aligment FDR threshold
      --DIAlignR_analyte_FDR            DIAlignR Analyte FDR threshold
      --DIAlignR_unalign_FDR            DIAlignR UnAligment FDR threshold
      --DIAlignR_align_FDR              DIAlignR Aligment FDR threshold
      --DIAlignR_query_FDR              DIAlignR Query FDR threshold
      --run_msstats                     Set flag if MSstats should be run
      --generate_plots                  Set flag if plots should be generated and included in the output
      --force_option                    Force the analysis despite severe warnings

    Other options:
      --outdir [file]                 The output directory where the results will be saved
      --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move (Default: copy)
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      --max_multiqc_email_size [str]  Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

/*
 * SET UP CONFIGURATION VARIABLES
 */


// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$baseDir/docs/images/", checkIfExists: true)
sample_sheet = file(params.input)
Channel
 .from( sample_sheet )
 .set { input_exp_design}

params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()

// DIA MS input
Channel.from( sample_sheet )
       .splitCsv(header: true, sep:'\t')
       .map { col -> tuple("${col.Sample}", "${col.BatchID}", "${col.MSstats_Condition}", file("${col.Spectra_Filepath}", checkifExists: true))}
       .flatMap{it -> [tuple(it[0],it[1].toString(),it[2],it[3])]}
       .set {input_branch}


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


// Validate inputs
sample_sheet = file(params.input)


if( params.align_libraries) {
    align_libraries = 'true'
} else {
    align_libraries = 'false'
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

    input_dda.branch {
        raw: hasExtension(it[2], 'raw')
        mzml: hasExtension(it[2], 'mzML')
        mzxml: hasExtension(it[2], 'mzXML')
        other: true 
    }.set{input_dda_ms_files}

    Channel
        .fromPath( params.unimod )
        .ifEmpty { exit 1, "params.unimod was empty - no unimod.xml supplied" }
        .set { input_unimod}

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
        .flatMap{it -> [tuple(it[0],it[1],it[2])]}
        .set {input_lib_nd}

    input_lib = Channel.empty()
    input_lib_1 = Channel.empty()
    input_dda = Channel.empty()
    input_dda.branch {
        raw: hasExtension(it[2], 'raw')
        mzml: hasExtension(it[2], 'mzML')
        mzxml: hasExtension(it[2], 'mzXML')
        other: true
    }.set{input_dda_ms_files}
    input_unimod = Channel.empty()

} else {

    // Spectral library input
    library_sheet = file(params.input_spectral_library)

    Channel.from( library_sheet )
        .splitCsv(header: true, sep:'\t')
        .map { col -> tuple("${col.Sample}", "${col.BatchID}", file("${col.Library_Filepath}", checkifExists: true))}
        .flatMap{it -> [tuple(it[0],it[1],it[2])]}
        .into {input_lib; input_lib_1 }

    input_lib_nd = Channel.empty()
    input_lib_nd_1 = Channel.empty()
    input_dda = Channel.empty()
    input_dda.branch {
        raw: hasExtension(it[2], 'raw')
        mzml: hasExtension(it[2], 'mzML')
        mzxml: hasExtension(it[2], 'mzXML')
        other: true
    }.set{input_dda_ms_files}
    input_unimod = Channel.empty()
}


if( !params.generate_pseudo_irts){
   // iRT library input
   irt_sheet = file(params.irts)

   Channel.from( irt_sheet )
       .splitCsv(header: true, sep:'\t')
       .map { col -> tuple("${col.BatchID}", file("${col.irt_Filepath}", checkifExists: true))}
       .flatMap{it -> [tuple(it[0],it[1])]}
       .set {input_irts}

} else {

    input_irts = Channel.empty()

}


// Force option
if (params.force_option){
    force_option='-force'
   } else {
    force_option=''
}


// Header log info
log.info nfcoreHeader()
def summary = [:]
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Run Name']         = custom_runName ?: workflow.runName
summary['Spectral Library']    = params.input_spectral_library
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
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

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
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}


/*
 * STEP 0 - Raw File Conversion
 */
process dda_raw_file_conversion {
    input:
     set val(id), val(Sample), file(raw_file), file(dda_id_file) from input_dda_ms_files.raw

    output:
     set val(id), val(Sample), file("${raw_file.baseName}.mzML"), file(dda_id_file) into converted_dda_input_mzmls

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
     set val(id), val(Sample), file(dda_mzml), file(dda_id_file) from input_dda_ms_files.mzml.mix(input_dda_ms_files.mzxml).mix(converted_dda_input_mzmls)

    output:
     set val(id), val(Sample), file(dda_mzml), file("${id}_${Sample}_peptide_ids.idXML") into input_dda_converted

    when:
     params.generate_spectral_library

    script:
     """
     IDFileConverter -in ${dda_id_file} \\
                     -out ${id}_${Sample}_peptide_ids.idXML
     """
}


/*
 * STEP 2 - Spectral Library Generation using EasyPQP
 */
process dda_library_generation {

    input:
     set val(id), val(Sample), file(dda_mzml_file), file(idxml_file) from input_dda_converted
     file unimod_file from input_unimod.first()

    output:
     set val(id), val(Sample), file("${id}_${Sample}_library.tsv") into input_lib_dda_nd

    when:
     params.generate_spectral_library

    script:
     """
     easypqp convert --unimod ${unimod_file} \\
                     --pepxml ${idxml_file} \\
                     --spectra ${dda_mzml_file} \\

     easypqp library --out ${dda_mzml_file.baseName}_run_peaks.tsv \\
                     --rt_psm_fdr_threshold ${params.library_rt_fdr} \\
                     --nofdr \\
                     ${dda_mzml_file.baseName}.psmpkl \\
                     ${dda_mzml_file.baseName}.peakpkl \\

     mv ${dda_mzml_file.baseName}_run_peaks.tsv ${id}_${Sample}_library.tsv
     """
}


/*
 * STEP 3 - Assay Generation for Spectral Library
 */
process assay_generation {

    input:
     set val(id), val(Sample), file(lib_file_na) from input_lib.mix(input_lib_dda_nd)

    output:
     set val(id), val(Sample), file("${id}_${Sample}_assay.tsv") into (input_lib_assay, input_lib_assay_for_irt, input_lib_assay_for_merging)

    when:
     !params.skip_decoy_generation

    script:
     """
     TargetedFileConverter -in ${lib_file_na} \\
                           -out ${lib_file_na.baseName}.tsv

     OpenSwathAssayGenerator -in ${lib_file_na.baseName}.tsv \\
                             -min_transitions ${params.min_transitions} \\
                             -max_transitions ${params.max_transitions} \\
                             -out ${id}_${Sample}_assay.tsv \\
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
     set val(id), val(Sample), file(lib_files_for_merging) from input_lib_assay_for_merging.groupTuple(by:1)

    output:
     set val(id), val(Sample), file("${Sample}_library_merged.tsv") into (input_lib_assay_merged, input_lib_assay_merged_for_irt)
     set val(id), val(Sample), file("*.png") optional true 

    when:
     params.merge_libraries

    script:
     """
     merge_and_align_libraries_from_easypqp.py --input_libraries ${lib_files_for_merging} --min_overlap ${params.min_overlap_for_merging} --rsq_threshold 0.75 --align ${align_libraries} --output ${Sample}_library_merged.tsv
     """
}


/*
 * STEP 5 - Pseudo iRT Library Generation
 */
process pseudo_irt_generation {
    publishDir "${params.outdir}/spectral_library_files"

    input:
     set val(id), val(Sample), file(lib_file_assay_irt) from input_lib_assay_for_irt.mix(input_lib_assay_merged_for_irt)

    output:
     set val(Sample), file("${lib_file_assay_irt.baseName}_pseudo_irts.pqp") into input_lib_assay_irt_2

    when:
     params.generate_pseudo_irts

    script:
     """
     select_pseudo_irts_from_lib.py --input_libraries ${lib_file_assay_irt} --min_rt 0 --n_irts ${params.n_irts} --max_rt 100 --output ${lib_file_assay_irt.baseName}_pseudo_irts.tsv \\

     TargetedFileConverter -in ${lib_file_assay_irt.baseName}_pseudo_irts.tsv \\
                           -out ${lib_file_assay_irt.baseName}_pseudo_irts.pqp
     """
}


/*
 * STEP 6 - Decoy Generation for Spectral Library
 */
process decoy_generation {
    publishDir "${params.outdir}/spectral_library_files"

    input:
     set val(id), val(Sample), file(lib_file_nd) from input_lib_assay.mix(input_lib_assay_merged)

    output:
     set val(id), val(Sample), file("${lib_file_nd.baseName}_decoy.pqp") into input_lib_decoy

    when:
     !params.skip_decoy_generation

    script:
     """
     TargetedFileConverter -in ${lib_file_nd} \\
                           -out ${lib_file_nd.baseName}.pqp \\

     OpenSwathDecoyGenerator -in ${lib_file_nd.baseName}.pqp \\
                             -method ${params.decoy_method} \\
                             -out ${lib_file_nd.baseName}_decoy.pqp \\
     """
}


/*
 * STEP 7 - DIA Raw File Conversion
 */
process dia_raw_file_conversion {

    input:
     set val(id), val(Sample), val(Condition), file(raw_file) from input_dia_ms_files.raw

    output:
     set val(id), val(Sample), val(Condition), file("${raw_file.baseName}.mzML") into converted_dia_input_mzmls

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

    input:
     set val(Sample), val(id), val(Condition), file(mzml_file), val(dummy_id), file(lib_file), file(irt_file) from converted_dia_input_mzmls.mix(input_dia_ms_files.mzml.mix(input_dia_ms_files.mzxml)).combine(input_lib_decoy.mix(input_lib_nd), by:1).combine(input_irts.mix(input_lib_assay_irt_2), by:0)

    output:
     set val(id), val(Sample), val(Condition), file("${mzml_file.baseName}_chrom.mzML") into chromatogram_files
     set val(id), val(Sample), val(Condition), file("${mzml_file.baseName}.osw") into osw_files
     set val(id), val(Sample), file(lib_file) into (input_lib_used, input_lib_used_I)

    script:
     """
     TargetedFileConverter -in ${lib_file} \\
                           -out ${lib_file.baseName}.pqp \\

     TargetedFileConverter -in ${irt_file} \\
                           -out ${irt_file.baseName}.pqp \\

     OpenSwathWorkflow -in ${mzml_file} \\
                       -tr ${lib_file.baseName}.pqp \\
                       -sort_swath_maps \\
                       -tr_irt ${irt_file.baseName}.pqp \\
                       -min_rsq ${params.irt_min_rsq} \\
                       -out_osw ${mzml_file.baseName}.osw \\
                       -out_chrom ${mzml_file.baseName}_chrom.mzML \\
                       -mz_extraction_window ${params.mz_extraction_window} \\
                       -mz_extraction_window_ms1 ${params.mz_extraction_window_ms1} \\
                       -mz_extraction_window_unit 'ppm' \\
                       -mz_extraction_window_ms1_unit 'ppm' \\
                       -rt_extraction_window ${params.rt_extraction_window} \\
                       -RTNormalization:alignmentMethod ${params.irt_alignment_method} \\
                       -RTNormalization:estimateBestPeptides \\
                       -RTNormalization:outlierMethod none \\
                       -RTNormalization:NrRTBins ${params.irt_n_bins} \\
                       -RTNormalization:MinBinsFilled ${params.irt_min_bins_covered} \\
                       -mz_correction_function quadratic_regression_delta_ppm \\
                       -use_ms1_traces \\
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
                       -Scoring:Scores:use_ms1_mi \\
                       -Scoring:Scores:use_mi_score \\
                       -batchSize 1000 \\
                       -Scoring:DIAScoring:dia_nr_isotopes 3 \\
                       -enable_uis_scoring \\
                       -Scoring:uis_threshold_sn -1 \\
                       -threads ${task.cpus} \\
                       ${force_option} \\  
     """
}


/*
 * STEP 9 - Pyprophet merging of OpenSwath results
 */
process dia_search_output_merging {

    input:
     set val(Sample), val(id), val(Condition), file(all_osws), val(dummy_id), file(lib_file_template) from osw_files.groupTuple(by:1).join(input_lib_used, by:1)

    output:
     set val(id), val(Sample), val(Condition), file("${Sample}_osw_file_merged.osw") into (merged_osw_file, merged_osw_file_for_global)

    script:
     """
     pyprophet merge --template=${lib_file_template} \\
                     --out=${Sample}_osw_file_merged.osw \\
                     --no-same_run \\
                     ${all_osws} \\
     """
}


/*
 * STEP 10 - Pyprophet FDR Scoring
 */
process false_discovery_rate_estimation {
    publishDir "${params.outdir}/pyprophet_output"

    input:
     set val(id), val(Sample), val(Condition), file(merged_osw) from merged_osw_file

    output:
     set val(id), val(Sample), val(Condition), file("${merged_osw.baseName}_scored_merged.osw") into merged_osw_scored_for_pyprophet
     set val(id), val(Sample), val(Condition), file("*.pdf") into target_decoy_score_plots

    when:
     params.pyprophet_global_fdr_level==''

    script:
    if (params.pyprophet_classifier=='LDA'){
     """
     pyprophet score --in=${merged_osw} \\
                     --level=${params.pyprophet_fdr_ms_level} \\
                     --out=${merged_osw.baseName}_scored_merged.osw \\
                     --classifier=${params.pyprophet_classifier} \\
                     --pi0_lambda ${params.pyprophet_pi0_start} ${params.pyprophet_pi0_end} ${params.pyprophet_pi0_steps} \\
                     --threads=${task.cpus} \\
     """
    } else {
     """
     pyprophet score --in=${merged_osw} \\
                     --level=${params.pyprophet_fdr_ms_level} \\
                     --out=${merged_osw.baseName}_scored_merged.osw \\
                     --classifier=${params.pyprophet_classifier} \\
                     --threads=${task.cpus} \\
     """
    }
}


/*
 * STEP 11 - Pyprophet global FDR Scoring
 */
process global_false_discovery_rate_estimation {
    publishDir "${params.outdir}/pyprophet_output"

    input:
     set val(id), val(Sample), val(Condition), file(scored_osw) from merged_osw_file_for_global

    output:
     set val(id), val(Sample), val(Condition), file("${scored_osw.baseName}_global_merged.osw") into merged_osw_scored_global_for_pyprophet
     set val(id), val(Sample), val(Condition), file("*.pdf") into target_decoy_global_score_plots

    when:
     params.pyprophet_global_fdr_level!=''

    script:
    if (params.pyprophet_classifier=='LDA'){
     """
     pyprophet score --in=${scored_osw} \\
                     --level=${params.pyprophet_fdr_ms_level} \\
                     --out=${scored_osw.baseName}_scored.osw \\
                     --classifier=${params.pyprophet_classifier} \\
                     --pi0_lambda ${params.pyprophet_pi0_start} ${params.pyprophet_pi0_end} ${params.pyprophet_pi0_steps} \\
                     --threads=${task.cpus} \\

     pyprophet ${params.pyprophet_global_fdr_level} --in=${scored_osw.baseName}_scored.osw \\
                                                    --out=${scored_osw.baseName}_global_merged.osw \\
                                                    --context=global \\
     """
    } else {
     """
     pyprophet score --in=${scored_osw} \\
                     --level=${params.pyprophet_fdr_ms_level} \\
                     --out=${scored_osw.baseName}_scored.osw \\
                     --classifier=${params.pyprophet_classifier} \\
                     --threads=${task.cpus} \\

     pyprophet ${params.pyprophet_global_fdr_level} --in=${scored_osw.baseName}_scored.osw \\
                                                    --out=${scored_osw.baseName}_global_merged.osw \\
                                                    --context=global \\
     """
    }
}


/*
 * STEP 12 - Pyprophet Export
 */
process export_of_scoring_results {
    publishDir "${params.outdir}/pyprophet_output"

    input:
     set val(id), val(Sample), val(Condition), file(global_osw) from merged_osw_scored_for_pyprophet.mix(merged_osw_scored_global_for_pyprophet)

    output:
     set val(id), val(Sample), val(Condition), file("*.tsv") into pyprophet_results
     set val(id), val(Sample), val(Condition), file(global_osw) into osw_for_dialignr

    script:
     """
     pyprophet export --in=${global_osw} \\
                      --max_rs_peakgroup_qvalue=${params.pyprophet_peakgroup_fdr} \\
                      --max_global_peptide_qvalue=${params.pyprophet_peptide_fdr} \\
                      --max_global_protein_qvalue=${params.pyprophet_protein_fdr} \\
                      --out=legacy.tsv \\
     """
}


/*
 * STEP 13 - Index Chromatogram mzMLs
 */
process chromatogram_indexing {

    input:
     set val(id), val(Sample), val(Condition), file(chrom_file_noindex) from chromatogram_files

    output:
     set val(id), val(Sample), val(Condition), file("${chrom_file_noindex.baseName.split('_chrom')[0]}.chrom.mzML") into chromatogram_files_indexed

    script:
     """
     FileConverter -in ${chrom_file_noindex} \\
                   -out ${chrom_file_noindex.baseName.split('_chrom')[0]}.chrom.mzML \\
     """
}


// Combine channels of osw files and osw chromatograms
osw_for_dialignr
 .transpose()
 .join(chromatogram_files_indexed, by:1)
 .groupTuple(by:0)
 // Channel contains now the following elements:
 // ([id, [Samples], [Conditions], [osw_files], id_2, condition_2, [chromatogram_files]])
 .flatMap{it -> [tuple(it[0],it[1].unique()[0],it[2].unique()[0],it[3].unique()[0],it[4],it[5],it[6])]}
 .set{osw_and_chromatograms_combined_by_condition}

/*
 * STEP 10 - Align DIA Chromatograms using DIAlignR
 */
process chromatogram_alignment {
    publishDir "${params.outdir}/"

    input:
     set val(Sample), val(id), val(Condition), file(pyresults), val(id_dummy), val(condition_dummy), file(chrom_files_index) from osw_and_chromatograms_combined_by_condition

    output:
     set val(id), val(Sample), val(Condition), file("${Sample}_peptide_quantities.csv") into (DIALignR_result, DIALignR_result_I)

    script:
     """
     mkdir osw
     mv ${pyresults} osw/ 
     mkdir mzml 
     mv *.chrom.mzML mzml/

     DIAlignR.R ${params.DIAlignR_global_align_FDR} ${params.DIAlignR_analyte_FDR} ${params.DIAlignR_unalign_FDR} ${params.DIAlignR_align_FDR} ${params.DIAlignR_query_FDR}

     mv DIAlignR.csv ${Sample}_peptide_quantities.csv
     """
}


/*
 * STEP 11 - Reformat output for MSstats: Combine with experimental design and missing columns from input library
 */
process reformatting {

   input:
    set val(id), val(Sample), val(Condition), file(dialignr_file) from DIALignR_result
    file exp_design from input_exp_design.first()
    set val(id), val(Sample), file(lib_file) from input_lib_used_I.first()

   output:
    set val(id), val(Sample), val(Condition), file("${Sample}_${Condition}.csv") into msstats_file

   when:
    params.run_msstats

   script:

    if (params.pyprophet_global_fdr_level==''){

    """
     TargetedFileConverter -in ${lib_file} \\
                           -out ${lib_file.baseName}.tsv

     reformat_output_for_msstats.py --input ${dialignr_file} --exp_design ${exp_design} --library ${lib_file.baseName}.tsv --fdr_level "none" --output "${Sample}_${Condition}.csv"
    """

    } else {

    """
     TargetedFileConverter -in ${lib_file} \\
                           -out ${lib_file.baseName}.tsv

     reformat_output_for_msstats.py --input ${dialignr_file} --exp_design ${exp_design} --library ${lib_file.baseName}.tsv --fdr_level ${params.pyprophet_global_fdr_level} --output "${Sample}_${Condition}.csv"
    """
    }
}


/*
 * STEP 12 - Run MSstats
 */
process statistical_post_processing {
   publishDir "${params.outdir}/"

   input:
    set val(id), val(Sample), val(Condition), file(csv) from msstats_file.groupTuple(by:1)

   output:
    file "*.pdf" // Output plots: 1) Comparative plots across pairwise conditions, 2) VolcanoPlot
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
 * STEP 13 - Generate plots describing output:
 * 1) BarChartProtein/Peptide Counts
 * 2) Pie Chart: Peptide Charge distribution
 * 3) Density Scatter: Library vs run RT deviations for all identifications
 * 4) Heatmap: Peptide quantities across MS runs
 * 5) Pyprophet score plots
 */
process output_visualization {
   publishDir "${params.outdir}/"

   input:
    set val(Sample), val(id), val(Condition), file(quantity_csv_file), val(dummy_id), val(dummy_Condition), file(pyprophet_tsv_file) from DIALignR_result_I.transpose().join(pyprophet_results, by:1)

   output:
    file "*.pdf" into output_plots

   when:
    params.generate_plots

   script:
    """
    plot_quantities_and_counts.R ${Sample}
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
    file "results_description.html"

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
    email_fields['runName'] = custom_runName ?: workflow.runName
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

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/diaproteomics] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/diaproteomics] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
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
            if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
              mail_cmd += [ '-A', mqc_report ]
            }
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


def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/diaproteomics v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
