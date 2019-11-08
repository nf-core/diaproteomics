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

    nextflow run nf-core/diaproteomics --dia_mzmls '*.mzML' --spectral_lib '*.pqp' --irts '*.pqp' --swath_windows '*.txt' -profile standard,docker

    Mandatory arguments:
      --dia_mzmls                       Path to input data (must be surrounded with quotes)
      --swath_windows                   Path to swath_windows.txt file, containing swath window mz ranges
      -profile                          Configuration profile to use. Can use multiple (comma separated)
                                        Available: standard, conda, docker, singularity, awsbatch, test

    DIA Mass Spectrometry Search:
      --spectral_lib                    Path to spectral library input file (pqp)
      --irts                            Path to internal retention time standards (pqp)
      --irt_min_rsq			Minimal rsq error for irt RT alignment (default=0.95)
      --generate_spectral_lib           Set flag if spectral lib should be generated from provided DDA data (pepXML and mzML)
      --dda_pepxmls                     Path to DDA pepXML input for library generation
      --dda_mzmls                       Path to DDA mzML input for library generation
      --skip_decoy_generation           Use a spectral library that already includes decoy sequences
      --decoy_method                    Method for generating decoys ('shuffle','pseudo-reverse','reverse','shift')
      --min_transitions                 Minimum peptide length for filtering
      --max_transitions                 Maximum peptide length for filtering
      --mz_extraction_window            Mass tolerance for transition extraction (ppm)
      --rt_extraction_window            RT window for transition extraction (seconds)
      --fdr_threshold                   Threshold for FDR filtering
      --fdr_level                       Level of FDR calculation ('ms1', 'ms2', 'transition')
      --prec_charge                     Precursor charge (eg. "2:3")
      --force_option                    Force the Analysis despite severe warnings

    Other options:
      --outdir                          The output directory where the results will be saved
      --email                           Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                             Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                        The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                       The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}


// Validate inputs
params.dia_mzmls = params.dia_mzmls ?: { log.error "No dia mzml data provided. Make sure you have used the '--dia_mzmls' option."; exit 1 }()
params.swath_windows = params.swath_windows ?: { log.error "No swath windows provided. Make sure you have used the '--swath_windows' option."; exit 1 }()
params.irts = params.irts ?: { log.error "No internal retention time standards provided. Make sure you have used the '--irts' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()


/*
 * Define the default parameters
 */

//MS params
params.generate_spectral_lib = false
params.skip_decoy_generation = false

params.irt_min_rsq = 0.95
params.min_transitions = 4
params.max_transitions = 6
params.mz_extraction_window = 30
params.rt_extraction_window = 600
params.fdr_threshold = 0.01
params.fdr_level = 'ms2'

params.number_mods = 3
params.num_hits = 1
params.prec_charge = '2:3'
params.variable_mods = 'Oxidation (M)'
params.decoy_method = 'shuffle'
params.spectrum_batch_size = 500

params.force_option = false
params.force_option ?: { force_option='-force'; force_option='' }()
/*
 * SET UP CONFIGURATION VARIABLES
 */


// Configurable variables
params.name = false
params.email = false
params.plaintext_email = false

output_docs = file("$baseDir/docs/output.md")


// AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}
//
// NOTE - THIS IS NOT USED IN THIS PIPELINE, EXAMPLE ONLY
// If you want to use the above in a process, define the following:
//   input:
//   file fasta from fasta
//

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


if( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

// Stage config files
ch_multiqc_config = Channel.fromPath(params.multiqc_config)
ch_output_docs = Channel.fromPath("$baseDir/docs/output.md")



Channel.fromPath( params.dia_mzmls )
        .ifEmpty { exit 1, "Cannot find any mzmls matching: ${params.dia_mzmls}\nNB: Path needs to be enclosed in quotes!" }
        .set { input_mzmls }

Channel.fromPath( params.swath_windows)
        .ifEmpty { exit 1, "Cannot find any swath_windows matching: ${params.swath_windows}\nNB: Path needs to be enclosed in quotes!" }
        .set { input_swath_windows }

Channel.fromPath( params.irts)
        .ifEmpty { exit 1, "Cannot find any irts matching: ${params.irts}\nNB: Path needs to be enclosed in quotes!" }
        .set { input_irts }

/*
 * Create a channel for input spectral library
 */
if( params.generate_spectral_lib) {

    input_spectral_lib = Channel.empty()

} else if( !params.skip_decoy_generation) {
    Channel
        .fromPath( params.spectral_lib )
        .ifEmpty { exit 1, "params.spectral_lib was empty - no input spectral library supplied" }
        .set { input_lib_nd }

    input_lib = Channel.empty()
    input_lib_1 = Channel.empty()

} else {
    Channel
        .fromPath( params.spectral_lib )
        .ifEmpty { exit 1, "params.spectral_lib was empty - no input spectral library supplied" }
        .into { input_lib; input_lib_1 }

    input_lib_decoy = Channel.empty()
    input_lib_decoy_1 = Channel.empty()

}



// Header log info
log.info nfcoreHeader()
def summary = [:]
summary['Pipeline Name']  = 'nf-core/diaproteomics'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['mzMLs']        = params.dia_mzmls
summary['Spectral Library']    = params.spectral_lib
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
summary['Config Profile'] = workflow.profile
if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if(params.email) {
  summary['E-mail Address']  = params.email
  summary['MultiQC maxsize'] = params.maxMultiqcEmailFileSize
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[2m----------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-diaproteomics-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/diaproteomics Workflow Summary'
    section_href: 'https://github.com/nf-core/diaproteomics'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/*
 * Parse software version numbers
 */
process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy',
    saveAs: {filename ->
        if (filename.indexOf(".csv") > 0) filename
        else null
    }
    
    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file "software_versions.csv"

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    FileInfo --help &> v_openms.txt
    pyprophet --version &> v_pyprophet.txt
    #python -c "import msproteomicstoolslib; print(msproteomicstoolslib.__version__)" &> v_msproteomicstools.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}


/*
 * STEP 0 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}


/*
 * STEP 0.5 - Decoy Generation for Spectral Library
 */
process generate_decoys_for_spectral_library {
    publishDir "${params.outdir}/"

    input:
     file lib_file_nd from input_lib_nd

    output:
     file "${lib_file_nd.baseName}_decoy.pqp" into (input_lib_decoy, input_lib_decoy_1)

    when:
     !params.skip_decoy_generation

    script:
     """
     OpenSwathDecoyGenerator -in ${lib_file_nd} \\
                             -method ${params.decoy_method} \\
                             -out "${lib_file_nd.baseName}_decoy.pqp" \\
     """
}


/*
 * STEP 1 - OpenSwathWorkFlow
 */
process run_openswathworkflow {
    publishDir "${params.outdir}/"

    input:
     file mzml_file from input_mzmls
     file swath_file from input_swath_windows.first()
     file lib_file from input_lib_decoy.mix(input_lib).first()
     file irt_file from input_irts.first()

    output:
     file "${mzml_file.baseName}_chrom.mzML" into chromatogram_files
     file "${mzml_file.baseName}.osw" into osw_files

    script:
     """
     OpenSwathWorkflow -in ${mzml_file} \\
                       -tr ${lib_file} \\
                       -swath_windows_file ${swath_file} \\
                       -tr_irt ${irt_file} \\
                       -min_rsq ${params.irt_min_rsq} \\
                       -out_osw ${mzml_file.baseName}.osw \\
                       -out_chrom ${mzml_file.baseName}_chrom.mzML \\
                       -mz_extraction_window ${params.mz_extraction_window} \\
                       -ppm \\
                       -rt_extraction_window ${params.rt_extraction_window} \\
                       -RTNormalization:alignmentMethod linear \\
                       -RTNormalization:outlierMethod none \\
                       -threads ${task.cpus} \\
                       ${force_option} \\
     """
}


/*
 * STEP 2 - Pyprophet merging of OpenSwath results
 */
process merge_openswath_output {
    publishDir "${params.outdir}/"

    input:
     file all_osws from osw_files.collect{it}
     file lib_file_1 from input_lib_decoy_1.mix(input_lib_1).first()

    output:
     file "merged_osw_file.osw" into merged_osw_file

    script:
     """
     pyprophet merge --template=${lib_file_1} \\
                     --out=merged_osw_file.osw \\
                     ${all_osws} \\
     """
}


/*
 * STEP 3 - Pyprophet FDR Scoring
 */
process run_fdr_scoring {
    publishDir "${params.outdir}/"

    input:
     file merged_osw from merged_osw_file

    output:
     file "${merged_osw.baseName}_scored.osw" into merged_osw_scored

    script:
     """
     pyprophet score --in=${merged_osw} \\
                     --level=${params.fdr_level} \\
                     --out=${merged_osw.baseName}_scored.osw \\
                     --threads=${task.cpus} \\
     """
}


/*
 * STEP 4 - Pyprophet Export
 */
process export_pyprophet_results {
    publishDir "${params.outdir}/"

    input:
     file scored_osw from merged_osw_scored

    output:
     file "*.tsv" into pyprophet_results

    script:
     """
     pyprophet export --in=${scored_osw} \\
                      --out=legacy.tsv \\
     """
}


/*
 * STEP 5 - Align DIA Chromatograms using TRIC
 */
process align_dia_runs {
    publishDir "${params.outdir}/"

    input:
     file pyresults from pyprophet_results

    output:
     file "aligned.tsv" into TRIC_result

    script:
     """
     feature_alignment.py --in ${pyresults}
                          --out aligned.tsv
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
    if (workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // TODO nf-core: If not using MultiQC, strip out this code (including params.maxMultiqcEmailFileSize)
    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = multiqc_report.getVal()
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
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
          if ( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/diaproteomics] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, email_address ].execute() << email_txt
          log.info "[nf-core/diaproteomics] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if (!output_d.exists()) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if (workflow.success) {
        log.info "${c_purple}[nf-core/diaproteomics]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[nf-core/diaproteomics]${c_red} Pipeline completed with errors${c_reset}"
    }

}


def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

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

def checkHostname(){
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
