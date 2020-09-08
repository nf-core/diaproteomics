# nf-core/diaproteomics: Usage

## Table of contents

* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Running the pipeline](#running-the-pipeline)
  * [Updating the pipeline](#updating-the-pipeline)
  * [Reproducibility](#reproducibility)
* [Main arguments](#main-arguments)
  * [`--input`](#--input)
  * [`--swath_windows`](#--swath_windows)
  * [`--spectral_lib`](#--spectral_lib)
  * [`--irts`](#--irts)
  * [`--irt_min_rsq`](#--irt_min_rsq)
  * [`--irt_alignment_method`](#--irt_alignment_method)
  * [`--generate_spectral_lib`](#--generate_spectral_lib)
  * [`--dda_id`](#--dda_id)
  * [`--dda_mzml`](#--dda_mzml)
  * [`--library_rt_fdr`](#--library_rt_fdr)
  * [`--unimod`](#--unimod)
  * [`--skip_decoy_generation`](#--skip_decoy_generation)
  * [`--decoy_method`](#--decoy_method)
  * [`--min_transitions`](#--min_transitions)
  * [`--max_transitions`](#--max_transitions)
  * [`--mz_extraction_window`](#--mz_extraction_window)
  * [`--rt_extraction_window`](#--rt_extraction_window)
  * [`--pyprophet_classifier`](#--pyprophet_classifier)
  * [`--pyprophet_fdr_ms_level`](#--pyprophet_fdr_ms_level)
  * [`--pyprophet_global_fdr_level`](#--pyprophet_global_fdr_level)
  * [`--pyprophet_peakgroup_fdr`](#--pyprophet_peakgroup_fdr)
  * [`--pyprophet_peptide_fdr`](#--pyprophet_peptide_fdr)
  * [`--pyprophet_protein_fdr`](#--pyprophet_protein_fdr)
  * [`--pyprophet_pi0_start`](#--pyprophet_pi0_start)
  * [`--pyprophet_pi0_end`](#--pyprophet_pi0_end)
  * [`--pyprophet_pi0_steps`](#--pyprophet_pi0_steps)
  * [`--DIAlignR_global_align_FDR`](#--DIAlignR_global_align_FDR)
  * [`--DIAlignR_analyte_FDR`](#--DIAlignR_analyte_FDR)
  * [`--prec_charge`](#--prec_charge)
  * [`--force_option`](#--force_option)
  * [`-profile`](#-profile)
* [Job resources](#job-resources)
  * [Automatic resubmission](#automatic-resubmission)
  * [Custom resource requests](#custom-resource-requests)
* [AWS Batch specific parameters](#aws-batch-specific-parameters)
  * [`--awsqueue`](#--awsqueue)
  * [`--awsregion`](#--awsregion)
  * [`--awscli`](#--awscli)
* [Other command line parameters](#other-command-line-parameters)
  * [`--outdir`](#--outdir)
  * [`--publish_dir_mode`](#--publish_dir_mode)
  * [`--email`](#--email)
  * [`--email_on_fail`](#--email_on_fail)
  * [`--max_multiqc_email_size`](#--max_multiqc_email_size)
  * [`-name`](#-name)
  * [`-resume`](#-resume)
  * [`-c`](#-c)
  * [`--custom_config_version`](#--custom_config_version)
  * [`--custom_config_base`](#--custom_config_base)
  * [`--max_memory`](#--max_memory)
  * [`--max_time`](#--max_time)
  * [`--max_cpus`](#--max_cpus)
  * [`--plaintext_email`](#--plaintext_email)
  * [`--monochrome_logs`](#--monochrome_logs)
  * [`--multiqc_config`](#--multiqc_config)

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/diaproteomics --input '*_R{1,2}.fastq.gz' -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/diaproteomics
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/diaproteomics releases page](https://github.com/nf-core/diaproteomics/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Main arguments

### `--input`

Use this to specify the location of your input dia raw files (mzML). For example:

```bash
--input 'path/to/data/sample_*.mzML'
```

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character

### `--swath_windows`

Path to swath_windows.txt file (tsv), containing swath window mz ranges

```bash
--swath_windows 'path/to/data/swath_windows.txt'
```

### `--spectral_lib`

Path to spectral library input file (pqp, TraML or tsv). The library should not contain decoys and can be unfiltered, as these are later steps within the pipeline.

```bash
--spectral_lib 'path/to/data/spectral_library.pqp'
```

### `--irts`

Path to internal retention time standards (pqp, TraML or tsv)

```bash
--irts 'path/to/data/irts.TraML'
```

### `--irt_min_rsq`

Minimal rsq error for irt RT alignment (default=0.95)

### `--irt_alignment_method`

Method for irt RT alignment for example 'linear' or 'lowess'.

### `--generate_spectral_lib`

Set this flag if the spectral library should be generated using EasyPQP from provided DDA data - identification search results and corresponding raw data.

### `--dda_id`

Path to mzid, idXML or other formats of DDA search results to use for spectral library generation

```bash
--dda_id 'path/to/data/dda_peptide_identifications.pepXML'
```

### `--dda_mzml`

Path to corresponding dda raw data to generate spectral library (mzML)

```bash
--dda_mzml 'path/to/data/dda_peptide_identifications.mzML'
```

### `--library_rt_fdr`

PSM fdr threshold to align peptide ids with reference run (default = 0.01)

### `--unimod`

Path to unimod.xml file describing modifications ("https://github.com/nf-core/test-datasets/tree/diaproteomics/unimod.xml")

### `--skip_decoy_generation`

Set this flag if using a spectral library that already includes decoy sequences and therefor skip assay and decoy generation.

### `--decoy_method`

Method for generating decoys ('shuffle','pseudo-reverse','reverse','shift')

### `--min_transitions`

Minimum number of transitions for assay

### `--max_transitions`

Maximum number of transitions for assay

### `--mz_extraction_window`

Mass tolerance for transition extraction (ppm)

### `--rt_extraction_window`

RT window for transition extraction (seconds)

### `--pyprophet_classifier`

Machine learning lassifier used for pyprophet target / decoy separation ('LDA','XGBoost')

### `--pyprophet_fdr_ms_level`

MS Level of pyprophet FDR calculation: 'ms1', 'ms2' or both 'ms1ms2'

### `--pyprophet_global_fdr_level`

Abstraction level of pyrophet FDR calculation ('peptide', 'protein')

### `--pyprophet_peakgroup_fdr`

Threshold for pyprophet FDR filtering on peakgroup abstraction level

### `--pyprophet_peptide_fdr`

Threshold for pyprophet FDR filtering on peptide abstraction level

### `--pyprophet_protein_fdr`

Threshold for pyprophet FDR filtering on protein abstraction level

### `--pyprophet_pi0_start`

Start for pyprophet non-parametric pi0 estimation

### `--pyprophet_pi0_end`

End for pyprophet non-parametric pi0 estimation

### `--pyprophet_pi0_steps`

Steps for pyprophet non-parametric pi0 estimation

### `--DIAlignR_global_align_FDR`

DIAlignR global Aligment FDR threshold

### `--DIAlignR_analyte_FDR`

DIAlignR Analyte FDR threshold

### `--DIAlignR_unalign_FDR`

DIAlignR UnAligment FDR threshold

### `--DIAlignR_align_FDR`

DIAlignR Aligment FDR threshold

### `--DIAlignR_query_FDR`

DIAlignR Query FDR threshold

### `--prec_charge`

Precursor charge eg. "2:3"

### `--force_option`

Force the analysis of the OpenSwathWorkflow despite severe warnings

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](https://docker.com/)
  * Pulls software from Docker Hub: [`nfcore/diaproteomics`](https://hub.docker.com/r/nfcore/diaproteomics/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
  * Pulls software from Docker Hub: [`nfcore/diaproteomics`](https://hub.docker.com/r/nfcore/diaproteomics/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use [`-profile awsbatch`](https://github.com/nf-core/configs/blob/master/conf/awsbatch.config) and then specify all of the following parameters.

### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

### `--awsregion`

The AWS region in which to run your job. Default is set to `eu-west-1` but can be adjusted to your needs.

### `--awscli`

The [AWS CLI](https://www.nextflow.io/docs/latest/awscloud.html#aws-cli-installation) path in your custom AMI. Default: `/home/ec2-user/miniconda/bin/aws`.

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### `--outdir`

The output directory where the results will be saved.

### `--publish_dir_mode`

Value passed to Nextflow [`publishDir`](https://www.nextflow.io/docs/latest/process.html#publishdir) directive for publishing results in the output directory. Available: 'symlink', 'rellink', 'link', 'copy', 'copyNoFollow' and 'move' (Default: 'copy').

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `--max_multiqc_email_size`

Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB).

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

### Running in the background

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default: `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.
