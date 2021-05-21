# Call Somatic Structural Variant Pipeline

* [Overview](#Overview)
* [How To Run](#How-To-Run)
* [Flow Diagram](#flow-diagram)
* [Pipeline Steps](#pipeline-steps)
    1. [Calling Structural Variants](#somatic-sv)
* [Inputs](#Inputs)
* [Outputs](#outputs)
* [Testing and Validation](#testing-and-validation)
    * [Test Data Set](#test-data-set)
    * [Performance Validation](#performance-validation)
    * [Validation Tool](#validation-tool)


## Overview:
The call-sSV pipeline calls structural variants utilizing [Delly](https://github.com/dellytools/delly). This pipeline requires at least one tumor sample and a matched control sample.
This pipeline is developed using Nextflow , docker and can run either on a single node linux machine or a multi-node HPC cluster (e.g. Slurm, SGE). Additionally, it has been validated with the SMC-HET dataset.

## How to Run:

1.	Make sure the pipeline is already downloaded to your machine. You can either download the stable release or the dev version by cloning the repo.
2.	Update the nextflow.config file for input, output, and parameters. An example can be found [here](https://github.com/uclahs-cds/pipeline-call-sSV/blob/main/pipeline/config/nextflow.config). See [Inputs](#inputs) for the description of each variables in the config file.
3.	Update the input CSV. See [Inputs](#inputs) for the columns needed. All columns must exist to run the pipeline. An example can be found [here](https://github.com/uclahs-cds/pipeline-call-sSV/blob/main/pipeline/input/paired_turmor_control_samples.csv).
4.	See the submission script, [here](https://github.com/uclahs-cds/tool-submit-nf), to submit your pipeline

## Flow Diagram:

![](https://github.com/uclahs-cds/pipeline-call-sSV/blob/dev-doc-ghouse/call-sSV-workflow.PNG)

## Pipeline Steps:

#### 1. Calling Single Sample Somatic Structural Variants
The first step of the pipeline requires an aligned and sorted tumor sample .bam file and a matched control sample as an input for variant calling with Delly.

#### 2. Pre-filter Somatic Structural Variants
The second step of the pipeline applies soamtic pre-filtering to filter out any germline variants from the .bcf file generated in Step 1.

#### 3. Genotype pre-filtered somatic sites
For each tumor sample, genotype pre-filtered somatic sites across a larger panel of control samples to efficiently filter false postives and germline SVs. For performance reasons, this can be run in parallel for each sample of the control panel and you may want to combine multiple pre-filtered somatic site lists from multiple tumor samples.

#### 4. Post-filter 
Post-filter for somatic SVs using all control samples

## Inputs

### Input CSV

The input CSV should have all columns below and in the same order. An example of the input CSV can be found [here](https://github.com/uclahs-cds/pipeline-call-sSV/blob/main/pipeline/input/paired_turmor_control_samples.csv).

| Field |	Type |	Description |
|--- | --- | --- |
|tumor_sample_name |	string |	The tumor sample name to be passed to the Delly call.No white space allowed. |
|tumor_sample_bam	| string	| Absolute path to the tumor sample .BAM file. |
|control_sample_name	| string	| The control sample name to be passed. No white space allowed. |
|control_sample_bam |	string	| Absolute path to the control sample .BAM file |

## Nextflow Config File Parameters
| Input Parameter |	Required |	Type |	Description |
| ------- |   --------- | ------ | -------------|
| dataset_id |	yes	| string |	Boutros lab dataset id |
| blcds_registered_dataset	| yes |	Boolean | Affirms if dataset should be registered in the Boutros Lab Data registry. Default value is false. |
| sge_scheduler	| yes	| boolean	| Affirms whether job will be executed on the SGE cluster. Default value is false. |
| input_bams |	yes |	string	| Absolute path to the input CSV file for the pipeline. |
| reference_fasta	| yes |	path	| Absolute path to the reference genome fasta file. The reference genome is used by Delly for structural variant calling. |
| exclusion_file |	yes	| path |	Absolute path to the delly reference genome exclusion file utilized to remove suggested regions for structural variant calling. On Slurm/SGE , an HG38 exclusion file is located at /[hot\|data]/ref/hg38/delly/human.hg38.excl.tsv
| save_intermediate_files |	yes	| boolean |	Optional parameter to indicate whether intermediate files will be saved. Default value is true. |
| output_dir |	yes |	path |	Absolute path to the directory where the output files to be saved. |
| temp_dir	| yes	| path |	Absolute path to the directory where the nextflow’s intermediate files are saved. |
| skip_regenotype |	yes |	Boolean	| Paramter to conform whether Regenotype step should be carried out or not|

## Outputs

| Output |	Output type |	Description |
| ---- | ----- | -------- |
| .bcf |	final	| Binary VCF output format with somatic structural variants if found. |
| .bcf.csi	| final	| CSI-format index for BCF files |
| report.html, timeline.html and trace.txt	| Log |	A Nextflow report, timeline and trace files. |
| \*.log.command.*	| log |	Process and sample specific logging files created by nextflow. |
| *.sha512 |	checksum |	Generates SHA-512 hash to validate file integrity. |


## Testing and Validation

### Test Data Set

Testing was performed leveraging aligned and sorted bams generated using bwa-mem2-2.1 against reference GRCh38 (SMC-HET was aligned against hs37d5):

* A-mini: BWA-MEM2-2.1_TEST0000000_TWGSAMIN000001-T001-S01-F.bam and bai
* A-partial: BWA-MEM2-2.1_TEST0000000_TWGSAPRT000001-T001-S01-F.bam and bai
* A-full: a-full-CPCG0196-B1.bam and bai
* SMC-HET: HG002.N.bam and bai

### Test runs for the A-mini/partial/full samples were performed using the following reference files

* reference_fasta: /hot/ref/hg38/genome/genome.fa
* reference_fasta_index: /hot/ref/hg38/genome/genome.fa.fai
* exclusion_file: /hot/ref/hg38/delly/human.hg38.excl.tsv

## Performance Validation

Testing was performed primarily in the Boutros Lab SLURM Development cluster and the SLURM Covid cluster. Metrics below will be updated where relevant with additional testing and tuning outputs.

|Test Case	| Test Date	| Node Type |	Duration	| CPU Hours	| Virtual Memory Usage (RAM)-peak rss|
|----- | -------| --------| ----------| ---------| --------|
|A-mini	| 2021-05-12 |	F2 |	27m 47s	| 0.5h	| 1.8 GB |

## Validation Tool

Included is a template for validating your input files. For more information on the tool checkout:
https://github.com/uclahs-cds/tool-validate-nf