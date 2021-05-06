#!/usr/bin/env nextflow

nextflow.enable.dsl=2

log.info """\
======================================
C A L L - G S V   N F  P I P E L I N E
======================================
Boutros Lab

Current Configuration:
- input:
    input_bams_bcfs: "${params.input_bams}"
    reference_fasta: "${params.reference_fasta}"
    reference_fasta_index: "${params.reference_fasta}.fai"
    reference_prefix: "${params.reference_prefix}"
    exclusion_file: "${params.exclusion_file}"

- output:
    output_dir: "${params.output_dir}"
    output_log_dir: "${params.output_log_dir}"
    temp_dir: "${params.temp_dir}"

- options:
    save_intermediate_files: ${params.save_intermediate_files}

- tools:
    delly: ${params.delly_version}
    bcftools: ${params.bcftools_version}

------------------------------------
Starting workflow...
------------------------------------
"""
.stripIndent()

include { validate_file } from './modules/validation'
include { delly_call_NT } from './modules/delly'
//include { delly_call_sv; delly_filter_sv } from './modules/delly'

/**
* Check the params
*/

if (!params.reference_fasta) {
    // error out - must provide a reference FASTA file
    error "***Error: You must specify a reference FASTA file***"
    }

if (!params.exclusion_file) {
    // error out - must provide exclusion file
    error "*** Error: You must provide an exclusion file***"
    }

reference_fasta_index = "${params.reference_fasta}.fai"

/**
* The input file "paired_turmor_control_samples.csv" looks as below:
* tumor_sample_name,tumor_sample_bam,control_sample_name,control_sample_bam
* S2_v1.1.5,/hot/users/ybugh/A-mini/0/output/S2.T-0.bam,HG002.N,/hot/users/ybugh/A-mini/0/output/HG002.N-0.bam
*
* Later, calling "delly filter -f somatic -s samples.tsv -o t1.pre.bcf t1.bcf" requires a samples.tsv, which should look like:
* HG002.N	control
* S2_v1.1.5	tumor	/hot/users/ybugh/A-mini/0/output/S2.T-0.bam
* 
* The pipeline will create such samples.tsv from the "paired_turmor_control_samples.csv" on the fly.
*
* Also calling "delly call -g hg19.fa -v t1.pre.bcf -o geno.bcf -x hg19.excl tumor1.bam control1.bam ... controlN.bam" needs all the control samples, 
* which will be collected from the "paired_turmor_control_samples.csv" too.
*/

/**
* Create validation_channel to validate the input bams
* tumor_sample_name,tumor_sample_bam,control_sample_name,control_sample_bam
*/
validation_channel = Channel
    .fromPath(params.input_bams, checkIfExists:true)
    .splitCsv(header:true)
    .map{
        row -> [
            row.tumor_sample_bam,
            row.control_sample_bam
            ]
        }
    .flatten()

validation_channel.view()

input_bams_ch = Channel
    .fromPath(params.input_bams, checkIfExists:true)
    .splitCsv(header:true)
    .map {
        row -> tuple(
            row.tumor_sample_name,
            row.tumor_sample_bam,
            "${row.tumor_sample_bam}.bai",
            row.control_sample_name,
            row.control_sample_bam,
            "${row.control_sample_bam}.bai"
            )
        }

input_bams_ch.view()

/**
* input_sv_bcfs_ch contains the list of original SV bcfs that will be merged to get the unified sites.
*/
/*
input_sv_bcfs_ch_toList = input_sv_bcfs_ch.toList()
*/

workflow {
    validate_file(validation_channel)
    delly_call_NT(input_bams_ch, params.reference_fasta, reference_fasta_index, params.exclusion_file, params.SINGLE_CTRL_SAMPLE)
    } 
