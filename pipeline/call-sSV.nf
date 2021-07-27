#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import java.nio.file.Paths;

log.info """\
======================================
C A L L - S S V   N F  P I P E L I N E
======================================
Boutros Lab

Current Configuration:
- pipeline:
    mainScript: ${workflow.manifest.mainScript}
    version: ${workflow.manifest.version}

- input:
    input_paired_bams: "${params.input_paired_bams}"
    input_control_bams = "${params.input_control_bams}"
    sample_types = "${params.sample_types}"
    reference_fasta: "${params.reference_fasta}"
    reference_fasta_index: "${params.reference_fasta}.fai"
    exclusion_file: "${params.exclusion_file}"

- output:
    output_dir: "${params.output_dir}"
    output_log_dir: "${params.output_log_dir}"
    temp_dir: "${params.temp_dir}"

- options:
    save_intermediate_files: ${params.save_intermediate_files}

- tools:
    delly: ${params.delly_version}

------------------------------------
Starting workflow...
------------------------------------
"""
.stripIndent()

include { validate_file } from './modules/validation'
include { query_sample_name_Bcftools } from './modules/bcftools'
include { call_sSV_Delly; regenotype_sSV_Delly; filter_sSV_Delly as filter_sSV_Delly_initialCall; filter_sSV_Delly as filter_sSV_Delly_regenotyped } from './modules/delly'
include { generate_sha512 } from './modules/sha512'

/**
* Check the params
*/

if (!params.reference_fasta){
    // error out - must provide a reference FASTA file
    error "***Error: You must specify a reference FASTA file***"
    }

if (!params.exclusion_file){
    // error out - must provide exclusion file
    error "*** Error: You must provide an exclusion file***"
    }

reference_fasta_index = "${params.reference_fasta}.fai"

/**
* The input file params.input_paired_bams looks as below:
* control_sample_bam, tumor_sample_bam
* /hot/users/ybugh/A-mini/0/output/HG002.N-0.bam, /hot/users/ybugh/A-mini/0/output/S2.T-0.bam
*
* Later, calling "delly call -g hg19.fa -v t1.pre.bcf -o geno.bcf -x hg19.excl tumor1.bam control1.bam ... controlN.bam" needs all the control samples, 
* which will be collected from params.input_control_bams.
*/

/**
* Create validation_channel to validate the input bams
*/
validation_channel = Channel
    .fromPath(params.input_paired_bams, checkIfExists:true)
    .splitCsv(header:true)
    .map{
        row -> [
            row.tumor_sample_bam,
            row.control_sample_bam
            ]
        }
    .flatten()

//validation_channel.view()

/**
* Create input_paired_bams_ch to get the paired turmor sample and control sample
*/
input_paired_bams_ch = Channel
    .fromPath(params.input_paired_bams, checkIfExists:true)
    .splitCsv(header:true)
    .map{
        row -> tuple(
            Paths.get(row.tumor_sample_bam).getFileName().toString().split('.bam')[0],
            row.tumor_sample_bam,
            "${row.tumor_sample_bam}.bai",
            Paths.get(row.control_sample_bam).getFileName().toString().split('.bam')[0],
            row.control_sample_bam,
            "${row.control_sample_bam}.bai"
            )
        }

input_paired_bams_ch.view()

/**
* Create tumor_bams_ch to only get the turmor samples.
* I tried to reuse input_paired_bams_ch, however, in that way, I have to filter the paired control sample out of the all_control_samples_bams_list,
* otherwise, nextflow complains a same control sample is declared twice.
*/
tumor_bams_ch = Channel
    .fromPath(params.input_paired_bams, checkIfExists:true)
    .splitCsv(header:true)
    .map{
        row -> tuple(
            Paths.get(row.tumor_sample_bam).getFileName().toString().split('.bam')[0],
            row.tumor_sample_bam,
            "${row.tumor_sample_bam}.bai"
            )
        }

tumor_bams_ch.view()

/**
* Create all_control_samples_bams_bais_ch.
* this is used to declare these files in dockers.
*/
all_control_samples_bams_bais_list = Channel
    .fromPath(params.input_control_bams, checkIfExists:true)
    .splitCsv(header:true)
    .map{
        row -> [
            row.control_sample_bam,
            "${row.control_sample_bam}.bai"
            ]
        }
    .flatten()
    .toList()

all_control_samples_bams_bais_list.view()

workflow{
    /**
    * Validate the input bams
    */
    validate_file(validation_channel)

    /**
    * Call "delly call -x hg19.excl -o t1.bcf -g hg19.fa tumor1.bam control1.bam" per paired (tumor sample, control sample)
    * The sv are stored in call_sSV_Delly.out.nt_call_bcf
    * also create call_sSV_Delly.out.samples per paired (tumor sample, control sample)
    */
    call_sSV_Delly(
        input_paired_bams_ch,
        params.reference_fasta,
        reference_fasta_index,
        params.exclusion_file
        )

    /**
    * calling "delly filter -f somatic -s samples.tsv -o t1.pre.bcf t1.bcf" requires a samples.tsv, which should look like:
    * sample_name   sample_type
    * S2_v1.1.5	tumor
    * HG002.N	control
    * 
    * Use bcftools query -l to get the sample names out of call_sSV_Delly.out.nt_call_bcf
    * Further generate ${control_sample_bam_name}_${tumor_sample_bam_name}_samples.tsv which will be used by delly filter
    * Note, the order of samples in call_sSV_Delly.out.nt_call_bcf is determined by the order of samples in delly call.
    * For example, 
    *    delly call \
    *    -g /tmp/ref/genome/genome.fa \
    *    -x /tmp/ref/delly/human.hg38.excl.tsv \
    *    -o /tmp/output/output.bcf \
    *    /tmp/bams/HG002.N-0.bam \
    *    /tmp/bams/S2.T-0.bam"
    * bcftools query -l /tmp/output/output.bcf will yield
    * HG002.N
    * S2_v1.1.5
    * If you put /tmp/bams/S2.T-0.bam in front of /tmp/bams/HG002.N-0.bam, bcftools query -l /tmp/output/output.bcf will yield
    * S2_v1.1.5 
    * HG002.N
    */
    query_sample_name_Bcftools(
        call_sSV_Delly.out.nt_call_bcf,
        call_sSV_Delly.out.samples,
        params.sample_types
    )

    /**
    * Call "delly filter -f somatic -o t1.pre.bcf -s samples.tsv t1.bcf"
    * by using the call_sSV_Delly.out.samples and call_sSV_Delly.out.nt_call_bcf
    */
    filter_sSV_Delly_initialCall(
        query_sample_name_Bcftools.out.samples,
        call_sSV_Delly.out.nt_call_bcf,
        call_sSV_Delly.out.nt_call_bcf_csi,
        params.SINGLE_CTRL_SAMPLE
        )

    if (!params.skip_regenotype) {
        /**
        * Genotype pre-filtered somatic sites across a larger panel of control samples.
        * If something is being seen in all samples then it's more probable that it's a false positive.
        */
        regenotype_sSV_Delly(
            tumor_bams_ch,
            params.reference_fasta,
            reference_fasta_index,
            params.exclusion_file,
            all_control_samples_bams_bais_list,
            filter_sSV_Delly_initialCall.out.filtered_somatic_bcf
        )

        /**
        * Call filter_sSV_Delly again to filter out germline SVs.
        */
        filter_sSV_Delly_regenotyped(
            query_sample_name_Bcftools.out.samples,
            regenotype_sSV_Delly.out.nt_regenotype_bcf,
            regenotype_sSV_Delly.out.nt_regenotype_bcf_csi,
            params.ALL_CTRL_SAMPLES
            )
    }

    /**
    * Generate sha512 checksum for the filtered_somatic_AllCtrlSamples.bcf.
    */
    generate_sha512(filter_sSV_Delly_regenotyped.out.filtered_somatic_bcf)
    }
