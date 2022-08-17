#!/usr/bin/env nextflow

nextflow.enable.dsl=2

import java.nio.file.Paths

log.info """\
======================================
C A L L - S S V   N F  P I P E L I N E
======================================
Boutros Lab

Current Configuration:
- pipeline:
    name: ${workflow.manifest.name}
    version: ${workflow.manifest.version}

- input:
    input_csv: "${params.input_csv}"
    reference_fasta: "${params.reference_fasta}"
    reference_fasta_index: "${params.reference_fasta}.fai"
    exclusion_file: "${params.exclusion_file}"
    filter_condition: "${params.filter_condition}"

- output:
    output_dir: "${params.output_dir}"
    log_output_dir: "${params.log_output_dir}"

- options:
    save_intermediate_files: ${params.save_intermediate_files}

- tools:
    delly: ${params.delly_version}

------------------------------------
Starting workflow...
------------------------------------
"""
.stripIndent()

include { run_validate_PipeVal } from './module/validation'
include { query_SampleName_BCFtools; filter_BCF_BCFtools } from './module/bcftools'
include { call_sSV_Delly; filter_sSV_Delly } from './module/delly'
include { generate_sha512 } from './module/sha512'

/**
* Check the params
*/

if (!params.input_csv){
    // error out - must provide an input_csv file
    error "***Error: You must specify an input_csv file***"
    }

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
* The input file params.input_csv looks as below:
* control_sample_bam, tumor_sample_bam
* /hot/users/ybugh/A-mini/0/output/HG002.N-0.bam, /hot/users/ybugh/A-mini/0/output/S2.T-0.bam
*
* Later, calling "delly call -g hg19.fa -v t1.pre.bcf -o geno.bcf -x hg19.excl tumor1.bam control1.bam ... controlN.bam" needs all the control samples, 
* which will be collected from params.input_control_bams.
*/

/**
* Create input_validation to validate the input bams
*/
input_validation = Channel
    .fromPath(params.input_csv, checkIfExists:true)
    .splitCsv(header:true)
    .map{
        row -> [
            row.tumor_sample_bam,
            row.control_sample_bam
            ]
        }
    .flatten()

if (params.verbose){
    input_validation.view()
    }

/**
* Create input_paired_bams_ch to get the paired turmor sample and control sample
*/
input_paired_bams_ch = Channel
    .fromPath(params.input_csv, checkIfExists:true)
    .splitCsv(header:true)
    .map{
        row -> tuple(
            Paths.get(row.tumor_sample_bam).getFileName().toString().split('.bam')[0],
            row.tumor_sample_bam,
            "${row.tumor_sample_bam}.bai",
            row.control_sample_bam,
            "${row.control_sample_bam}.bai"
            )
        }

if (params.verbose){
    input_paired_bams_ch.view()
    }

/**
* Create tumor_bams_ch to only get the turmor samples.
* I tried to reuse input_paired_bams_ch, however, in that way, I have to filter the paired control sample out of the all_control_samples_bams_list,
* otherwise, nextflow complains a same control sample is declared twice.
*/
tumor_bams_ch = Channel
    .fromPath(params.input_csv, checkIfExists:true)
    .splitCsv(header:true)
    .map{
        row -> tuple(
            Paths.get(row.tumor_sample_bam).getFileName().toString().split('.bam')[0],
            row.tumor_sample_bam,
            "${row.tumor_sample_bam}.bai"
            )
        }

if (params.verbose){
    tumor_bams_ch.view()
    }

workflow{
    /**
    * Validate the input bams
    */
    run_validate_PipeVal(input_validation)
    // Collect and store input validation output
    run_validate_PipeVal.out.val_file.collectFile(
      name: 'input_validation.txt',
      storeDir: "${params.output_dir}/validation"
      )

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
    query_SampleName_BCFtools(
        call_sSV_Delly.out.nt_call_bcf,
        call_sSV_Delly.out.samples
    )

    /**
    * Call "delly filter -f somatic -o t1.pre.bcf -s samples.tsv t1.bcf"
    * by using the call_sSV_Delly.out.samples and call_sSV_Delly.out.nt_call_bcf
    */
    filter_sSV_Delly(
        query_SampleName_BCFtools.out.samples,
        call_sSV_Delly.out.nt_call_bcf,
        call_sSV_Delly.out.nt_call_bcf_csi,
        call_sSV_Delly.out.tumor_sample_name
        )

    /**
    * Filter the output bcf from filter_sSV_Delly.
    * The default filter_condition is "FILTER=='PASS'", which filters out NonPass calls.
    */
    filter_BCF_BCFtools(filter_sSV_Delly.out.somatic_bcf, params.filter_condition)

    /**
    * Generate one sha512 checksum for the output files.
    */
    generate_sha512(filter_BCF_BCFtools.out.nonPassCallsFiltered.mix(filter_BCF_BCFtools.out.nonPassCallsFiltered_csi))

    }