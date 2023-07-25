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
    samples: ${params.samples_to_process}
    reference_fasta: "${params.reference_fasta}"
    reference_fasta_index: "${params.reference_fasta}.fai"
    exclusion_file: "${params.exclusion_file}"
    filter_condition: "${params.filter_condition}"

- output:
    output_dir: "${params.output_dir}"
    log_output_dir: "${params.log_output_dir}"

- options:
    save_intermediate_files: ${params.save_intermediate_files}
    SV caller(s): ${params.algorithm}

- tools:
    DELLY: ${params.delly_version}
    BCFtools: ${params.bcftools_version}
    Manta: ${params.manta_version}
    PipeVal: ${params.pipeval_version}

------------------------------------
Starting workflow...
------------------------------------
"""
.stripIndent()

include { run_validate_PipeVal } from "./external/pipeline-Nextflow-module/modules/PipeVal/validate/main.nf" addParams(
    options: [ docker_image_version: params.pipeval_version ]
    )
include { query_SampleName_BCFtools; filter_BCF_BCFtools } from './module/bcftools' addParams(
    workflow_output_dir: "${params.output_dir_base}/DELLY-${params.delly_version}"
    )
include { call_sSV_Delly; filter_sSV_Delly } from './module/delly' addParams(
    workflow_output_dir: "${params.output_dir_base}/DELLY-${params.delly_version}"
    )
include { call_sSV_Manta } from './module/manta' addParams(
    workflow_output_dir: "${params.output_dir_base}/Manta-${params.manta_version}"
    )
include { generate_sha512 as generate_sha512_BCFtools } from './module/sha512' addParams(
    workflow_output_dir: "${params.output_dir_base}/DELLY-${params.delly_version}"
    )
include { generate_sha512 as generate_sha512_Manta } from './module/sha512' addParams(
    workflow_output_dir: "${params.output_dir_base}/Manta-${params.manta_version}"
    )

// Returns the index file for the given bam
def indexFile(bam) {
    if (bam.endsWith('.bam')) {
        return "${bam}.bai"
    } else {
        throw new Exception("Index file for ${bam} file type not supported. Use .bam!")
    }
}

Channel.from(params.samples_to_process)
    .map{ sample -> ['index': indexFile(sample.path)] + sample }
    .set{ input_ch_samples_with_index }

input_ch_samples_with_index
    .map{ sample -> [sample.path, sample.index] }
    .flatten()
    .set{ input_validation }
if (params.verbose){
    input_validation.view()
    }

tumor_id_bam_bai = input_ch_samples_with_index
    .filter{ it.sample_type == 'tumor' }
    .map{ it -> [it.id, it.path, it.index] }
normal_bam_bai = input_ch_samples_with_index
    .filter{ it.sample_type == 'normal' }
    .map{ it -> [it.path, it.index] }

input_paired_bams_ch = tumor_id_bam_bai.combine(normal_bam_bai)
if (params.verbose){
    input_paired_bams_ch.view()
    }

tumor_bams_ch = tumor_id_bam_bai
if (params.verbose){
    tumor_bams_ch.view()
    }

reference_fasta_index = "${params.reference_fasta}.fai"

workflow {
    /**
    * Validate the input bams
    */
    run_validate_PipeVal(input_validation)
    // Collect and store input validation output
    run_validate_PipeVal.out.validation_result.collectFile(
        name: 'input_validation.txt',
        storeDir: "${params.output_dir_base}/validation/run_validate_PipeVal"
        )

    /**
    * Call "delly call -x hg19.excl -o t1.bcf -g hg19.fa tumor1.bam normal1.bam" per paired (tumor sample, normal sample)
    * The sv are stored in call_sSV_Delly.out.nt_call_bcf
    * also create call_sSV_Delly.out.samples per paired (tumor sample, normal sample)
    */
    if ('delly' in params.algorithm) {
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
        * HG002.N	normal
        *
        * Use bcftools query -l to get the sample names out of call_sSV_Delly.out.nt_call_bcf
        * Further generate BCFtools_${params.bcftools_version}_${params.dataset_id}_${tumor_id}_query-tumor-normal-name.tsv which will be used by delly filter
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
            call_sSV_Delly.out.samples,
            call_sSV_Delly.out.tumor_id
        )

        /**
        * Call "delly filter -f somatic -o t1.pre.bcf -s samples.tsv t1.bcf"
        * by using the call_sSV_Delly.out.samples and call_sSV_Delly.out.nt_call_bcf
        */
        filter_sSV_Delly(
            query_SampleName_BCFtools.out.samples,
            call_sSV_Delly.out.nt_call_bcf,
            call_sSV_Delly.out.nt_call_bcf_csi,
            call_sSV_Delly.out.tumor_id
            )

        /**
        * Filter the output bcf from filter_sSV_Delly.
        * The default filter_condition is "FILTER=='PASS'", which filters out NonPass calls.
        */
        filter_BCF_BCFtools(
            filter_sSV_Delly.out.somatic_bcf,
            params.filter_condition,
            call_sSV_Delly.out.tumor_id
            )

        /**
        * Generate one sha512 checksum for DELLY's output files.
        */
        generate_sha512_BCFtools(
            filter_BCF_BCFtools.out.nonPassCallsFiltered_and_csi.flatten()
            )
        }
    if ('manta' in params.algorithm) {
        call_sSV_Manta(
            input_paired_bams_ch,
            params.reference_fasta,
            reference_fasta_index
            )
        /**
        * Generate one sha512 checksum for the output files.
        */
        generate_sha512_Manta(
            call_sSV_Manta.out.manta_vcfs.flatten()
            )
        }
    }
