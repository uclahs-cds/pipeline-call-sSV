#!/usr/bin/env nextflow

log.info """\
------------------------------------
             M A N T A
------------------------------------
Docker Images:
- docker_image_manta: ${params.docker_image_manta}
"""

include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

process call_sSV_Manta {
    container params.docker_image_manta

    publishDir "$params.output_dir/${params.docker_image_manta.split("/")[1].replace(':', '-').capitalize()}/output",
        pattern: "*vcf.gz*",
        mode: "copy"

    publishDir "$params.output_dir/${params.docker_image_manta.split("/")[1].replace(':', '-').capitalize()}/QC",
        pattern: "*Stats*",
        mode: "copy"

    publishDir "$params.log_output_dir/process-log/${params.docker_image_manta.split("/")[1].replace(':', '-').capitalize()}",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}-${task.index}/log${file(it).getName()}" }

    input:
        tuple(val(tumor_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai))
        path reference_fasta
        path reference_fasta_fai

    output:
        path("MantaWorkflow/results/variants/candidateSmallIndels.vcf.gz"), emit: vcf_small_indel_sv_file
        path("MantaWorkflow/results/variants/candidateSmallIndels.vcf.gz.tbi"), emit: vcf_small_indel_sv_tbi
        path("MantaWorkflow/results/variants/diploidSV.vcf.gz"), emit: vcf_diploid_sv_file
        path("MantaWorkflow/results/variants/diploidSV.vcf.gz.tbi"), emit: vcf_diploid_sv_tbi
        path("MantaWorkflow/results/variants/candidateSV.vcf.gz"), emit: vcf_candidate_sv_file
        path("MantaWorkflow/results/variants/candidateSV.vcf.gz.tbi"), emit: vcf_candidate_sv_tbi
        path "*vcf.gz*"
        path "*Stats*"
        path ".command.*"
        val tumor_id, emit: tumor_id
    
    script:
        output_filename = generate_standard_filename(
            "Manta-${params.manta_version}",
            params.dataset_id,
            tumor_id,
            [:]
            )
        """
        set -euo pipefail
        configManta.py \
            --normalBam "$normal_bam" \
            --tumorBam "$tumor_bam"
            --referenceFasta "$reference_fasta" \
            --runDir MantaWorkflow
        MantaWorkflow/runWorkflow.py

        cp MantaWorkflow/results/variants/* ./
        cp MantaWorkflow/results/stats/* ./
        """
    }
