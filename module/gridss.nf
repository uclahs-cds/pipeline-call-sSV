#!/usr/bin/env nextflow

log.info """\
------------------------------------
             G R I D S S 2
------------------------------------
Docker Images:
- docker_image_gridss: ${params.docker_image_gridss}
"""

include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

process preprocess_BAM_GRIDSS {
    container params.docker_image_gridss

    publishDir "${params.workflow_output_dir}/output",
        pattern: "*vcf.gz*",
        mode: "copy"

    publishDir "${params.workflow_output_dir}/QC",
        pattern: "*Stats*",
        mode: "copy"

    publishDir "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
        tuple(val(sample_id), path(sample_bam), path(sample_bai))
        path reference_fasta
        path reference_fasta_fai
        path gridss_jar

    output:
        path "*.bam", emit: gridss_preprocess_bam
        path ".command.*"
        val sample_id, emit: sample_id

    script:
        output_filename = generate_standard_filename(
            "GRIDSS2-${params.gridss_version}",
            params.dataset_id,
            tumor_id,
            [:]
            )
        """
        set -euo pipefail
        gridss \
            -r ${reference_fasta}
            -j ${gridss_jar}
            -s preprocess \
            -t ${task.cpus} \
            ${bam}
        """
    }
