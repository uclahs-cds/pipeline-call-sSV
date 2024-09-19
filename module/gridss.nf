#!/usr/bin/env nextflow

log.info """\
------------------------------------
             G R I D S S 2
------------------------------------
Docker Images:
- docker_image_gridss: ${params.docker_image_gridss}
"""

include { generate_standard_filename; sanitize_string } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

process preprocess_BAM_GRIDSS {
    container params.docker_image_gridss

    publishDir "${params.workflow_output_dir}/intermediate",
        pattern: "${bam_name}.gridss.working/*",
        mode: "copy",
        saveAs: {
            "${output_filename}.${sanitize_string(file(it).getName().replace("${bam_name}.", ""))}"
            }

    publishDir "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
        tuple(val(sample_id), path(sample_bam), path(sample_index))
        path(gridss_reference_fasta)
        path(gridss_reference_files)

    output:
        path "${bam_name}.gridss.working/*", emit: gridss_preprocess
        path ".command.*"

    script:
        gridss_mem = "${task.memory.toGiga()}g"
        gridss_jar = "/usr/local/share/gridss-${params.gridss_version}-1/gridss.jar"
        bam_name = file(sample_bam).getName()
        output_filename = generate_standard_filename(
            "GRIDSS2-${params.gridss_version}",
            params.dataset_id,
            sample_id,
            [:]
            )

        """
        set -euo pipefail
        gridss \
            -r ${gridss_reference_fasta} \
            -j ${gridss_jar} \
            -s preprocess \
            -t ${task.cpus} \
            --jvmheap ${gridss_mem} \
            ${sample_bam}
        """
    }
