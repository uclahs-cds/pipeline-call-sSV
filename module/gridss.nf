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

    publishDir "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
        enabled: params.save_intermediate_files,
        pattern: "${bam_name}.gridss.working/*",
        mode: "copy",
        saveAs: {
            "${output_filename}.gridss.working/${output_filename}.${sanitize_string(file(it).getName().replace("${bam_name}.", ""))}"
            }

    publishDir "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/${task.process}-${task.index}/log${file(it).getName()}" }

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

process run_assembly_GRIDSS {
    container params.docker_image_gridss

    publishDir "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
        enabled: params.save_intermediate_files,
        pattern: "${tumor_id}.assembly.bam",
        mode: "copy",
        saveAs: {
            "${output_filename}_${sanitize_string(file(it).getName().replace("${tumor_id}.", ""))}"
            }

    publishDir "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
        enabled: params.save_intermediate_files,
        pattern: "${tumor_id}.assembly.bam.gridss.working/*",
        mode: "copy",
        saveAs: {
            "${output_filename}.assembly.bam.gridss.working/${output_filename}_${sanitize_string(file(it).getName().replace("${tumor_id}.", ""))}"
            }

    publishDir "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
        tuple(val(tumor_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai))
        path(gridss_preprocess_dir)
        path(gridss_reference_fasta)
        path(gridss_reference_files)
        path(gridss_blacklist)

    output:
        path "${tumor_id}.assembly.bam", emit: gridss_assembly_bam
        path "${tumor_id}.assembly.bam.gridss.working/*", emit: gridss_assembly
        path ".command.*"

    script:
        otherjvmheap = params.other_jvm_heap
        gridss_otherjvmheap = "${otherjvmheap.toGiga()}g"
        gridss_jvmheap = "${(task.memory - otherjvmheap).toGiga()}g"
        gridss_jar = "/usr/local/share/gridss-${params.gridss_version}-1/gridss.jar"
        output_filename = generate_standard_filename(
            "GRIDSS2-${params.gridss_version}",
            params.dataset_id,
            tumor_id,
            [:]
            )

        """
        set -euo pipefail
        gridss \
            -r ${gridss_reference_fasta} \
            -j ${gridss_jar} \
            -s assemble \
            -t ${task.cpus} \
            --jvmheap ${gridss_jvmheap} \
            --otherjvmheap ${gridss_otherjvmheap} \
            -b ${gridss_blacklist} \
            -a ${tumor_id}.assembly.bam \
            ${normal_bam} \
            ${tumor_bam}
        """
    }

process call_sSV_GRIDSS {
    container params.docker_image_gridss

    publishDir "${params.workflow_output_dir}/output/",
        pattern: "${tumor_id}.{vcf,vcf.idx}",
        mode: "copy",
        saveAs: {
            "${output_filename}.${sanitize_string(file(it).getName().replace("${tumor_id}.", ""))}"
            }

    publishDir "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
        enabled: params.save_intermediate_files,
        pattern: "${tumor_id}.vcf.gridss.working/*",
        mode: "copy",
        saveAs: {
            "${output_filename}.vcf.gridss.working/${output_filename}.${sanitize_string(file(it).getName().replace("${tumor_id}.", ""))}"
            }

    publishDir "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
        tuple(val(tumor_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai))
        path(gridss_preprocess_dir)
        path(gridss_assembly_dir)
        path(gridss_assembly_bam)
        path(gridss_reference_fasta)
        path(gridss_reference_files)
        path(gridss_blacklist)

    output:
        path "${tumor_id}.vcf", emit: gridss_vcf
        path "${tumor_id}.vcf.idx", emit: gridss_vcf_idx
        path "${tumor_id}.vcf.gridss.working/*", emit: gridss_vcf_dir
        path ".command.*"

    script:
        otherjvmheap = params.other_jvm_heap
        gridss_otherjvmheap = "${otherjvmheap.toGiga()}g"
        gridss_jvmheap = "${(task.memory - otherjvmheap).toGiga()}g"
        gridss_jar = "/usr/local/share/gridss-${params.gridss_version}-1/gridss.jar"
        output_filename = generate_standard_filename(
            "GRIDSS2-${params.gridss_version}",
            params.dataset_id,
            tumor_id,
            [:]
            )

        """
        set -euo pipefail
        gridss \
            -r ${gridss_reference_fasta} \
            -j ${gridss_jar} \
            -s call \
            -t ${task.cpus} \
            --jvmheap ${gridss_jvmheap} \
            --otherjvmheap ${gridss_otherjvmheap} \
            -b ${gridss_blacklist} \
            -a ${gridss_assembly_bam} \
            --output ${tumor_id}.vcf \
            ${normal_bam} \
            ${tumor_bam}
        """
    }

process filter_sSV_GRIDSS {
    container params.docker_image_gridss

    publishDir "${params.workflow_output_dir}/output/",
        pattern: "${tumor_id}.{vcf,vcf.idx}",
        mode: "copy",
        saveAs: {
            "${output_filename}.${sanitize_string(file(it).getName().replace("${tumor_id}.", ""))}"
            }

    publishDir "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
        enabled: params.save_intermediate_files,
        pattern: "${tumor_id}.vcf.gridss.working/*",
        mode: "copy",
        saveAs: {
            "${output_filename}.vcf.gridss.working/${output_filename}.${sanitize_string(file(it).getName().replace("${tumor_id}.", ""))}"
            }

    publishDir "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
        tuple(val(tumor_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai))
        path(gridss_preprocess_dir)
        path(gridss_assembly_dir)
        path(gridss_assembly_bam)
        path(gridss_reference_fasta)
        path(gridss_reference_files)
        path(gridss_blacklist)

    output:
        path "${tumor_id}.vcf", emit: gridss_vcf
        path "${tumor_id}.vcf.idx", emit: gridss_vcf_idx
        path "${tumor_id}.vcf.gridss.working/*", emit: gridss_vcf_dir
        path ".command.*"

    script:
        otherjvmheap = params.other_jvm_heap
        gridss_otherjvmheap = "${otherjvmheap.toGiga()}g"
        gridss_jvmheap = "${(task.memory - otherjvmheap).toGiga()}g"
        gridss_jar = "/usr/local/share/gridss-${params.gridss_version}-1/gridss.jar"
        output_filename = generate_standard_filename(
            "GRIDSS2-${params.gridss_version}",
            params.dataset_id,
            tumor_id,
            [:]
            )

        """
        set -euo pipefail
        /usr/local/share/gridss-2.13.2-1/gridss_somatic_filter \
            --pondir ${params.gridss_pon_dir} \
            --scriptdir /usr/local/share/gridss-2.13.2-1/ \
            --input \
            --output \
            --fulloutput \
            -n 1 \
            -t 2 \
            -j ${gridss_jar} \
            -t ${task.cpus} \
            --jvmheap ${gridss_jvmheap} \
            --otherjvmheap ${gridss_otherjvmheap} \
            -b ${gridss_blacklist} \
            -a ${gridss_assembly_bam} \
            --output ${tumor_id}.vcf \
            ${normal_bam} \
            ${tumor_bam}
        """
    }
