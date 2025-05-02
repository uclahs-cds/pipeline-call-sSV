#!/usr/bin/env nextflow

log.info """\
------------------------------------
             G R I D S S 2
------------------------------------
Docker Images:
- docker_image_gridss2: ${params.docker_image_gridss2}
"""

include { generate_standard_filename; sanitize_string } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

process preprocess_BAM_GRIDSS2 {
    container params.docker_image_gridss2

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
        path(gridss2_reference_fasta)
        path(gridss2_reference_files)

    output:
        path "${bam_name}.gridss.working/*", emit: gridss2_preprocess
        path ".command.*"

    script:
        gridss2_mem = "${task.memory.toGiga()}g"
        gridss2_jar = "/usr/local/share/gridss-${params.gridss2_version}-1/gridss.jar"
        bam_name = file(sample_bam).getName()
        output_filename = generate_standard_filename(
            "GRIDSS2-${params.gridss2_version}",
            params.dataset_id,
            sample_id,
            [:]
            )

        """
        set -euo pipefail
        gridss \
            -r ${gridss2_reference_fasta} \
            -j ${gridss2_jar} \
            -s preprocess \
            -t ${task.cpus} \
            --jvmheap ${gridss2_mem} \
            ${sample_bam}
        """
    }

process run_assembly_GRIDSS2 {
    container params.docker_image_gridss2

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
        path(gridss2_preprocess_dir)
        path(gridss2_reference_fasta)
        path(gridss2_reference_files)
        path(gridss2_blacklist)

    output:
        path "${tumor_id}.assembly.bam", emit: gridss2_assembly_bam
        path "${tumor_id}.assembly.bam.gridss.working/*", emit: gridss2_assembly
        path ".command.*"

    script:
        otherjvmheap = params.other_jvm_heap
        gridss2_otherjvmheap = "${otherjvmheap.toGiga()}g"
        gridss2_jvmheap = "${(task.memory - otherjvmheap).toGiga()}g"
        gridss2_jar = "/usr/local/share/gridss-${params.gridss2_version}-1/gridss.jar"
        output_filename = generate_standard_filename(
            "GRIDSS2-${params.gridss2_version}",
            params.dataset_id,
            tumor_id,
            [:]
            )

        """
        set -euo pipefail
        gridss \
            -r ${gridss2_reference_fasta} \
            -j ${gridss2_jar} \
            -s assemble \
            -t ${task.cpus} \
            --jvmheap ${gridss2_jvmheap} \
            --otherjvmheap ${gridss2_otherjvmheap} \
            -b ${gridss2_blacklist} \
            -a ${tumor_id}.assembly.bam \
            ${normal_bam} \
            ${tumor_bam}
        """
    }

process call_sSV_GRIDSS2 {
    container params.docker_image_gridss2

    publishDir "${params.workflow_output_dir}/output/",
        pattern: "${output_filename}.{vcf,vcf.idx}",
        mode: "copy"

    publishDir "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
        enabled: params.save_intermediate_files,
        pattern: "${output_filename}.vcf.gridss.working/*",
        mode: "copy"

    publishDir "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
        tuple(val(tumor_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai))
        path(gridss2_preprocess_dir)
        path(gridss2_assembly_dir)
        path(gridss2_assembly_bam)
        path(gridss2_reference_fasta)
        path(gridss2_reference_files)
        path(gridss2_blacklist)

    output:
        path "${output_filename}.vcf", emit: gridss2_vcf
        path "${output_filename}.vcf.idx", emit: gridss2_vcf_idx
        path "${output_filename}.vcf.gridss.working/*", emit: gridss2_vcf_dir
        path ".command.*"

    script:
        otherjvmheap = params.other_jvm_heap
        gridss2_otherjvmheap = "${otherjvmheap.toGiga()}g"
        gridss2_jvmheap = "${(task.memory - otherjvmheap).toGiga()}g"
        gridss2_jar = "/usr/local/share/gridss-${params.gridss2_version}-1/gridss.jar"
        output_filename = generate_standard_filename(
            "GRIDSS2-${params.gridss2_version}",
            params.dataset_id,
            tumor_id,
            [:]
            )

        """
        set -euo pipefail
        gridss \
            -r ${gridss2_reference_fasta} \
            -j ${gridss2_jar} \
            -s call \
            -t ${task.cpus} \
            --jvmheap ${gridss2_jvmheap} \
            --otherjvmheap ${gridss2_otherjvmheap} \
            -b ${gridss2_blacklist} \
            -a ${gridss2_assembly_bam} \
            --output ${output_filename}.vcf \
            ${normal_bam} \
            ${tumor_bam}
        """
    }

process filter_sSV_GRIDSS2 {
    container params.docker_image_gridss2

    publishDir "${params.workflow_output_dir}/output/",
        pattern: "*.vcf.bgz*",
        mode: "copy"

    publishDir "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
        val(tumor_id)
        path(gridss2_vcf)
        path(gridss2_pon_dir)

    output:
        path "*vcf.bgz*", emit: gridss2_filter_vcf_files
        path ".command.*"

    script:
        gridss2_script_dir = "/usr/local/share/gridss-${params.gridss2_version}-1"
        output_filename = generate_standard_filename(
            "GRIDSS2-${params.gridss2_version}",
            params.dataset_id,
            tumor_id,
            [:]
            )

        """
        set -euo pipefail
        ${gridss2_script_dir}/gridss_somatic_filter \
            --pondir ${gridss2_pon_dir} \
            --scriptdir ${gridss2_script_dir}/ \
            --input ${gridss2_vcf} \
            --output ${output_filename}_high-confidence-somatic.vcf \
            --fulloutput ${output_filename}_high-low-confidence-somatic.vcf \
            --normalordinal 1 \
            --tumourordinal 2
        """
    }
