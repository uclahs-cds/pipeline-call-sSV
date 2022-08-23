#!/usr/bin/env nextflow

log.info """\
------------------------------------
             D E L L Y
------------------------------------
Docker Images:
- docker_image_delly: ${params.docker_image_delly}
"""

include { generate_standard_filename } from '${projectDir}/external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

process call_sSV_Delly {
    container params.docker_image_delly

    publishDir "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
        enabled: params.save_intermediate_files,
        pattern: "DELLY-*.bcf*",
        mode: "copy"

    publishDir "$params.log_output_dir/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
        tuple(val(tumor_sample_name), path(tumor_sample_bam), path(tumor_sample_bai), path(control_sample_bam), path(control_sample_bai))
        path reference_fasta
        path reference_fasta_fai
        path exclusion_file

    output:
        path "${output_filename}.bcf", emit: nt_call_bcf
        path "${output_filename}.bcf.csi", emit: nt_call_bcf_csi
        path "${tumor_sample_name}", emit: samples
        path ".command.*"
        val tumor_sample_name, emit: tumor_sample_name

    script:
        output_filename = generate_standard_filename(
            "DELLY-${params.delly_version}",
            params.dataset_id,
            tumor_sample_name,
            [:]
            )
        """
        set -euo pipefail
        delly call \
            --genome "$reference_fasta" \
            --exclude "$exclusion_file" \
            --map-qual "${params.map_qual}" \
            --min-clique-size "${params.min_clique_size}" \
            --mad-cutoff "${params.mad_cutoff}" \
            --outfile "${output_filename}.bcf" \
            "$tumor_sample_bam" \
            "$control_sample_bam"

        touch "${tumor_sample_name}"
        """
    }

process filter_sSV_Delly {
    container params.docker_image_delly

    publishDir "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
        enabled: params.save_intermediate_files,
        pattern: "${filename_base}_somatic.bcf*",
        mode: "copy"

    publishDir "$params.log_output_dir/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
        path samples
        path input_bcf
        path input_bcf_csi
        val tumor_sample_name

    output:
        path "${filename_base}_somatic.bcf", emit: somatic_bcf
        path "${filename_base}_somatic.bcf.csi", emit: somatic_bcf_csi
        path ".command.*"
 
    script:
        filename_base = file(input_bcf).baseName
        """
        set -euo pipefail
        delly filter -f somatic -s ${samples} -o ${filename_base}_somatic.bcf "$input_bcf"
        """
    }
