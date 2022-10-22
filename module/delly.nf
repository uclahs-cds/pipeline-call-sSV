#!/usr/bin/env nextflow

log.info """\
------------------------------------
             D E L L Y
------------------------------------
Docker Images:
- docker_image_delly: ${params.docker_image_delly}
"""

include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

process call_sSV_Delly {
    container params.docker_image_delly

    publishDir "$params.output_dir_base/intermediate/${task.process.replace(':', '/')}",
        enabled: params.save_intermediate_files,
        pattern: "DELLY-*.bcf*",
        mode: "copy"

    publishDir "$params.log_output_dir/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
        tuple(val(tumor_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai))
        path reference_fasta
        path reference_fasta_fai
        path exclusion_file

    output:
        path "${output_filename}.bcf", emit: nt_call_bcf
        path "${output_filename}.bcf.csi", emit: nt_call_bcf_csi
        path "${tumor_id}", emit: samples
        path ".command.*"
        val tumor_id, emit: tumor_id

    script:
        output_filename = generate_standard_filename(
            "DELLY-${params.delly_version}",
            params.dataset_id,
            tumor_id,
            [additional_information: "unfiltered"]
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
            "$tumor_bam" \
            "$normal_bam"

        touch "${tumor_id}"
        """
    }

process filter_sSV_Delly {
    container params.docker_image_delly

    publishDir "$params.output_dir_base/intermediate/${task.process.replace(':', '/')}",
        enabled: params.save_intermediate_files,
        pattern: "${output_filename}.bcf*",
        mode: "copy"

    publishDir "$params.log_output_dir/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
        path samples
        path input_bcf
        path input_bcf_csi
        val tumor_id

    output:
        path "${output_filename}.bcf", emit: somatic_bcf
        path "${output_filename}.bcf.csi", emit: somatic_bcf_csi
        path ".command.*"
 
    script:
        output_filename = generate_standard_filename(
            "DELLY-${params.delly_version}",
            params.dataset_id,
            tumor_id,
            [additional_information: "somatic-filtered"]
            )
        """
        set -euo pipefail
        delly filter -f somatic -s ${samples} -o ${output_filename}.bcf "$input_bcf"
        """
    }
