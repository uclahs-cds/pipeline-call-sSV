#!/usr/bin/env nextflow

log.info """\
------------------------------------
             D E L L Y
------------------------------------
Docker Images:
- docker_image_delly: ${params.docker_image_delly}
"""

process call_sSV_Delly{
    container params.docker_image_delly

    publishDir "$params.output_dir/output",
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
        path "DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_name}.bcf", emit: nt_call_bcf
        path "DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_name}.bcf.csi", emit: nt_call_bcf_csi
        path "${tumor_sample_name}", emit: samples
        path ".command.*"
        val tumor_sample_name, emit: tumor_sample_name

    script:
        """
        set -euo pipefail
        delly call \
            --genome "$reference_fasta" \
            --exclude "$exclusion_file" \
            --outfile "DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_name}.bcf" \
            "$tumor_sample_bam" \
            "$control_sample_bam"

        touch "${tumor_sample_name}"
        """
    }

process filter_sSV_Delly{
    container params.docker_image_delly

    publishDir "$params.output_dir/output",
        pattern: "filtered_somatic_${tumor_sample_name}.bcf*",
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
        path "filtered_somatic_${tumor_sample_name}.bcf", emit: filtered_somatic_bcf
        path "filtered_somatic_${tumor_sample_name}.bcf.csi", emit: filtered_somatic_bcf_csi
        path ".command.*"
 
    script:
        """
        set -euo pipefail
        delly filter -f somatic -s ${samples} -o filtered_somatic_${tumor_sample_name}.bcf "$input_bcf"
        """
    }
