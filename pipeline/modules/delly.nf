#!/usr/bin/env nextflow

def docker_image_delly = "blcdsdockerregistry/delly:${params.delly_version}"

log.info """\
------------------------------------
             D E L L Y
------------------------------------
Docker Images:
- docker_image_delly:   ${docker_image_delly}
"""

process call_sSV_Delly{
    container docker_image_delly

    publishDir params.output_dir,
        enabled: params.save_intermediate_files,
        pattern: "*.bcf*",
        mode: "copy",
        saveAs: { "delly-${params.delly_version}/${file(it).getName()}" }

    publishDir params.output_log_dir,
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "call_sSV_Delly/log${file(it).getName()}" }

    input:
        tuple(val(tumor_sample_bam_name), path(tumor_sample_bam), path(tumor_sample_bai), val(control_sample_bam_name), path(control_sample_bam), path(control_sample_bai))
        path reference_fasta
        path reference_fasta_fai
        path exclusion_file

    output:
        path "DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_bam_name}.bcf", emit: nt_call_bcf
        path "DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_bam_name}.bcf.csi", emit: nt_call_bcf_csi
        path "${tumor_sample_bam_name}_samples", emit: samples
        path ".command.*"

    script:
        """
        set -euo pipefail
        delly call \
            --genome "$reference_fasta" \
            --exclude "$exclusion_file" \
            --outfile "DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_bam_name}.bcf" \
            "$tumor_sample_bam" \
            "$control_sample_bam"

        touch "${tumor_sample_bam_name}_samples"
        """
    }

process filter_sSV_Delly{
    container docker_image_delly

    publishDir params.output_dir,
        enabled: params.save_intermediate_files,
        pattern: "filtered_somatic.bcf*",
        mode: "copy",
        saveAs: { "delly-${params.delly_version}/${file(it).getName()}" }

    publishDir params.output_log_dir,
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "filter_sSV_Delly/log${file(it).getName()}" }

    input:
        path samples
        path input_bcf
        path input_bcf_csi

    output:
        path "filtered_somatic.bcf", emit: filtered_somatic_bcf
        path "filtered_somatic.bcf.csi", emit: filtered_somatic_bcf_csi
        path ".command.*"
 
    script:
        """
        set -euo pipefail
        delly filter -f somatic -s ${samples} -o filtered_somatic.bcf "$input_bcf"
        """
    }
