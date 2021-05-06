#!/usr/bin/env nextflow

def docker_image_delly = "blcdsdockerregistry/call-gsv:delly-${params.delly_version}"

log.info """\
------------------------------------
             D E L L Y
------------------------------------
Docker Images:
- docker_image_delly:   ${docker_image_delly}
"""

process delly_call_NT {
    container docker_image_delly

    publishDir params.output_dir,
        enabled: params.save_intermediate_files,
        pattern: "*.bcf*",
        mode: "copy",
        saveAs: { "delly-${params.delly_version}/${file(it).getName()}" }

    publishDir params.output_log_dir,
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "delly_call_${mode}/log${file(it).getName()}" }

    input:
        tuple val(tumor_sample_name), path(tumor_sample_bam), path(tumor_sample_bai), val(control_sample_name), path(control_sample_bam), path(control_sample_bai)
        path reference_fasta
        path reference_fasta_fai
        path exclusion_file
        //path control_sample_bams
        //val control_sample_bams_list
        val mode

    output:
        path "DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_name}_${mode}.bcf"
        path "DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_name}_${mode}.bcf.csi"
        path ".command.*"

    script:
        //bams_concat_string = bams_list.join(' ')
        if (params.SINGLE_CTRL_SAMPLE == mode)
            """
            set -euo pipefail
            delly call \
                --genome $reference_fasta \
                --exclude $exclusion_file \
                --outfile DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_name}_${mode}.bcf \
                $tumor_sample_bam \
                $control_sample_bam
            """
        else if (params.MULTIPLE_CTRL_SAMPLES == mode)
            """
            set -euo pipefail
            delly call \
                -v $somatic_sites \
                --genome $reference_fasta \
                --exclude $exclusion_file \
                --outfile DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_name}_${mode}.bcf \
                $tumor_sample_bam \
                $control_sample_bams
            """
        else
            error "Invalid mode in delly_call_sv: ${mode}"
    }

process delly_filter_sv {
    container docker_image_delly

    publishDir params.output_dir,
        enabled: params.save_intermediate_files,
        pattern: "filtered_germline_${mode}.bcf*",
        mode: "copy",
        saveAs: { "delly-${params.delly_version}/${file(it).getName()}" }

    publishDir params.output_log_dir,
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "delly_filter/log${file(it).getName()}" }

    input:
        path merged
        path merged_index
        val mode

    output:
        path "filtered_germline_${mode}.bcf", emit: filtered_germline_bcf
        path "filtered_germline_${mode}.bcf.csi"
        path ".command.*"
 
    script:    
        if (params.GSV == mode)
            """
            set -euo pipefail
            delly filter -f germline -o filtered_germline_${mode}.bcf "$merged"
            """
        else if (params.GCNV == mode)
            """
            set -euo pipefail
            delly classify -f germline -o filtered_germline_${mode}.bcf "$merged"
            """
        else
            error "Invalid mode in delly_filter: ${mode}"    
    }
