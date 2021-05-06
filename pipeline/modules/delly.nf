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

    publishDir params.output_dir,
        enabled: params.save_intermediate_files,
        pattern: "*_samples.tsv",
        mode: "copy",
        saveAs: { "delly-${params.delly_version}/${file(it).getName()}" }

    input:
        tuple val(tumor_sample_name), path(tumor_sample_bam), path(tumor_sample_bai), val(control_sample_name), path(control_sample_bam), path(control_sample_bai)
        path reference_fasta
        path reference_fasta_fai
        path exclusion_file
        //path control_sample_bams
        //val control_sample_bams_list
        val mode

    output:
        //path "DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_name}_${control_sample_name}.bcf"
        //path "DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_name}_${control_sample_name}.bcf.csi"
        path ".command.*"
        path "*_samples.tsv"

    script:
        //bams_concat_string = bams_list.join(' ')
        if (params.SINGLE_CTRL_SAMPLE == mode)
            """
            set -euo pipefail
            #delly call \
            #    --genome $reference_fasta \
            #    --exclude $exclusion_file \
            #    --outfile DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_name}_${control_sample_name}.bcf \
            #    $tumor_sample_bam \
            #    $control_sample_bam

            echo -e "${control_sample_name}\tcontrol" > ${tumor_sample_name}_samples.tsv
            echo -e "${tumor_sample_name}\ttumor" >> ${tumor_sample_name}_samples.tsv
            """
        else if (params.MULTIPLE_CTRL_SAMPLES == mode)
            """
            set -euo pipefail
            delly call \
                -v $somatic_sites \
                --genome $reference_fasta \
                --exclude $exclusion_file \
                --outfile DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_name}_${control_sample_name}.bcf \
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
        pattern: "${mode}_filtered_somatic.bcf*",
        mode: "copy",
        saveAs: { "delly-${params.delly_version}/${file(it).getName()}" }

    publishDir params.output_log_dir,
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "delly_filter_sv/log${file(it).getName()}" }

    input:
        path samples
        path input_bcf
        val mode

    output:
        path "filtered_germline_${mode}.bcf", emit: filtered_germline_bcf
        path "filtered_germline_${mode}.bcf.csi"
        path ".command.*"
 
    script:    
        """
        set -euo pipefail
        delly filter -f somatic -s ${samples} -o ${mode}_filtered_somatic.bcf $input_bcf
        """
    }
