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
        saveAs: { "delly_call/log${file(it).getName()}" }

    publishDir params.output_dir,
        enabled: params.save_intermediate_files,
        pattern: "*_samples.tsv",
        mode: "copy",
        saveAs: { "delly-${params.delly_version}/${file(it).getName()}" }

    input:
        tuple(val(tumor_sample_name), path(tumor_sample_bam), path(tumor_sample_bai), val(control_sample_name), path(control_sample_bam), path(control_sample_bai))
        path reference_fasta
        path reference_fasta_fai
        path exclusion_file

    output:
        path "DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_name}_${control_sample_name}.bcf", emit: nt_call_bcf
        path "DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_name}_${control_sample_name}.bcf.csi", emit: nt_call_bcf_csi
        path ".command.*"
        path "${tumor_sample_name}_samples.tsv", emit: samples

    script:
        """
        set -euo pipefail
        delly call \
            --genome "$reference_fasta" \
            --exclude "$exclusion_file" \
            --outfile "DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_name}_${control_sample_name}.bcf" \
            "$tumor_sample_bam" \
            "$control_sample_bam"

        echo -e "${control_sample_name}\tcontrol" > "${tumor_sample_name}_samples.tsv"
        echo -e "${tumor_sample_name}\ttumor" >> "${tumor_sample_name}_samples.tsv"
        """
    }

process delly_regenotype_NT {
    container docker_image_delly

    publishDir params.output_dir,
        enabled: params.save_intermediate_files,
        pattern: "*.bcf*",
        mode: "copy",
        saveAs: { "delly-${params.delly_version}/${file(it).getName()}" }

    publishDir params.output_log_dir,
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "delly_regenotype/log${file(it).getName()}" }

    input:
        tuple(val(tumor_sample_name), path(tumor_sample_bam), path(tumor_sample_bai))
        path reference_fasta
        path reference_fasta_fai
        path exclusion_file
        path control_samples_bams_bais_list
        //val control_samples_bams_list
        path somatic_sites

    output:
        path "DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_name}_all_control_samples.bcf", emit: nt_regenotype_bcf
        path "DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_name}_all_control_samples.bcf.csi", emit: nt_regenotype_bcf_csi
        path ".command.*"

    script:
        //control_samples_bams_concat_string = control_samples_bams_list.join(' ')
        control_samples_bams_concat_string = control_samples_bams_bais_list.findAll{!it.toString().contains("bai")}.join(' ')
        log.info("control_samples_bams_concat_string: $control_samples_bams_concat_string")
        """
        set -euo pipefail
        delly call \
            -v $somatic_sites \
            --genome $reference_fasta \
            --exclude $exclusion_file \
            --outfile "DELLY-${params.delly_version}_${params.dataset_id}_${tumor_sample_name}_all_control_samples.bcf" \
            "$tumor_sample_bam" \
            "$control_samples_bams_concat_string"
        """
    }

process delly_filter_NT {
    container docker_image_delly

    publishDir params.output_dir,
        enabled: params.save_intermediate_files,
        pattern: "filtered_somatic_${tag}.bcf*",
        mode: "copy",
        saveAs: { "delly-${params.delly_version}/${file(it).getName()}" }

    publishDir params.output_log_dir,
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "delly_filter_NT_${tag}/log${file(it).getName()}" }

    input:
        path samples
        path input_bcf
        path input_bcf_csi
        val tag

    output:
        path "filtered_somatic_${tag}.bcf", emit: filtered_somatic_bcf
        path "filtered_somatic_${tag}.bcf.csi"
        path ".command.*"
 
    script:    
        """
        set -euo pipefail
        delly filter -f somatic -s ${samples} -o filtered_somatic_${tag}.bcf "$input_bcf"
        """
    }
