#!/usr/bin/env nextflow

def docker_image_bcftools = "blcdsdockerregistry/call-gsv:bcftools-${params.bcftools_version}"

log.info """\
------------------------------------
          B C F T O O L S
------------------------------------
Docker Images:
- docker_image_bcftools:   ${docker_image_bcftools}
"""

process bcftools_merge_bcfs {
    container docker_image_bcftools

    publishDir params.output_dir,
        enabled: params.save_intermediate_files,
        pattern: "merged_${mode}.bcf*",
        mode: "copy",
        saveAs: { "bcftools_merge_${mode}_bcfs-${params.bcftools_version}/${file(it).getName()}" }

    publishDir params.output_log_dir,
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "bcftools_merge_${mode}_bcfs/log${file(it).getName()}" }

    input:
        path bcfs
        val bcfs_list
        val mode
        val force_merge_duplicated_samples

    output:
        path "merged_${mode}.bcf", emit: merged
        path "merged_${mode}.bcf.csi", emit: merged_index
        path ".command.*"

    script:
        bcfs_concat_string = bcfs_list.join(' ')
        /**
        * Multiple bcfs might come from a same sample, normally this should not happen, but provide an option to allow this.
        */
        if (force_merge_duplicated_samples)
            fs = "--force-samples"
        else
            fs = ""
        """
        set -euo pipefail
        bcftools merge ${fs} -m id -O b -o merged_${mode}.bcf ${bcfs_concat_string}
        bcftools index merged_${mode}.bcf
        """
    }
