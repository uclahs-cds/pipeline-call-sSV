#!/usr/bin/env nextflow

def docker_image_bcftools = "blcdsdockerregistry/bcftools:${params.bcftools_version}"

log.info """\
------------------------------------
          B C F T O O L S
------------------------------------
Docker Images:
- docker_image_bcftools:   ${docker_image_bcftools}
"""

process query_sample_name_Bcftools {
    container docker_image_bcftools

    publishDir params.output_dir,
        enabled: params.save_intermediate_files,
        pattern: "${tmp_samples}.tsv",
        mode: "copy",
        saveAs: { "bcftools-${params.bcftools_version}/${file(it).getName()}" }

    publishDir params.output_log_dir,
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "query_sample_name_Bcftools/log${file(it).getName()}" }

    input:
        path input_bcf
        path tmp_samples
        path sample_types

    output:
        path ".command.*"
        path "${tmp_samples}.tsv", emit: samples

    script:
        """
        set -euo pipefail

        bcftools query -l $input_bcf > ${tmp_samples}

        paste ${tmp_samples} ${sample_types} > "${tmp_samples}.tsv"
        """
    }
