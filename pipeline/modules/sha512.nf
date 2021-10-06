#!/usr/bin/env nextflow

def docker_image_sha512 = "blcdsdockerregistry/validate:${params.validate_version}"

log.info """\
------------------------------------
          S H A - 5 1 2
------------------------------------
Docker Images:
- docker_image_sha512:   ${docker_image_sha512}
"""

process generate_sha512 {
    container docker_image_sha512

    publishDir "$params.output_dir",
        enabled: params.save_intermediate_files,
        pattern: "*.vcf.sha512",
        mode: "copy",
        saveAs: { "bcftools-${params.bcftools_version}/output/${file(it).getName()}" }

    publishDir "$params.output_dir",
        enabled: params.save_intermediate_files,
        pattern: "*.bcf.sha512",
        mode: "copy",
        saveAs: { "delly-${params.delly_version}/output/${file(it).getName()}" }

    publishDir "$params.output_log_dir/process-log",
        enabled: params.save_intermediate_files,
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "generate_sha512/log${file(it).getName()}" }

    input:
        path input_checksum_file

    output:
        path "${input_checksum_file}.sha512"
        path ".command.*"

    script:
        """
        set -euo pipefail
        python -m validate -t sha512-gen $input_checksum_file > ${input_checksum_file}.sha512
        """
    }
