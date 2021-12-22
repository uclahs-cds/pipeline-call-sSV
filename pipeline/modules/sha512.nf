#!/usr/bin/env nextflow

log.info """\
------------------------------------
          S H A - 5 1 2
------------------------------------
Docker Images:
- docker_image_validate:   ${params.docker_image_validate}
"""

process generate_sha512_NonPassFiltering {
    container params.docker_image_validate

    publishDir "$params.output_dir/output",
        pattern: "*.sha512",
        mode: "copy"

    publishDir "$params.log_output_dir/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process}/${task.process}-${task.index}/log${file(it).getName()}" }

    input:
        path input_checksum_file

    output:
        path "${input_checksum_file}.sha512"
        path ".command.*"

    script:
        """
        set -euo pipefail
        sha512sum $input_checksum_file > ${input_checksum_file}.sha512
        """
    }

process generate_sha512_intermediate {
    container params.docker_image_validate

    publishDir "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
        enabled: params.save_intermediate_files,
        pattern: "*.sha512",
        mode: "copy"

    publishDir "$params.log_output_dir/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process}/${task.process}-${task.index}/log${file(it).getName()}" }

    input:
        path input_checksum_file

    output:
        path "${input_checksum_file}.sha512"
        path ".command.*"

    script:
        """
        set -euo pipefail
        sha512sum $input_checksum_file > ${input_checksum_file}.sha512
        """
    }
