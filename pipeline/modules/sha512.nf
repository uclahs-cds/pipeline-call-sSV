#!/usr/bin/env nextflow

log.info """\
------------------------------------
          S H A - 5 1 2
------------------------------------
Docker Images:
- docker_image_validate:   ${params.docker_image_validate}
"""

process generate_sha512 {
    container params.docker_image_validate

    publishDir "$params.output_dir/output",
        pattern: "*.bcf.sha512",
        mode: "copy"

    publishDir "$params.log_output_dir/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

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
