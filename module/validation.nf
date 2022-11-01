#!/usr/bin/env nextflow

log.info """\
------------------------------------
         V A L I D A T I O N
------------------------------------
Docker Images:
- docker_image_validate: ${params.docker_image_validate}
"""

process run_validate_PipeVal {
    container params.docker_image_validate

    publishDir path: "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process}/${task.process}-${task.index}/log${file(it).getName()}" }

    input:
        path(file_to_validate)

    output:
        path ".command.*"
        path "input_validation.txt"

    script:
        """
        set -euo pipefail
        python -m validate -t file-input ${file_to_validate} > 'input_validation.txt'
        """
    }
