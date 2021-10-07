#!/usr/bin/env nextflow

log.info """\
------------------------------------
         V A L I D A T I O N
------------------------------------
Docker Images:
- docker_image_validate: ${params.docker_image_validate}
"""

process validate_file {
    container params.docker_image_validate

    input:
        path(file_to_validate)

    script:
        """
        set -euo pipefail
        python -m validate -t file-input ${file_to_validate}
        """
    }
