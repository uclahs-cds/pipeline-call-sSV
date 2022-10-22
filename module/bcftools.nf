#!/usr/bin/env nextflow

log.info """\
------------------------------------
          B C F T O O L S
------------------------------------
Docker Images:
- docker_image_bcftools: ${params.docker_image_bcftools}
"""

include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

process query_SampleName_BCFtools {
    container params.docker_image_bcftools

    publishDir "$params.output_dir_base/intermediate/${task.process.replace(':', '/')}",
        enabled: params.save_intermediate_files,
        pattern: "${output_filename}.tsv",
        mode: "copy"

    publishDir "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
        path input_bcf
        path tmp_samples
        val tumor_id

    output:
        path ".command.*"
        path "${output_filename}.tsv", emit: samples

    script:
        output_filename = generate_standard_filename(
            "BCFtools-${params.bcftools_version}",
            params.dataset_id,
            tumor_id,
            [additional_information: "query-tumor-normal-name"]
            )
        """
        set -euo pipefail

        echo -e "tumor\ncontrol" > samples_type

        bcftools query -l $input_bcf > ${tmp_samples}

        paste ${tmp_samples} samples_type > "${output_filename}.tsv"
        """
    }

process filter_BCF_BCFtools {
    container params.docker_image_bcftools

    publishDir "$params.output_dir_base/output",
        pattern: "${output_filename}.bcf*",
        mode: "copy"

    publishDir "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
        path input_bcf
        val filter_condition
        val tumor_id

    output:
        path ".command.*"
        path "${output_filename}.bcf*", emit: nonPassCallsFiltered_and_csi

    script:
        output_filename = generate_standard_filename(
            "DELLY-${params.delly_version}",
            params.dataset_id,
            tumor_id,
            [:]
            )
        """
        set -euo pipefail

        bcftools view -i "$filter_condition" -O b -o "${output_filename}.bcf" $input_bcf
        
        bcftools index "${output_filename}.bcf"
        """
    }
