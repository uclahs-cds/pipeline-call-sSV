#!/usr/bin/env nextflow

log.info """\
------------------------------------
          B C F T O O L S
------------------------------------
Docker Images:
- docker_image_bcftools: ${params.docker_image_bcftools}
"""

process query_SampleName_BCFtools {
    container params.docker_image_bcftools

    publishDir "${params.output_dir}/intermediate/${task.process.replace(':', '/')}",
        enabled: params.save_intermediate_files,
        pattern: "${tmp_samples}.tsv",
        mode: "copy"

    publishDir "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
        path input_bcf
        path tmp_samples

    output:
        path ".command.*"
        path "${tmp_samples}.tsv", emit: samples

    script:
        """
        set -euo pipefail

        echo -e "tumor\ncontrol" > samples_type

        bcftools query -l $input_bcf > ${tmp_samples}

        paste ${tmp_samples} samples_type > "${tmp_samples}.tsv"
        """
    }

process filter_BCF_BCFtools {
    container params.docker_image_bcftools

    publishDir "${params.output_dir}/output",
        pattern: "${filename_base}_nonPassCallsFiltered.bcf*",
        mode: "copy"

    publishDir "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
        path input_bcf

    output:
        path ".command.*"
        path "${filename_base}_nonPassCallsFiltered.bcf", emit: nonPassCallsFiltered
        path "${filename_base}_nonPassCallsFiltered.bcf.csi", emit: nonPassCallsFiltered_csi

    script:
        filename_base = file(input_bcf).baseName
        """
        set -euo pipefail

        bcftools view -i 'FILTER=="PASS"' -O b -o "${filename_base}_nonPassCallsFiltered.bcf" $input_bcf
        
        bcftools index "${filename_base}_nonPassCallsFiltered.bcf"
        """
    }
