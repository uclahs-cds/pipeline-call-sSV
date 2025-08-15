include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

process call_SV_SVision {
    container params.docker_image_svision

    tag "${sample_id}"

    errorStrategy 'ignore'

    publishDir "${params.workflow_output_dir}/intermediate/${task.process.replace(':', '/')}",
        enabled: params.save_intermediate_files,
        pattern: "*.vcf",
        mode: "copy"

    publishDir "${params.log_output_dir}/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
        tuple val(sample_id), path(bam), path(bam_index)
        path svision_model
        path reference_fasta
        path reference_fasta_fai

    output:
        tuple val(sample_id), path("SVision*.vcf"), emit: vcf
        path ".command.*"

    script:
    standard_filename = generate_standard_filename(
        "SVision-${params.svision_version}",
        params.dataset_id,
        sample_id,
        [:]
    )
    """
    set -euo pipefail

    export TF_XLA_FLAGS=--tf_xla_cpu_global_jit

    SVision \
        -o `pwd` \
        -b `realpath ${bam}` \
        -m `realpath ${svision_model}` \
        -g `realpath ${reference_fasta}` \
        -n "${sample_id}" \
        --min_mapq ${params.svision_min_mapq}

    for i in `ls *svision*.vcf 2>/dev/null`
    do
        mv \$i ${standard_filename}.vcf
    done
    """
}
