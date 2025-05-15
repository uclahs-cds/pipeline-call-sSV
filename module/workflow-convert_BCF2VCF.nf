nextflow.enable.dsl=2

include { convert_BCF2VCF_BCFtools } from "../external/pipeline-Nextflow-module/modules/BCFtools/convert_BCF2VCF_BCFtools/main.nf" addParams(
    options: [
        docker_image: params.docker_image_bcftools,
        output_dir: "${params.workflow_output_dir}/output",
        log_output_dir: "${params.log_output_dir}/process-log/",
        ]
    )
include { compress_VCF } from "./workflow-compress_VCF.nf" addParams(
    options: [
        docker_image: params.docker_image_bcftools,
        output_dir: "${params.workflow_output_dir}",
        log_output_dir: "${params.log_output_dir}/process-log/",
        ]
    )

workflow convert_BCF2VCF {
    take:
    sample_name
    variant_file
    variant_file_index

    main:

    bcf2vcf_channel = sample_name
        .combine(variant_file)
        .combine(variant_file_index)

    convert_BCF2VCF_BCFtools(bcf2vcf_channel)

    compress_VCF(sample_name, convert_BCF2VCF_BCFtools.out.vcf)

    emit:
    gzvcf = compress_VCF.out.gzvcf
    idx = compress_VCF.out.idx
}
