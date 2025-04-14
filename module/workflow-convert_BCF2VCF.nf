nextflow.enable.dsl=2

include { convert_BCF2VCF_BCFtools } from "../external/pipeline-Nextflow-module/modules/BCFtools/convert_BCF2VCF_BCFtools/main.nf" addParams(
    options: [
        docker_image: params.docker_image_bcftools,
        main_process: "DELLY-${params.delly_version}",
        output_dir: "${params.workflow_output_dir}/output",
        log_output_dir: "${params.log_output_dir}",
        ]
    )
include { compress_index_VCF } from "../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf" addParams(
    options: [
        docker_image: params.docker_image_bcftools,
        output_dir: "${params.workflow_output_dir}/output",
        log_output_dir: "${params.log_output_dir}/process-log/DELLY-${params.delly_version}",
        ]
    )

workflow convert_BCF2VCF {
    take:
    sample_name
    bcf
    index

    main:

    bcf2vcf_channel = sample_name
        .combine(bcf)
        .combine(index)

    bcf2vcf_channel.view()

    convert_BCF2VCF_BCFtools(bcf2vcf_channel)

    index_channel = sample_name.combine(convert_BCF2VCF_BCFtools.out.vcf)
    index_channel.view()

    compress_index_VCF(index_channel)

    emit:
    gzvcf = compress_index_VCF.out.index_out.map{ it -> ["${it[1]}"] }
    idx = compress_index_VCF.out.index_out.map{ it -> ["${it[2]}"] }
}
