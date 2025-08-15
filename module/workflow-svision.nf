nextflow.enable.dsl=2

include { call_SV_SVision } from "./svision.nf"
include { compress_index_VCF } from "../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf" addParams(
    options: [
        docker_image: params.docker_image_bcftools,
        output_dir: "${params.workflow_output_dir}",
        log_output_dir: "${params.log_output_dir}/process-log/",
        ]
    )
include { generate_sha512 as generate_sha512_SVision } from "./sha512.nf"

workflow workflow_SVision {
    take:
    bams_to_call

    main:
    call_SV_SVision(
        bams_to_call,
        params.svision_cnn_model,
        params.reference_fasta,
        "${params.reference_fasta}.fai"
    )

    compress_index_VCF(call_SV_SVision.out.vcf)

    compress_index_VCF.out.index_out.map{ it -> ["${it[1]}"] }
        .mix(compress_index_VCF.out.index_out.map{ it -> ["${it[2]}"] })
        .set{ files_for_checksum }

    generate_sha512_SVision(files_for_checksum)
}
