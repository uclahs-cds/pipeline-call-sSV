nextflow.enable.dsl=2

include { compress_index_VCF } from "../external/pipeline-Nextflow-module/modules/common/index_VCF_tabix/main.nf" addParams(
    options: [
        docker_image: params.docker_image_bcftools,
        output_dir: "${params.workflow_output_dir}",
        log_output_dir: "${params.log_output_dir}/process-log/",
        ]
    )

workflow compress_VCF {
    take:
    sample_name
    variant_file

    main:

    index_channel = sample_name.combine(variant_file)
    compress_index_VCF(index_channel)

    emit:
    gzvcf = compress_index_VCF.out.index_out.map{ it -> ["${it[1]}"] }
    idx = compress_index_VCF.out.index_out.map{ it -> ["${it[2]}"] }
    }
