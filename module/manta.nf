#!/usr/bin/env nextflow

log.info """\
------------------------------------
             M A N T A
------------------------------------
Docker Images:
- docker_image_manta: ${params.docker_image_manta}
"""

include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

process call_sSV_Manta {
    container params.docker_image_manta

    publishDir "$params.output_dir/${params.docker_image_manta.split("/")[1].replace(':', '-').toUpperCase()}/output",
        pattern: "*vcf.gz*",
        mode: "copy"

    publishDir "$params.output_dir/${params.docker_image_manta.split("/")[1].replace(':', '-').toUpperCase()}/QC",
        pattern: "*Stats*",
        mode: "copy"

    publishDir "$params.log_output_dir/process-log",
        pattern: ".command.*",
        mode: "copy",
        saveAs: { "${task.process.replace(':', '/')}/log${file(it).getName()}" }

    input:
        tuple(val(tumor_id), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai))
        path reference_fasta
        path reference_fasta_fai

    output:
        path "${output_filename}_candidateSmallIndels.vcf.gz", emit: vcf_small_indel_sv_file
        path "${output_filename}_candidateSmallIndels.vcf.gz.tbi", emit: vcf_small_indel_sv_tbi
        path "${output_filename}_diploidSV.vcf.gz", emit: vcf_diploid_sv_file
        path "${output_filename}_diploidSV.vcf.gz.tbi", emit: vcf_diploid_sv_tbi
        path "${output_filename}_candidateSV.vcf.gz", emit: vcf_candidate_sv_file
        path "${output_filename}_candidateSV.vcf.gz.tbi", emit: vcf_candidate_sv_tbi
        path "${output_filename}_somaticSV.vcf.gz", emit: vcf_somatic_sv_file
        path "${output_filename}_somaticSV.vcf.gz.tbi", emit: vcf_somatic_sv_tbi
        path "${output_filename}_*Stats*"
        path ".command.*"
        val tumor_id, emit: tumor_id
    
    script:
        output_filename = generate_standard_filename(
            "Manta-${params.manta_version}",
            params.dataset_id,
            tumor_id,
            [:]
            )
        """
        set -euo pipefail
        configManta.py \
            --normalBam "$normal_bam" \
            --tumorBam "$tumor_bam" \
            --referenceFasta "$reference_fasta" \
            --runDir MantaWorkflow
        MantaWorkflow/runWorkflow.py

        # re-name Manta outputs based on output file name standardization - `output_filename`
        for variant_file in `ls MantaWorkflow/results/variants/*`
            do
                variant_file_base_name=`basename \${variant_file}`
                mv \${variant_file} ./${output_filename}_\${variant_file_base_name}
            done
        
        for stats_file in `ls MantaWorkflow/results/stats/*`
            do
                stats_file_base_name=`basename \${stats_file}`
                mv \${stats_file} ./${output_filename}_\${stats_file_base_name}
            done
        
        """
    }
