# Call-sSV Pipeline

* [Overview](#Overview)
* [How To Run](#How-To-Run)
* [Flow Diagram](#flow-diagram)
* [Pipeline Steps](#pipeline-steps)
* [Inputs](#Inputs)
* [Outputs](#outputs)
* [Testing and Validation](#testing-and-validation)
    * [Test Data Set](#test-data-set)
    * [Performance Validation](#performance-validation)
    * [Validation Tool](#validation-tool)
* [References](#references)
* [License](#license)

## Overview:
The call-sSV pipeline calls somatic structural variants utilizing [DELLY](https://github.com/dellytools/delly) and [Manta](https://github.com/Illumina/manta). This pipeline requires at least one tumor sample and a matched normal sample.
This pipeline is developed using Nextflow, docker and can run either on a single node linux machine or a multi-node HPC Slurm cluster.

## How to Run:

Below is a summary of how to run the pipeline.  See [here](https://confluence.mednet.ucla.edu/pages/viewpage.action?spaceKey=BOUTROSLAB&title=How+to+run+a+nextflow+pipeline) for full instructions.

Pipelines should be run **WITH A SINGLE SAMPLE AT A TIME**. Otherwise resource allocation and Nextflow errors could cause the pipeline to fail.

1. The recommended way of running the pipeline is to directly use the source code located here: `/hot/software/pipeline/pipeline-call-sSV/Nextflow/release/`, rather than cloning a copy of the pipeline.

    * The source code should never be modified when running our pipelines

2. Create a config file for input, output, and parameters. An example for a config file can be found [here](config/template.config). See [Nextflow Config File Parameters](#Nextflow-Config-File-Parameters) for the detailed description of each variable in the config file.

    * Do not directly modify the source `template.config`, but rather you should copy it from the pipeline release folder to your project-specific folder and modify it there

3. Create the input YAML using the [template](input/call-sSV-input.yaml).See [Input YAML](#Input-YAML) for detailed description of each column.

   * Again, do not directly modify the source template input YAML file.  Instead, copy it from the pipeline release folder to your project-specific folder and modify it there.

4. The pipeline can be executed locally using the command below:

- YAML input
```bash
nextflow run path/to/main.nf -config path/to/sample-specific.config -params-file path/to/input.yaml
```

* For example, `path/to/main.nf` could be: `/hot/software/pipeline/pipeline-call-sSV/Nextflow/release/6.0.0/main.nf`
* `path/to/sample-specific.config` is the path to where you saved your project-specific copy of [template.config](config/template.config)
* `path/to/input.yaml` is the path to where you saved your sample-specific copy of [input-sSV.yaml](input/call-sSV-input.yaml)

To submit to UCLAHS-CDS's Azure cloud, use the submission script [here](https://github.com/uclahs-cds/tool-submit-nf) with the command below:

```bash
python path/to/submit_nextflow_pipeline.py \
    --nextflow_script path/to/main.nf \
    --nextflow_config path/to/sample-specific.config \
    --nextflow_yaml path/to/input.yaml \
    --pipeline_run_name <sample_name> \
    --partition_type F16 \
    --email <your UCLA email, jdoe@ucla.edu>
```
In the above command, the partition type can be changed based on the size of the dataset. At this point, node F16 is generally recommended for larger datasets like A-full and node F2 for smaller datasets like A-mini.

\* Manta SV calling wouldn't work on an F2 node due to incompatible resources. In order to test the pipeline for tasks not relevant to Manta, please set `algorithm = ['delly']` in the sample specific [config](config/template.config) file.

> **Note**: Because this pipeline uses an image stored in the GitHub Container Registry, you must follow the steps listed in the [Docker Introduction](https://confluence.mednet.ucla.edu/display/BOUTROSLAB/Docker+Introduction#DockerIntroduction-GitHubContainerRegistryGitHubContainerRegistry|Setup) on Confluence to set up a PAT for your GitHub account and log into the registry on the cluster before running this pipeline.

## Flow Diagram:

![call-sSV flow diagram](call-sSV-workflow.svg?raw=true)

## Pipeline Steps:

### Call Somatic Structural Variants - DELLY workflow:

#### 1. Calling Single Sample Somatic Structural Variants
```script
delly call --genome hg38.fa --exclude hg38.excl --map-qual 20 --min-clique-size 5 --mad-cutoff 15 --outfile t1.bcf tumor1.bam normal1.bam
```
This step requires an aligned and sorted tumor sample BAM file and a matched normal sample as an input for variant calling with DELLY.
The stringent filters (`--map-qual 20` `--min-clique-size 5` `--mad-cutoff 15`) are added, which can drastically reduce the runtime, especially when the input BAMs are big. In the pipeline, these filters are specified in the NextFlow input parameters [config file](config/template.config). If need be, these stringent filters can be adjusted in the config file.

#### 2. Query the generated bcfs to get the sample names, which will be used in step 3.
```script
echo -e "tumor\ncontrol" > samples_type
bcftools query -l t1.bcf > samples_name
paste samples_name samples_type > samples.tsv
```

#### 3. Somatic Filtering
```script
delly filter -f somatic -o t1.pre.bcf -s samples.tsv t1.bcf
```
This step applies somatic filtering against the `.bcf` file generated in Step 1.

Note: cohort based false positive filtering is compuationally heavy and not implemented in this pipeline.

### Call Somatic Structural Variants - Manta workflow:

#### 1. Calling Single Sample Somatic Structural Variants
```script
configManta.py --normalBam "${normal_bam}" --tumorBam "${tumor_bam}" --referenceFasta "${reference_fasta}" --runDir MantaWorkflow
MantaWorkflow/runWorkflow.py
```
This step requires an aligned and sorted tumor sample BAM file and a matched normal sample as an input for variant calling with Manta.

## Inputs

### Input YAML

| Field | Type | Description |
|:------|:-----|:------------|
| sample_id | string | Tumor ID |
| normal | path | Set to absolute path to normal BAM |
| tumor | path | Set to absolute path to tumour BAM |

```
---
sample_id: "sample_id"
input:
  BAM:
    normal:
      - "/absolute/path/to/BAM"
      - "/absolute/path/to/BAM"
    tumor:
      - "/absolute/path/to/BAM"
      - "/absolute/path/to/BAM"

```

## Nextflow Config File Parameters

| Field |	Required |	Type |	Description |
| ------- |   --------- | ------ | -------------|
| dataset_id |	yes	| string |	Boutros Lab dataset id |
| blcds_registered_dataset	| yes |	boolean | Affirms if dataset should be registered in the Boutros Lab Data registry. Default value is `false`. |
| algorithm | yes | list | List containing a combination of SV callers `delly`, `manta`. List can contain a single caller of choice.  |
| reference_fasta	| yes |	path	| Absolute path to the reference genome FASTA file. The reference genome is used by Delly for structural variant calling. GRCh37 - /hot/ref/reference/GRCh37-EBI-hs37d5/hs37d5.fa, GRCh38 - /hot/ref/reference/GRCh38-BI-20160721/Homo_sapiens_assembly38.fasta |
| exclusion_file |	yes	| path | Absolute path to the Delly reference genome exclusion file utilized to remove suggested regions for structural variant calling. GRCh37 - /hot/ref/tool-specific-input/Delly/GRCh37-EBI-hs37d/human.hs37d5.excl.tsv, GRCh38 - /hot/ref/tool-specific-input/Delly/hg38/human.hg38.excl.tsv |
| map_qual | yes | integer | Minimum paired-end (PE) mapping quality (MAPQ) for Delly. Default set to 20.|
| min_clique_size | yes | integer | Minimum number of supporting PE or split-read (SR) alignments required for a clique to be identified as a structural variant by Delly. Adjust this parameter to control the sensitivity and specificity of Delly variant calling. Default set to 5.|
| mad_cutoff | yes | integer | Insert size cutoff, median+s*MAD (deletions only) for Delly. Default set to 15.|
| save_intermediate_files |	yes	| boolean |	Optional parameter to indicate whether intermediate files will be saved. Default value is `false`. |
| output_dir |	yes |	path |	Absolute path to the directory where the output files to be saved. |
| work_dir	| no	| path |	Path of working directory for Nextflow. When included in the sample config file, Nextflow intermediate files and logs will be saved to this directory. With `ucla_cds`, the default is `/scratch` and should only be changed for testing/development. Changing this directory to `/hot` or `/tmp` can lead to high server latency and potential disk space limitations, respectively. |
| verbose |	yes |	boolean	| If set to `true`, the values of input channels will be printed, can be used for debugging|
| `docker_container_registry` | optional | string | Registry containing tool Docker images. Default: `ghcr.io/uclahs-cds` |

An example of the NextFlow Input Parameters Config file can be found [here](config/template.config).

## Outputs

| Output |	Description |
| ---- | -------- |
| .bcf | Binary VCF output format from DELLY with somatic structural variants if found. |
| .bcf.csi | CSI-format index for BCF files from DELLY. |
| .vcf.gz | zipped VCF output format from Manta with somatic structural variants if found. |
| .vcf.gz.tbi | TBI-format index for zipped VCF files from Manta. |
| report.html, timeline.html and trace.txt | A Nextflow report, timeline and trace files. |
| \*.log.command.* | Process and sample specific logging files created by nextflow. |
| *.sha512 | Generates SHA-512 hash to validate file integrity. |


## Testing and Validation

### Test Data Set

| Data Set | Run Configuration | YAML input | Output Dir | Normal Sample | Tumor Sample |
| ------ | ------ | ------- | ------ | ------- | ------- |
| A-mini TWGSAMIN000001-T003-S03-F | /hot/software/pipeline/pipeline-call-sSV/Nextflow/development/unreleased/mmootor-replace-input-csv-with-yaml/TWGSAMIN000001-T003-S03-F.config | /hot/software/pipeline/pipeline-call-sSV/Nextflow/development/input/yaml/call-sSV-input-TWGSAMIN000001-T003-S03-F.yaml | /hot/software/pipeline/pipeline-call-sSV/Nextflow/development/unreleased/mmootor-replace-input-csv-with-yaml/TWGSAMIN000001-T003-S03-F/  | /hot/software/pipeline/pipeline-call-sSV/Nextflow/development/input/data/TWGSAMIN000001-N003-S03-F.bam | /hot/software/pipeline/pipeline-call-sSV/Nextflow/development/input/data/TWGSAMIN000001-T003-S03-F.bam |
| ILHNLNEV000009 | /hot/software/pipeline/pipeline-call-sSV/Nextflow/development/unreleased/mmootor-replace-input-csv-with-yaml/ILHNLNEV000009-T002-L01-F_F32.config | /hot/software/pipeline/pipeline-call-sSV/Nextflow/development/unreleased/mmootor-replace-input-csv-with-yaml/ILHNLNEV000009-T002-L01-F_F32.yaml | /hot/software/pipeline/pipeline-call-sSV/Nextflow/development/unreleased/mmootor-replace-input-csv-with-yaml/ILHNLNEV000009-T002-L01-F/ | /hot/project/disease/HeadNeckTumor/HNSC-000084-LNMEvolution/pipelines/call-gSNP/2020-12-22/ILHNLNEV000009-T002-L01-F//gSNP/2021-01-22_11.01.06/ILHNLNEV000009/SAMtools-1.10_Picard-2.23.3/recalibrated_reheadered_bam_and_bai/ILHNLNEV000009-N001-B01-F_realigned_recalibrated_reheadered.bam | /hot/project/disease/HeadNeckTumor/HNSC-000084-LNMEvolution/pipelines/call-gSNP/2020-12-22/ILHNLNEV000009-T002-L01-F//gSNP/2021-01-22_11.01.06/ILHNLNEV000009/SAMtools-1.10_Picard-2.23.3/recalibrated_reheadered_bam_and_bai/ILHNLNEV000009-T002-L01-F_realigned_recalibrated_reheadered.bam |
| DTB-266_WCDT | /hot/software/pipeline/pipeline-call-sSV/Nextflow/development/5.0.0/mmootor-release-5-0-0-rc-1/config/DTB-266_WCDT_F72.config | - |  /hot/software/pipeline/pipeline-call-sSV/Nextflow/development/unreleased/mmootor-release-5-0-0-rc-1/DTB-266_WCDT/ | /hot/data/unregistered/Quigley-Gebo-PRAD-SVMW/processed/output_call-gSNP/call-gSNP-DSL2-0.0.1/DTB-266/GATK-4.1.9.0/output/DTB-266_DNA_N_realigned_recalibrated_merged.bam | /hot/data/unregistered/Quigley-Gebo-PRAD-SVMW/processed/output_call-gSNP/call-gSNP-DSL2-0.0.1/DTB-266/GATK-4.1.9.0/output/DTB-266_DNA_T_realigned_recalibrated_merged.bam |

## Performance Validation

Testing was performed primarily in the Boutros Lab SLURM Development cluster. Metrics below will be updated where relevant with additional testing and tuning outputs.

| Test Case | Test Date | Node Type | Duration | CPU Hours | Peak RSS (RAM) |
|----- | -------| --------| ----------| ---------| --------|
| TWGSAMIN000001-T003-S03-F | 2023-01-19 | F16 | 41m 35s | 0.7 | 1.8 GB |
| ILHNLNEV000009-T002-L01-F | 2023-01-20 | F32 | 1d 23h 10m 46s | 63.3 | 12.1 GB |
| DTB-266_WCDT | 2023-01-19 | F72 | 22h 55m 17s | 45.1 | 13.2 GB |
|ILHNLNEV000005-T002-L01-F (with stringent filters. See [#10](https://github.com/uclahs-cds/pipeline-call-sSV/issues/10) [2f72de1](https://github.com/uclahs-cds/pipeline-call-sSV/commit/2f72de1ba190623e4344f144a12cc315fda1dd18)) | 2021-10-02 | F72 | 1d 10h 55m 13s | 2'478.4h | 11.590 GB |


## References

* [DELLY Structural Variant Caller](https://github.com/dellytools/delly)
* [Manta Structural Variant Caller](https://github.com/Illumina/manta)


## License

Authors: Yu Pan (YuPan@mednet.ucla.edu), Ghouse Mohammed (GMohammed@mednet.ucla.edu), Mohammed Faizal Eeman Mootor (MMootor@mednet.ucla.edu)

`call-sSV` is licensed under the GNU General Public License version 2. See the file LICENSE for the terms of the GNU GPL license.

`call-sSV` takes BAM files and utilizes DELLY and Manta to call somatic structural variants.

Copyright (C) 2021-2023 University of California Los Angeles ("Boutros Lab") All rights reserved.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
