# Call-sSV Pipeline

[![GitHub release](https://img.shields.io/github/v/release/uclahs-cds/pipeline-call-sSV)](https://github.com/uclahs-cds/pipeline-call-sSV/actions/workflows/prepare-release.yaml)

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
The call-sSV pipeline calls somatic structural variants utilizing [DELLY](https://github.com/dellytools/delly), [Manta](https://github.com/Illumina/manta), [GRIDSS2](https://github.com/PapenfussLab/gridss), and [SVision](https://github.com/xjtu-omics/SVision). This pipeline requires at least one tumor sample and a matched normal sample.
This pipeline is developed using Nextflow, docker and can run either on a single node linux machine or a multi-node HPC Slurm cluster.

## How to Run:

### Requirements
Currently supported Nextflow versions: `v23.04.2`

### Run steps
Below is a summary of how to run the pipeline.  See [here](https://uclahs-cds.atlassian.net/wiki/spaces/BOUTROSLAB/pages/3197004/How+to+run+a+nextflow+pipeline) for full instructions.

Pipelines should be run **WITH A SINGLE SAMPLE AT A TIME**. Otherwise resource allocation and Nextflow errors could cause the pipeline to fail.

1. Create a config file for input, output, and parameters. An example for a config file can be found [here](config/template.config). See [Nextflow Config File Parameters](#Nextflow-Config-File-Parameters) for the detailed description of each variable in the config file.

    * Do not directly modify the source `template.config`, but rather you should copy it from the pipeline release folder to your project-specific folder and modify it there

2. Create the input YAML using the [template](input/call-sSV-input.yaml).See [Input YAML](#Input-YAML) for detailed description.

   * Again, do not directly modify the source template input YAML file.  Instead, copy it from the pipeline release folder to your project-specific folder and modify it there.

3. The pipeline can be executed locally using the command below:

```bash
nextflow run path/to/main.nf -config path/to/sample-specific.config -params-file path/to/input.yaml
```

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

> **Note**: Because this pipeline uses an image stored in the GitHub Container Registry, you must follow the steps listed in the [Docker Introduction](https://uclahs-cds.atlassian.net/wiki/spaces/BOUTROSLAB/pages/3190419/Docker+Introduction#DockerIntroduction-HowtosetupPATandlogintotheregistryHowtosetupPATandlogintotheregistry) on Confluence to set up a PAT for your GitHub account and log into the registry on the cluster before running this pipeline.

## Flow Diagram:

![call-sSV flow diagram](docs/call-sSV.svg)

## Pipeline Steps:

### Call Somatic Structural Variants - DELLY workflow:

#### 1. Calling Single Sample Somatic Structural Variants
```script
delly call --genome hg38.fa --exclude hg38.excl --map-qual 20 --min-clique-size 5 --mad-cutoff 15 --outfile t1.bcf tumor1.bam normal1.bam
```
This step performs somatic variant calling using DELLY.
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
This step performs somatic variant calling using Manta.

### Call Somatic Structural Variants - GRIDSS2 workflow:

#### 1. Preprocess input BAMs
```script
gridss -r gridss2_reference.fasta -j gridss.jar -s preprocess input.bam
```
This step preprocesses input tumor/normal BAMs.

#### 2. Breakend Assembly
```script
gridss -r gridss2_reference.fasta -j gridss.jar -s assemble -b gridss2_blacklist.bed -a tumor_breakend_assembly.bam normal.bam tumor.bam
```
This step performs GRIDSS2 breakend assembly to produce a `tumor_breakend_assembly.bam`.
The preprocessed input tumor/normal BAMs need to be in the same directory while running breakend assembly, although the preprocced BAMs are not explicitly defined as input arguments here.

#### 3. Calling Single Sample Somatic Structural Variants
```script
gridss -r gridss2_reference.fasta -j gridss2.jar -s call -b gridss_blacklist.bed -a tumor_breakend_assembly.bam --output gridss2_tumor.vcf normal.bam tumor.bam
```
This step performs GRIDSS2 somatic variant calling.

#### 4. Somatic Filtering
```script
gridss_somatic_filter --pondir gridss2_pon_dir/ --scriptdir gridss_script_dir/ --input gridss2_tumor.vcf --output gridss2_high-confidence-somatic.vcf --fulloutput gridss2_high-low-confidence-somatic.vcf --normalordinal 1 --tumourordinal 2
```
This step performs somatic filtering to retain only somatic SVs in the filtered outputs, `gridss2_high-confidence-somatic.vcf` and `gridss2_high-low-confidence-somatic.vcf`.

### Call Structural Variants - SVision workflow:

#### 1. Call structural variants using SVision
This step runs the SVision workflow to identify structural variants in a sample.

## Inputs

### Input YAML

| Field | Type | Description |
|:------|:-----|:------------|
| sample_id | string | Tumor ID |
| normal | path | Set to absolute path to normal BAM |
| tumor | path | Set to absolute path to tumor BAM |

```
---
sample_id: "sample_id"
input:
  BAM:
    normal:
      - "/path/to/BAM"
    tumor:
      - "/path/to/BAM"

```

## Nextflow Config File Parameters

| Field |	Required |	Type |	Description |
| ------- |   --------- | ------ | -------------|
| dataset_id |	yes	| string |	Boutros Lab dataset id |
| blcds_registered_dataset	| yes |	boolean | Affirms if dataset should be registered in the Boutros Lab Data registry. Default value is `false`. |
| `genome_build` | no | string | Genome build for circos plot, `hg19` or `hg38`. Default is set to `hg38` |
| algorithm | yes | list | List containing a combination of SV callers `delly`, `manta`, `gridss2`. List can contain a single caller of choice.  |
| reference_fasta	| yes |	path	| Absolute path to the reference genome FASTA file. The reference genome is used by Delly for structural variant calling. |
| save_intermediate_files |	yes	| boolean |	Optional parameter to indicate whether intermediate files will be saved. Default value is `false`. |
| output_dir |	yes |	path |	Absolute path to the directory where the output files to be saved. |
| work_dir	| no	| path | The path to a temporary working directory for Nextflow, storing intermediate files and logs. It is recommended to use fast, local storage with high I/O performance. |
| verbose |	yes |	boolean	| If set to `true`, the values of input channels will be printed, can be used for debugging|
| `docker_container_registry` | optional | string | Registry containing tool Docker images. Default: `ghcr.io/uclahs-cds` |

An example of the NextFlow Input Parameters Config file can be found [here](config/template.config).

### DELLY Specific Parameters
| Field |	Required |	Type |	Description |
| ------- |   --------- | ------ | -------------|
| exclusion_file |	yes	| path | Absolute path to the Delly reference genome exclusion file utilized to remove suggested regions for structural variant calling. |
| map_qual | yes | integer | Minimum paired-end (PE) mapping quality (MAPQ) for Delly. Default set to 20.|
| min_clique_size | yes | integer | Minimum number of supporting PE or split-read (SR) alignments required for a clique to be identified as a structural variant by Delly. Adjust this parameter to control the sensitivity and specificity of Delly variant calling. Default set to 5.|
| mad_cutoff | yes | integer | Insert size cutoff, median+s*MAD (deletions only) for Delly. Default set to 15.|

### GRIDSS2 Specific Parameters
| Field |	Required |	Type |	Description |
| ------- |   --------- | ------ | -------------|
| gridss2_blacklist | yes | path | Path to GRIDSS2 blacklist BED file. |
| gridss2_reference_fasta | yes | path | Path to GRIDSS2 reference FASTA file. |
| gridss2_pon_dir | yes | path | Path to GRIDSS2 Panel Of Normals (PON) directory. |
| other_jvm_heap | no | string | Update `other_jvm_heap` if GRIDSS2 errors OutOfMemory. Default is `4.GB` |

### SVision Specific Parameters
| Field |	Required |	Type |	Description |
| ------- |   --------- | ------ | -------------|
| svision_cnn_model | yes | path | Path to trained CNN SVision model. |
| svision_min_mapq | yes | integer | Minimum mapping quality to consider for SVision. Default is 10.|

### Base resource allocation updaters
To optionally update the base resource (cpus or memory) allocations for processes, use the following structure and add the necessary parts to the [input.config](config/template.config) file. The default allocations can be found in the [node-specific config files](./config/)

```Nextflow
base_resource_update {
    memory = [
        [['process_name', 'process_name2'], <multiplier for resource>],
        [['process_name3', 'process_name4'], <different multiplier for resource>]
    ]
    cpus = [
        [['process_name', 'process_name2'], <multiplier for resource>],
        [['process_name3', 'process_name4'], <different multiplier for resource>]
    ]
}
```
> **Note** Resource updates will be applied in the order they're provided so if a process is included twice in the memory list, it will be updated twice in the order it's given.
Examples:

- To double memory of all processes:
```Nextflow
base_resource_update {
    memory = [
        [[], 2]
    ]
}
```
- To double memory for `call_sSV_Delly` and triple memory for `run_validate_PipeVal` and `call_sSV_Manta`:
```Nextflow
base_resource_update {
    memory = [
        ['call_sSV_Delly', 2],
        [['run_validate_PipeVal', 'call_sSV_Manta'], 3]
    ]
}
```
- To double CPUs and memory for `call_sSV_Manta` and double memory for `run_validate_PipeVal`:
```Nextflow
base_resource_update {
    cpus = [
        ['call_sSV_Manta', 2]
    ]
    memory = [
        [['call_sSV_Manta', 'run_validate_PipeVal'], 2]
    ]
}
```

## Outputs

| Output |	Description |
| ---- | -------- |
| .bcf | Binary VCF output format from DELLY with somatic structural variants if found. |
| .bcf.csi | CSI-format index for BCF files from DELLY. |
| .vcf | Uncompressed VCF output from GRIDSS2 |
| .vcf.idx | Index file for GRIDSS2 VCF |
| .vcf.bgz | Block compressed gzip VCF outputs from GRIDSS2 somatic filtering |
| .vcf.bgz.tbi | Index files for GRIDSS2 somatic filtered outputs |
| .vcf.gz | Compressed VCF output format with somatic structural variants if found. |
| .vcf.gz.tbi | TBI-format index for compressed VCF files. |
| .png | SV Circos plot for individual SV callers (DELLY and Manta) as QC output. |
| report.html, timeline.html and trace.txt | A Nextflow report, timeline and trace files. |
| \*.log.command.* | Process and sample specific logging files created by nextflow. |
| *.sha512 | Generates SHA-512 hash to validate file integrity. |

## Intermediates
| Output |	Description |
| ---- | -------- |
| *assembly.bam | Breakend assembly BAM from GRIDSS2 |
| *assembly.bam.gridss.working | Directory containing BAM metrics from GRIDSS2 breakend assembly |
| *<normal_sample>.gridss.working | Directory containing BAM metrics from GRIDSS2 preprocessing of normal BAM |
| *<tumor_sample>.gridss.working | Directory containing BAM metrics from GRIDSS2 preprocessing of tumor BAM |

## Testing and Validation

### Test Data Set

| Data Set | Normal Sample | Tumor Sample |
| ------ | ------- | ------- |
| A-mini TWGSAMIN000001-T003-S03-F | TWGSAMIN000001-N003-S03-F | TWGSAMIN000001-T003-S03-F |
| ILHNLNEV000009 | ILHNLNEV000009-N001-B01-F | ILHNLNEV000009-T002-L01-F |
| DTB-266_WCDT | DTB-266_DNA_N | DTB-266_DNA_T |

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
* [GRIDSS2 Structural Variant Caller](https://github.com/PapenfussLab/gridss)
* [SVision Structural Variant Caller](https://github.com/xjtu-omics/SVision)


## License

Authors: Yu Pan (YuPan@mednet.ucla.edu), Ghouse Mohammed (GMohammed@mednet.ucla.edu), Mohammed Faizal Eeman Mootor (MMootor@mednet.ucla.edu)

`call-sSV` is licensed under the GNU General Public License version 2. See the file LICENSE for the terms of the GNU GPL license.

`call-sSV` takes BAM files and utilizes DELLY, Manta and GRIDSS2 to call somatic structural variants.

Copyright (C) 2021-2025 University of California Los Angeles ("Boutros Lab") All rights reserved.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
