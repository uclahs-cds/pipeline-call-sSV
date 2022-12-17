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
The call-sSV pipeline calls somatic structural variants utilizing [Delly](https://github.com/dellytools/delly). This pipeline requires at least one tumor sample and a matched normal sample.
This pipeline is developed using Nextflow, docker and can run either on a single node linux machine or a multi-node HPC Slurm cluster.

## How to Run:

Below is a summary of how to run the pipeline.  See [here](https://confluence.mednet.ucla.edu/pages/viewpage.action?spaceKey=BOUTROSLAB&title=How+to+run+a+nextflow+pipeline) for full instructions.

Pipelines should be run **WITH A SINGLE SAMPLE AT A TIME**. Otherwise resource allocation and Nextflow errors could cause the pipeline to fail.

1. The recommended way of running the pipeline is to directly use the source code located here: `/hot/software/pipeline/pipeline-call-sSV/Nextflow/release/`, rather than cloning a copy of the pipeline.

    * The source code should never be modified when running our pipelines

2. Create a config file for input, output, and parameters. An example for a config file can be found [here](config/template.config). See [Nextflow Config File Parameters](#Nextflow-Config-File-Parameters) for the detailed description of each variable in the config file.

    * Do not directly modify the source `template.config`, but rather you should copy it from the pipeline release folder to your project-specific folder and modify it there

3. Create the input CSV using the [template](input/call-sSV-input.csv).See [Input CSV](#Input-CSV) for detailed description of each column. All columns must exist and should be comma separated in order to run the pipeline successfully.
   
   * Again, do not directly modify the source template CSV file.  Instead, copy it from the pipeline release folder to your project-specific folder and modify it there.

4. The pipeline can be executed locally using the command below:

```bash
nextflow run path/to/main.nf -config path/to/sample-specific.config
```

* For example, `path/to/main.nf` could be: `/hot/software/pipeline/pipeline-call-sSV/Nextflow/release/3.0.0/main.nf`
* `path/to/sample-specific.config` is the path to where you saved your project-specific copy of [template.config](config/template.config) 

To submit to UCLAHS-CDS's Azure cloud, use the submission script [here](https://github.com/uclahs-cds/tool-submit-nf) with the command below:

```bash
python path/to/submit_nextflow_pipeline.py \
    --nextflow_script path/to/main.nf \
    --nextflow_config path/to/sample-specific.config \
    --pipeline_run_name <sample_name> \
    --partition_type F16 \
    --email <your UCLA email, jdoe@ucla.edu>
```
In the above command, the partition type can be changed based on the size of the dataset. At this point, node F16 is generally recommended for larger datasets like A-full and node F2 for smaller datasets like A-mini.

> **Note**: Because this pipeline uses an image stored in the GitHub Container Registry, you must follow the steps listed in the [Docker Introduction](https://confluence.mednet.ucla.edu/display/BOUTROSLAB/Docker+Introduction#DockerIntroduction-GitHubContainerRegistryGitHubContainerRegistry|Setup) on Confluence to set up a PAT for your GitHub account and log into the registry on the cluster before running this pipeline.

## Flow Diagram:

![](https://github.com/uclahs-cds/pipeline-call-sSV/blob/yupan-documentation/call-sSV-workflow.svg)

## Pipeline Steps:

### Call Somatic Structural Variants:

#### 1. Calling Single Sample Somatic Structural Variants
```script
delly call --genome hg38.fa --exclude hg38.excl --map-qual 20 --min-clique-size 5 --mad-cutoff 15 --outfile t1.bcf tumor1.bam normal1.bam
```
This step requires an aligned and sorted tumor sample BAM file and a matched normal sample as an input for variant calling with Delly.
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


## Inputs

### Input CSV

The input CSV should have each of the input fields listed below as separate columns, using the same order and comma as column separator. An example of the input CSV can be found [here](input/call-sSV-input.csv).

| Field |	Type |	Description |
|--- | --- | --- |
|normal_bam |	string	| Absolute path to the normal sample `.bam` file |
|tumor_bam	| string	| Absolute path to the tumor sample `.bam` file. |

## Nextflow Config File Parameters

| Input Parameter |	Required |	Type |	Description |
| ------- |   --------- | ------ | -------------|
| dataset_id |	yes	| string |	Boutros Lab dataset id |
| blcds_registered_dataset	| yes |	boolean | Affirms if dataset should be registered in the Boutros Lab Data registry. Default value is `false`. |
| algorithm | yes | list | List containing a combination of SV callers `delly`, `manta`. List can contain a single caller of choice.  | 
| `run_delly` | true | boolean | Whether or not the workflow should run Delly (either run_delly or run_manta must be set to `true`) |
| `run_manta` | true | boolean | Whether or not the workflow should run Manta (either run_delly or run_manta must be set to `true`) |
| input_csv |	yes |	string	| Absolute path to the input CSV file for the pipeline. |
| reference_fasta	| yes |	path	| Absolute path to the reference genome FASTA file. The reference genome is used by Delly for structural variant calling. GRCh37 - /hot/ref/reference/GRCh37-EBI-hs37d5/hs37d5.fa, GRCh38 - /hot/ref/reference/GRCh38-BI-20160721/Homo_sapiens_assembly38.fasta |
| exclusion_file |	yes	| path | Absolute path to the delly reference genome exclusion file utilized to remove suggested regions for structural variant calling. GRCh37 - /hot/ref/tool-specific-input/Delly/GRCh37-EBI-hs37d/human.hs37d5.excl.tsv, GRCh38 - /hot/ref/tool-specific-input/Delly/hg38/human.hg38.excl.tsv |
| map_qual | yes | integer | Min. paired-end (PE) mapping quality |
| min_clique_size | yes | integer | Min. PE/SR clique size |
| mad_cutoff | yes | integer | Insert size cutoff, median+s*MAD (deletions only) |
| save_intermediate_files |	yes	| boolean |	Optional parameter to indicate whether intermediate files will be saved. Default value is `false`. |
| output_dir |	yes |	path |	Absolute path to the directory where the output files to be saved. |
| work_dir	| no	| path |	Path of working directory for Nextflow. When included in the sample config file, Nextflow intermediate files and logs will be saved to this directory. With `ucla_cds`, the default is `/scratch` and should only be changed for testing/development. Changing this directory to `/hot` or `/tmp` can lead to high server latency and potential disk space limitations, respectively. |
| verbose |	yes |	boolean	| If set to `true`, the values of input channels will be printed, can be used for debugging|
| `docker_container_registry` | optional | string | Registry containing tool Docker images. Default: `ghcr.io/uclahs-cds` |

An example of the NextFlow Input Parameters Config file can be found [here](config/template.config).

## Outputs

| Output |	Output type |	Description |
| ---- | ----- | -------- |
| .bcf |	final	| Binary VCF output format with somatic structural variants if found. |
| .bcf.csi	| final	| CSI-format index for BCF files |
| report.html, timeline.html and trace.txt	| log |	A Nextflow report, timeline and trace files. |
| \*.log.command.*	| log |	Process and sample specific logging files created by nextflow. |
| *.sha512 |	checksum |	Generates SHA-512 hash to validate file integrity. |


## Testing and Validation

### Test Data Set

| Data Set | Run Configuration | Output Dir | Normal Sample | Tumor Sample |  
| ------ | ------ | ------- | ------ | ------- |
| A-mini | /hot/software/pipeline/pipeline-call-sSV/Nextflow/development/3.0.0/mmoootor-upgrade-delly-0.9.1-to-1.0.3/config/A-mini-hg38.config | /hot/software/pipeline/pipeline-call-sSV/Nextflow/development/3.0.0/mmoootor-upgrade-delly-0.9.1-to-1.0.3/A-mini/call-sSV-2.0.0/S2.T-0/DELLY-1.0.3/output/ | /hot/resource/SMC-HET/normal/bams/A-mini/0/output/HG002.N-0.bam | /hot/resource/SMC-HET/tumours/A-mini/bams/0/output/S2.T-0.bam |
| A-full | /hot/software/pipeline/pipeline-call-sSV/Nextflow/development/3.0.0/mmoootor-upgrade-delly-0.9.1-to-1.0.3/config/A-full-F72-hg19.config | /hot/software/pipeline/pipeline-call-sSV/Nextflow/development/3.0.0/mmoootor-upgrade-delly-0.9.1-to-1.0.3/A-full-F72/call-sSV-2.0.0/T5.T.sorted_py/DELLY-1.0.3/output/ | /hot/resource/SMC-HET/normal/bams/HG002.N.bam | /hot/software/pipeline/pipeline-call-sSV/Nextflow/development/input/data/T5.T.sorted_py.bam |
| ILHNLNEV000001-T001-P01-F | /hot/software/pipeline/pipeline-call-sSV/Nextflow/development/3.0.0/mmoootor-upgrade-delly-0.9.1-to-1.0.3/config/ILHNLNEV000001-T001-P01-F-F32-hg38.config | /hot/software/pipeline/pipeline-call-sSV/Nextflow/development/3.0.0/mmoootor-upgrade-delly-0.9.1-to-1.0.3/ILHNLNEV000001-T001-P01-F-F32/call-sSV-2.0.0/ILHNLNEV000001-T001-P01-F_realigned_recalibrated_reheadered/DELLY-1.0.3/output/ | /hot/user/rhughwhite/ILHNLNEV/call-gSNP/output/2020-12-22/ILHNLNEV000001-T001-P01-F/gSNP/2021-01-05_22.01.08/ILHNLNEV000001/SAMtools-1.10_Picard-2.23.3/recalibrated_reheadered_bam_and_bai/ILHNLNEV000001-N001-B01-F_realigned_recalibrated_reheadered.bam | /hot/user/rhughwhite/ILHNLNEV/call-gSNP/output/2020-12-22/ILHNLNEV000001-T001-P01-F/gSNP/2021-01-05_22.01.08/ILHNLNEV000001/SAMtools-1.10_Picard-2.23.3/recalibrated_reheadered_bam_and_bai/ILHNLNEV000001-T001-P01-F_realigned_recalibrated_reheadered.bam |
| ILHNLNEV000005-T002-L01-F | /hot/user/rhughwhite/ILHNLNEV/call-sSV/inputs_configs/2021-09-10/ILHNLNEV000005-T002-L01-F/nextflow.config | /hot/user/rhughwhite/ILHNLNEV/call-sSV/output/ILHNLNEV000005-T002-L01-F_testing/call-sSV-20210930-180357/ | /hot/user/rhughwhite/ILHNLNEV/call-gSNP/output/2020-12-22/ILHNLNEV000005-T002-L01-F/gSNP/2021-01-08_17.01.47/ILHNLNEV000005/SAMtools-1.10_Picard-2.23.3/recalibrated_reheadered_bam_and_bai/ILHNLNEV000005-N001-B01-F_realigned_recalibrated_reheadered.bam | /hot/user/rhughwhite/ILHNLNEV/call-gSNP/output/2020-12-22/ILHNLNEV000005-T002-L01-F/gSNP/2021-01-08_17.01.47/ILHNLNEV000005/SAMtools-1.10_Picard-2.23.3/recalibrated_reheadered_bam_and_bai/ILHNLNEV000005-T002-L01-F_realigned_recalibrated_reheadered.bam |

## Performance Validation

Testing was performed primarily in the Boutros Lab SLURM Development cluster. Metrics below will be updated where relevant with additional testing and tuning outputs.

## with Delly v1.0.3 and newer versions
|Test Case      | Test Date     | Node Type |   Duration        | CPU Hours     | Virtual Memory Usage (RAM)-peak rss|
|----- | -------| --------| ----------| ---------| --------|
|A-mini(with stringent filters) | 2022-07-01 | F2 | 20m 32s | 18m | 1.8 GB |
|A-full(with stringent filters) | 2022-07-10 | F16 | 17h 53m 49s | 17h 54m | 15.1 GB |
|A-full(with stringent filters) | 2021-09-20 | F32 | 20h 14m 1s | 20h 12m | 15.1 GB |
|A-full(with stringent filters) | 2022-07-10 | F72 | 18h 16m 15s | 18h 18m | 15.1 GB |
|ILHNLNEV000001-T001-P01-F (with stringent filters)     | 2022-07-10 |  F16 | 8h 46m 31s | 8h 48m | 4.4 GB |
|ILHNLNEV000001-T001-P01-F (with stringent filters)     | 2022-07-10 |  F32 | 9h 38m 29s | 9h 36m | 4.4 GB |

## with Delly v0.9.1 and older versions
|Test Case	| Test Date	| Node Type |	Duration	| CPU Hours	| Virtual Memory Usage (RAM)-peak rss|
|----- | -------| --------| ----------| ---------| --------|
|A-mini(with default filters) | 2021-09-20 | F2 | 16m 24s | 16m 1s	| 1.7 GB |
|A-mini(with stringent filters)	| 2021-10-14 | F2 | 15m 13s | 15m	| 1.7 GB |
|A-full(with default filters) | 2021-09-20 | F72 | 19h 54m 5s | 19h 53m 56s | 15.1 GB |
|ILHNLNEV000001-T001-P01-F (with default filters) | 2021-09-20 | F72 | 22h 30m 16s | 22h 30m 6s | 4.5 GB |
|ILHNLNEV000001-T001-P01-F (with stringent filters)	| 2021-10-14 |	F32 | 8h 37m 34s | 8h 37m 29s | 4.4 GB |
|ILHNLNEV000005-T002-L01-F (with default filters) | 2021-09-20 | F72 | 6d 22h 10m 42s | 11'797.8h | 11.733 GB |
|ILHNLNEV000005-T002-L01-F (with stringent filters. See [#10](https://github.com/uclahs-cds/pipeline-call-sSV/issues/10) [2f72de1](https://github.com/uclahs-cds/pipeline-call-sSV/commit/2f72de1ba190623e4344f144a12cc315fda1dd18)) | 2021-10-02 | F72 | 1d 10h 55m 13s | 2'478.4h | 11.590 GB |


## References

* [Delly Structural Variant Calling](https://github.com/dellytools/delly)


## License

Authors: Yu Pan (YuPan@mednet.ucla.edu), Ghouse Mohammed (GMohammed@mednet.ucla.edu), Mohammed Faizal Eeman Mootor (MMootor@mednet.ucla.edu)

Call-sSV is licensed under the GNU General Public License version 2. See the file LICENSE for the terms of the GNU GPL license.

Call-sSV takes BAM files and utilizes Delly to call somatic structural variants.

Copyright (C) 2021-2022 University of California Los Angeles ("Boutros Lab") All rights reserved.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
