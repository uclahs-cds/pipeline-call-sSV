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
The call-sSV pipeline calls somatic structural variants utilizing [Delly](https://github.com/dellytools/delly). This pipeline requires at least one tumor sample and a matched control sample.
This pipeline is developed using Nextflow, docker and can run either on a single node linux machine or a multi-node HPC cluster (e.g. Slurm, SGE).

## How to Run:

1.	Check if the pipeline and the desired version is in /hot/pipeline/release.
2.	Create an instance of the template nextflow.config file. Give the nextflow.config instance with an easy-to-understand name. Modify the nextflow.config instance as needed(the input, the output, and the parameters). Later, you might need to check in the nextflow.config instance. An example can be found [Here](https://github.com/uclahs-cds/pipeline-call-sSV/blob/main/pipeline/config/nextflow.config). See [Inputs](#inputs) for the description of each variables in the config file.
3.	Update the input CSV. See [Inputs](#inputs) for the columns needed. All columns must exist to run the pipeline. An example can be found [Here](https://github.com/uclahs-cds/pipeline-call-sSV/blob/main/pipeline/input/tumor_control_pair_0.csv).
4.	See the submission script, [Here](https://github.com/uclahs-cds/tool-submit-nf), to submit your pipeline

## Flow Diagram:

![](https://github.com/uclahs-cds/pipeline-call-sSV/blob/yupan-documentation/call-sSV-workflow.svg)

## Pipeline Steps:

### Call Somatic Structural Variants:

#### 1. Calling Single Sample Somatic Structural Variants
```script
delly call --genome hg38.fa --exclude hg38.excl --map-qual 20 --min-clique-size 5 --mad-cutoff 15 --outfile t1.bcf tumor1.bam control1.bam
```
This step requires an aligned and sorted tumor sample BAM file and a matched control sample as an input for variant calling with Delly.
The stringent filters (--map-qual 20 --min-clique-size 5 --mad-cutoff 15) are added, which can drastically reduce the runtime, especially when the input bams are big.

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
This step applies somatic filtering against the .bcf file generated in Step 1.

Note: cohort based false positive filtering is compuationally heavy and not implemented in this pipeline.


## Inputs

### Input CSV

The input CSV should have each of the input fields listed below as separate columns, using the same order. An example of the input CSV can be found [here](https://github.com/uclahs-cds/pipeline-call-sSV/blob/main/pipeline/input/tumor_control_pair_0.csv).

| Field |	Type |	Description |
|--- | --- | --- |
|control_sample_bam |	string	| Absolute path to the control sample .BAM file |
|tumor_sample_bam	| string	| Absolute path to the tumor sample .BAM file. |

## Nextflow Config File Parameters
| Input Parameter |	Required |	Type |	Description |
| ------- |   --------- | ------ | -------------|
| dataset_id |	yes	| string |	Boutros Lab dataset id |
| blcds_registered_dataset	| yes |	boolean | Affirms if dataset should be registered in the Boutros Lab Data registry. Default value is false. |
| sge_scheduler	| yes	| boolean	| Affirms whether job will be executed on the SGE cluster. Default value is false. |
| input_csv |	yes |	string	| Absolute path to the input CSV file for the pipeline. |
| reference_fasta	| yes |	path	| Absolute path to the reference genome fasta file. The reference genome is used by Delly for structural variant calling. |
| exclusion_file |	yes	| path |	Absolute path to the delly reference genome exclusion file utilized to remove suggested regions for structural variant calling. |
| map_qual | yes | integer | Min. paired-end (PE) mapping quality |
| min_clique_size | yes | integer | Min. PE/SR clique size |
| mad_cutoff | yes | integer | Insert size cutoff, median+s*MAD (deletions only) |
| save_intermediate_files |	yes	| boolean |	Optional parameter to indicate whether intermediate files will be saved. Default value is true. |
| output_dir |	yes |	path |	Absolute path to the directory where the output files to be saved. |
| temp_dir	| yes	| path |	Absolute path to the directory where the nextflow’s intermediate files are saved. |
| verbose |	false |	boolean	| If set to true, the values of input channels will be printed, can be used for debugging|

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

| Data Set | Run Configuration | Output Dir | Control Sample | Tumor Sample |  
| ------ | ------ | ------- | ------ | ------- |
| A-mini | /hot/pipelines/development/slurm/call-sSV/config/nextflow_amini.config | /hot/pipelines/development/slurm/call-sSV/output_amini/call-sSV-20210920-135858 | /hot/resources/SMC-HET/normal/bams/A-mini/0/output/HG002.N-0.bam | /hot/resources/SMC-HET/tumours/A-mini/bams/0/output/S2.T-0.bam |
| A-full | /hot/pipelines/development/slurm/call-sSV/config/nextflow_afull.config | /hot/pipelines/development/slurm/call-sSV/output_afull/call-sSV-20210921-162552 | /hot/resources/SMC-HET/normal/bams/HG002.N.bam | /hot/pipelines/development/slurm/call-sSV/input/T5.T.sorted_py.bam |
| ILHNLNEV000001-T001-P01-F | /hot/pipelines/development/slurm/call-sSV/config/nextflow_ILHNLNEV000001-T001-P01-F.config | /hot/pipelines/development/slurm/call-sSV/output_ILHNLNEV000001-T001-P01-F_realigned_recalibrated_reheadered/call-sSV-1.0.0 | /hot/users/rhughwhite/ILHNLNEV/call-gSNP/output/2020-12-22/ILHNLNEV000001-T001-P01-F/gSNP/2021-01-05_22.01.08/ILHNLNEV000001/SAMtools-1.10_Picard-2.23.3/recalibrated_reheadered_bam_and_bai/ILHNLNEV000001-N001-B01-F_realigned_recalibrated_reheadered.bam | /hot/users/rhughwhite/ILHNLNEV/call-gSNP/output/2020-12-22/ILHNLNEV000001-T001-P01-F/gSNP/2021-01-05_22.01.08/ILHNLNEV000001/SAMtools-1.10_Picard-2.23.3/recalibrated_reheadered_bam_and_bai/ILHNLNEV000001-T001-P01-F_realigned_recalibrated_reheadered.bam |
| ILHNLNEV000005-T002-L01-F | /hot/users/rhughwhite/ILHNLNEV/call-sSV/inputs_configs/2021-09-10/ILHNLNEV000005-T002-L01-F/nextflow.config | /hot/users/rhughwhite/ILHNLNEV/call-sSV/output/ILHNLNEV000005-T002-L01-F_testing/call-sSV-20210930-180357 | /hot/users/rhughwhite/ILHNLNEV/call-gSNP/output/2020-12-22/ILHNLNEV000005-T002-L01-F/gSNP/2021-01-08_17.01.47/ILHNLNEV000005/SAMtools-1.10_Picard-2.23.3/recalibrated_reheadered_bam_and_bai/ILHNLNEV000005-N001-B01-F_realigned_recalibrated_reheadered.bam | /hot/users/rhughwhite/ILHNLNEV/call-gSNP/output/2020-12-22/ILHNLNEV000005-T002-L01-F//gSNP/2021-01-08_17.01.47/ILHNLNEV000005/SAMtools-1.10_Picard-2.23.3/recalibrated_reheadered_bam_and_bai/ILHNLNEV000005-T002-L01-F_realigned_recalibrated_reheadered.bam | 


### Test runs for the A-mini/partial/full samples were performed using the following reference files

* reference_fasta: /hot/ref/reference/GRCh38-BI-20160721/Homo_sapiens_assembly38.fasta
* reference_fasta_index: /hot/ref/reference/GRCh38-BI-20160721/Homo_sapiens_assembly38.fasta.fai
* exclusion_file: /hot/ref/tool-specific-input/Delly/GRCh38/human.hg38.excl.tsv

## Performance Validation

Testing was performed primarily in the Boutros Lab SLURM Development cluster. Metrics below will be updated where relevant with additional testing and tuning outputs.

|Test Case	| Test Date	| Node Type |	Duration	| CPU Hours	| Virtual Memory Usage (RAM)-peak rss|
|----- | -------| --------| ----------| ---------| --------|
|A-mini(with default filters) | 2021-09-20 | F2 | 16m 24s | 16m 1s	| 1.7 GB |
|A-mini(with stringent filters)	| 2021-10-14 | F2 | 15m 13s | 15m	| 1.7 GB |
|A-full(with default filters) | 2021-09-20 | F72 | 19h 54m 5s | 19h 53m 56s | 15.1 GB |
|ILHNLNEV000001-T001-P01-F (with default filters) | 2021-09-20 | F72 | 22h 30m 16s | 22h 30m 6s | 4.5 GB |
|ILHNLNEV000001-T001-P01-F (with stringent filters)	| 2021-10-14 |	F32 | 8h 37m 34s | 8h 37m 29s | 4.4 GB |
|ILHNLNEV000005-T002-L01-F (with default filters) | 2021-09-20 | F72 | 6d 22h 10m 42s | 11'797.8h | 11.733 GB |
|ILHNLNEV000005-T002-L01-F (with stringent filters. See #10 2f72de1) | 2021-10-02 | F72 | 1d 10h 55m 13s | 2'478.4h | 11.590 GB |


## References

* [Delly Structural Variant Calling](https://github.com/dellytools/delly)


## License

Authors: Yu Pan (YuPan@mednet.ucla.edu), Ghouse Mohammed (GMohammed@mednet.ucla.edu)

Call-sSV is licensed under the GNU General Public License version 2. See the file LICENSE for the terms of the GNU GPL license.

Call-sSV takes BAM files and utilizes Delly to call somatic structural variants.

Copyright (C) 2021 University of California Los Angeles ("Boutros Lab") All rights reserved.

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.