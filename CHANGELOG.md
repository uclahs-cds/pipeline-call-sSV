# Changelog
All notable changes to the pipeline-name pipeline.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

---

## [5.0.0] - 2023-01-25
### Changed
- Update `README.md` for release `5.0.0`
- Update `CHANGELOG.md`
- Perform additional test to release `5.0.0`

---

## [5.0.0-rc.1] - 2023-01-24
### Changed
- Update `README.md` for release `5.0.0-rc.1`
- Move param checking to `methods.config` using `schema.config`
- Parameterize Docker registry
- Use `ghcr.io/uclahs-cds` as default registry
- Simplify output directory declaration for each process using `addParams` and use `output_dir_base`

### Added
- Add `PipeVal:v3.0.0` using pipeline-Nextflow-module to validate the input CSV file
- Parameter `algorithm` to select the SV caller(s) of choice in `template.config`
- Add Manta SV caller

### Removed
- Remove `module/validation.nf` as PipeVal sub-module is used

---

## [4.0.0] - 2022-10-13
### Added
- Retry mechanism to retry a process with higher memory - PR #69

### Changed
- Standardize the variable names in the pipeline and change input CSV column names from `control_sample_bam` to `normal_bam` and `tumor_sample_bam` to `tumor_bam`
- Standardize output filename using [generate_standardized_filename](https://github.com/uclahs-cds/pipeline-Nextflow-module/tree/main/modules/common/generate_standardized_filename)
- Upgrade Delly 1.0.3 to 1.1.3

---

## [3.0.0] - 2022-07-25
### Added
- Add F16.config to enable F16 node compatibility
### Changed
- Rename `call-sSV.nf` to `main.nf` and replace `mainScript` with `name` in `main.nf` line 15
- Upgrade BCTools 1.12 to 1.15.1
- Remove `save_intermediate_files` argument in default.config as it is already specified in template.config
- Standardize the repo structure based on NF pipeline template
- Update README.md with resource locations and runtime of Delly v1.0.3 on different datasets and node types
- Change Delly version from 0.9.1 to 1.0.3 in default.config
### Fixed
- Fix memory allocation for process 'query_SampleName_BCFtools' in F2.config - Discussion #28 and PR #33
- Fix F2 detection in methods.config - PR #32

---

## [2.0.0] - 2022-02-17
### Added
- Add the process "filter_BCF_BCFtools" to filter Non-Pass calls.

---

## [1.0.1] - 2021-12-09
### Changed
- Change Delly version from 0.8.7 to 0.9.1 in nextflow.config and metadata.yaml
- Use sha512sum to replace 'python -m validate -t sha512-gen'

---

## [1.0.0] - 2021-10-29
### Changed
- Minor changes of wording in CHANGELOG.md and metadata.yaml

## [1.0.0-rc1] - 2021-10-18
### Added
- Use Delly to call the somatic structural variants (SVs) of one tumor sample and a matched control sample, further filter the somatic SVs.
- Add option to specify `map_qual`, `min_clique_size`, and `mad_cutoff` Delly filters to reduce the runtime.
- Use BCFtools to automatically generate the samples.tsv that is needed by the sSV filtering (see https://github.com/dellytools/delly).

---
