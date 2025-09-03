# Changelog

All notable changes to the pipeline-name pipeline.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- `SVision` caller

### Changed

- Update Delly `v1.2.6` to `v1.5.0`
- Remove cluster paths from README

## [7.0.0] - 2025-05-16

### Added

- Add workflows to compress VCF and convert BCF2VCF
- Add GRIDSS2 somatic filter
- Add Resource Handler
- Add GRIDSS2 variant calling
- Add GRIDSS2 assembly
- Add GRIDSS2 preprocessing
- Add supported Nextflow version to `README.md`
- Add PlantUML diagram
- Add call to `store_object_as_json` to store parameters
- CIRCOS plotting for Manta

### Changed

- Update README to document GRIDSS2
- Update pipeline DAG to include GRIDSS2
- Enable CIRCOS plotting for DELLY
- Update PlantUML action to `v1.0.1`
- Update memory allocations in `M64.config`
- Re-enable configuration testing workflows

### Fixed

- Fix NFTest cases after GRIDSS2 additions
- Fix configuration tests after GRIDSS2 additions

## [6.1.0] - 2024-03-12

### Added

- `CODEOWNERS` file

### Removed

- Local versions of modularized functions

## [6.0.0] - 2024-02-23

### Changed

- Perform additional pipeline runs towards release `v6.0.0`

## [6.0.0-rc.2] - 2024-02-06

### Added

- Add NFTest

### Changed

- Update DELLY `v1.1.3` to `v1.2.6`
- Update submodules

### Removed

- Remove `set_env` method code in `default/methods.config`

## [6.0.0-rc.1] - 2023-08-02

### Added

- YAML input

### Changed

- Update to latest commit of \`pipeline-Nextflow-config\`\`
- Update pipeline SVG
- Update YAML input and tests in `README.md`
- Update PipeVal to `v4.0.0-rc.2`
- Parse sample ID from tumor BAM for output directory naming
- Update `README.md` to clarify adjustable parameters and note lab default values

### Removed

- CSV input

## [5.0.0] - 2023-01-27

### Changed

- Update `README.md` for release `5.0.0`
- Perform additional test to release `5.0.0`

## [5.0.0-rc.1] - 2023-01-24

### Added

- Add `PipeVal:v3.0.0` using pipeline-Nextflow-module to validate the input CSV file
- Parameter `algorithm` to select the SV caller(s) of choice in `template.config`
- Add Manta SV caller

### Changed

- Update `README.md` for release `5.0.0-rc.1`
- Move param checking to `methods.config` using `schema.config`
- Parameterize Docker registry
- Use `ghcr.io/uclahs-cds` as default registry
- Simplify output directory declaration for each process using `addParams` and use `output_dir_base`

### Removed

- Remove `module/validation.nf` as PipeVal sub-module is used

## [4.0.0] - 2022-10-13

### Added

- Retry mechanism to retry a process with higher memory - PR #69

### Changed

- Standardize the variable names in the pipeline and change input CSV column names from `control_sample_bam` to `normal_bam` and `tumor_sample_bam` to `tumor_bam`
- Standardize output filename using [generate_standardized_filename](https://github.com/uclahs-cds/pipeline-Nextflow-module/tree/main/modules/common/generate_standardized_filename)
- Upgrade Delly 1.0.3 to 1.1.3

## [3.0.0] - 2022-07-25

### Added

- Add F16.config to enable F16 node compatibility

### Changed

- Rename `call-sSV.nf` to `main.nf` and replace `mainScript` with `name` in `main.nf` line 15
- Upgrade BCTools 1.12 to 1.15.1
- Remove `save_intermediate_files` argument in default.config as it is already specified in template.config
- Standardize the repo structure based on NF pipeline template
- Update <http://README.md> with resource locations and runtime of Delly v1.0.3 on different datasets and node types
- Change Delly version from 0.9.1 to 1.0.3 in default.config

### Fixed

- Fix memory allocation for process 'query_SampleName_BCFtools' in F2.config - Discussion #28 and PR #33
- Fix F2 detection in methods.config - PR #32

## [2.0.0] - 2022-02-17

### Added

- Add the process "filter_BCF_BCFtools" to filter Non-Pass calls.

## [1.0.1] - 2021-12-09

### Changed

- Change Delly version from 0.8.7 to 0.9.1 in nextflow.config and metadata.yaml
- Use sha512sum to replace 'python -m validate -t sha512-gen'

## [1.0.0] - 2021-10-29

### Changed

- Minor changes of wording in <http://CHANGELOG.md> and metadata.yaml

## [1.0.0-rc1] - 2021-10-18

### Added

- Use Delly to call the somatic structural variants (SVs) of one tumor sample and a matched control sample, further filter the somatic SVs.
- Add option to specify `map_qual`, `min_clique_size`, and `mad_cutoff` Delly filters to reduce the runtime.
- Use BCFtools to automatically generate the samples.tsv that is needed by the sSV filtering (see <https://github.com/dellytools/delly>).

[1.0.0]: https://github.com/uclahs-cds/pipeline-call-sSV/compare/v1.0.0-rc1...v1.0.0
[1.0.0-rc1]: https://github.com/uclahs-cds/pipeline-call-sSV/releases/tag/v1.0.0-rc1
[1.0.1]: https://github.com/uclahs-cds/pipeline-call-sSV/compare/v1.0.0...v1.0.1
[2.0.0]: https://github.com/uclahs-cds/pipeline-call-sSV/compare/v1.0.1...v2.0.0
[3.0.0]: https://github.com/uclahs-cds/pipeline-call-sSV/compare/v2.0.0...v3.0.0
[4.0.0]: https://github.com/uclahs-cds/pipeline-call-sSV/compare/v3.0.0...v4.0.0
[5.0.0]: https://github.com/uclahs-cds/pipeline-call-sSV/compare/v5.0.0-rc.1...v5.0.0
[5.0.0-rc.1]: https://github.com/uclahs-cds/pipeline-call-sSV/compare/v4.0.0...v5.0.0-rc.1
[6.0.0]: https://github.com/uclahs-cds/pipeline-call-sSV/compare/v6.0.0-rc.2...v6.0.0
[6.0.0-rc.1]: https://github.com/uclahs-cds/pipeline-call-sSV/compare/v5.0.0...v6.0.0-rc.1
[6.0.0-rc.2]: https://github.com/uclahs-cds/pipeline-call-sSV/compare/v6.0.0-rc.1...v6.0.0-rc.2
[6.1.0]: https://github.com/uclahs-cds/pipeline-call-sSV/compare/v6.0.0...v6.1.0
[7.0.0]: https://github.com/uclahs-cds/pipeline-call-sSV/compare/v6.1.0...v7.0.0
