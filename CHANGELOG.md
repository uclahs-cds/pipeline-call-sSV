# Changelog
All notable changes to the pipeline-name pipeline.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

---

## [1.0.0] - 2020-10-29
### Added
- Use Delly to call the SV of one tumor sample and a matched control sample, further filter the somatic SV.
- Added option to specify map_qual, min_clique_size, and mad_cutoff Delly filters to reduce the runtime.
- The sSV filtering needs a samples.tsv (see https://github.com/dellytools/delly). Use BCFtools to automatically generate it.
- Implemented the cohort based false positive filtering. However, it turned out the cohort based false positive filtering is compuationally heavy, so it has been excluded from this pipeline.

---