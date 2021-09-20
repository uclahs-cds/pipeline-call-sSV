# Changelog
All notable changes to the pipeline-name pipeline.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]
### Added
- Use Delly to call the SV of one tumor sample and a matched control sample, further filter the somatic SV.
- The sSV filtering needs a samples.tsv (see https://github.com/dellytools/delly). Use bcftools to automatically generate it.
- Implemented the cohort based false positive filtering. However, it turned out this is computationally infeasible. After discussing this with Taka, we decide to drop this step. 

---


