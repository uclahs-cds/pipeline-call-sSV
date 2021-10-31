# Changelog
All notable changes to the pipeline-name pipeline.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

---

## [1.0.0] - 2020-10-29
### Added
- Use Delly to call the somatic structural variants (SVs) of one tumor sample and a matched control sample, further filter the somatic SVs.
- Add option to specify `map_qual`, `min_clique_size`, and `mad_cutoff` Delly filters to reduce the runtime.
- Use BCFtools to automatically generate the samples.tsv that is needed by the sSV filtering (see https://github.com/dellytools/delly).

---