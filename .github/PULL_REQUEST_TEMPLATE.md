<!--- Please read each of the following items and confirm by replacing
 !--the [ ] with a [X] --->

- [ ] I have read the [code review guidelines](https://confluence.mednet.ucla.edu/display/BOUTROSLAB/Code+Review+Guidelines) and the [code review best practice on GitHub check-list](https://confluence.mednet.ucla.edu/pages/viewpage.action?pageId=84091668).

- [ ] I have set up the branch protection rule following the [github standards](https://confluence.mednet.ucla.edu/pages/viewpage.action?spaceKey=BOUTROSLAB&title=GitHub+Standards#GitHubStandards-Branchprotectionrule) before opening this pull request, or the branch protection rule has already been set up.

- [ ] I have added my name to the contributors listings in the
``metadata.yaml`` and the ``manifest`` block in the config as part of this pull request, am listed
already, or do not wish to be listed. (*This acknowledgement is optional.*)

- [ ] I have added the changes included in this pull request to the `CHANGELOG.md` under the next release version or unreleased, and updated the date.

- [ ] I have updated the version number in the `metadata.yaml` and config file following [semver](https://semver.org/), or the version number has already been updated. (*Leave it unchecked if you are unsure about new version number and discuss it with the infrastructure team in this PR.*)

- [ ] I have tested the pipeline on at least one DNA A-mini sample, one A-full sample, and one Rupert WGS sample. The paths to the test config files and output directories were attached below.

<!--- Briefly describe the changes included in this pull request and the paths to the test cases below
 !--- starting with 'Closes #...' if appropriate --->

Closes #...

**Test Results**

- DNA A-mini
	- sample:    <!-- e.g. A-mini S2.T-0, HG002.N-0 -->
	- input:     <!-- /hot/pipelines/development/slurm/call-sSV/input/tumor_control_pair_0.csv -->
	- config:    <!-- /hot/pipelines/development/slurm/call-sSV/config/nextflow_amini.config -->
	- output:    <!-- /hot/pipelines/development/slurm/call-sSV/output_amini/call-sSV-20210920-135858 --> 

- DNA A-full
	- sample:    <!-- e.g. T5.T, HG002.N -->
	- input:     <!-- /hot/pipelines/development/slurm/call-sSV/input/tumor_control_pair_afull.csv -->
	- config:    <!-- /hot/pipelines/development/slurm/call-sSV/config/nextflow_afull.config -->
	- output:    <!-- /hot/pipelines/development/slurm/call-sSV/output_afull/call-sSV-20210921-162552 --> 

- Rupert WGS sample
	- sample:    <!-- e.g. ILHNLNEV000001-T001-P01-F, ILHNLNEV000001-N001-B01-F -->
	- input:     <!-- /hot/pipelines/development/slurm/call-sSV/input/tumor_control_pair_rupert_WGS_real_sample.csv -->
	- config:    <!-- /hot/pipelines/development/slurm/call-sSV/config/nextflow_rupert_WGS_real_sample.config -->
	- output:    <!-- /hot/pipelines/development/slurm/call-sSV/output_rupert_WGS_real_sample/call-sSV-20210920-153803 --> 