# covid_rnaseq
Next generation sequencing for SARS‑CoV‑2 RNA from adenocarcinomic human alveolar basal epithelial cells (a595) and normal human bronchial epithelial cells (NHBE). These cell lines were merged in the download of the raw data.
Reference genome taken from Gene Expression Omnibus. Raw data can be found [here](https://www.ncbi.nlm.nih.gov/Traces/study/?query_key=1&WebEnv=MCID_64f899260fed995d09a01357&f=organism_s%3An%3Ahomo%2520sapiens%3Bcell_line_sam_ss_dpl110_ss%3An%3Aa549%2Cnhbe%3Ac&o=acc_s%3Aa).

## Usage
Running `covid_comparison-analysis.R` in its entirety will run the other two R files. It graphs the difference (on the basis of fold change) between the two normalised sets of count data and conducts a correlation test.

`rna-seq.R` separates the samples by cell line. It creates, cleans, and normalises counts for the NHBE site. It then does multiple quality tests, tests for GO enrichment, then writes the annotated results to a file.

`covid-a549-analysis.R` creates, cleans, and normalises counts for the NHBE site. It then does multiple quality tests, tests for GO enrichment.

Running `rna-seq.R` and `covid-a549-analysis.R` individually is encouraged as they provide great insight to the samples, but please be aware the output is quite noisy (many plots).
