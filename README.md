# PrimalSeq Snakemake Pipeline & Utilities
This repository contains a snakemake pipeline for processing amplicon-based sequencing reads created with the 
[PrimalSeq](https://www.protocols.io/view/primalseq-generation-of-tiled-virus-amplicons-for-bez7jf9n) protocol. 

## Pipeline
The pipeline aligns reads to a reference genome, trims primer sequences and low quality bases, and generates a consensus 
sequence for each sample. A number of alignment statistics are also provided, if specified. Currently, the pipeline is 
used to process data produced as part of the WestNile 4K project, so consensus sequences are also renamed to meet the 
projects specifications. Pipeline parameters can be edited within the `snakemake_config.json` file. 

## Scripts
As of version 4, the PrimalSeq protocol now includes methods for measuring within batch contamination through the use of
sample-specific barcoded spike-ins. This repository contains scripts for estimating the amount of contamination present 
in a sample and identifying and visualizing the sources of contamination, which may be used as part of the snakemake 
pipeline or on their own. 

- **graph_coverage.py**: Plots per base read coverage across the entire reference genome for a given sample.
- **contamination.py**: Identifies barcode reads which have accurately aligned to reference barcode sequences. A wrapper
for `samtools idxstats` which provides ability to filter by alignment score.
- **graph_contamination.py**: Calculates sources of contamination from output of **contamination.py** and visualizes. 