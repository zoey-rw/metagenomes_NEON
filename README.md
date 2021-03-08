# metagenomes_NEON

This repository contains code for downloading and processing shotgun metagenomics data from the National Ecological Observatory Network (NEON). This accompanies a tutorial for processing the raw reads, accessible on F1000 Research.

Two BLAST analyses are described in the tutorial: one focuses on N-cycling genes, and the other on antibiotic resistance genes. There is an example script here for parsing the BLAST output (in Python) and writing it to a CSV file, `blast_parser.py`. Then, the scripts beginning with `normalize_gene_counts` can be used to normalize data using the DeSeq R package and visualize the results in R. 

To update `metadata_metagenome.csv` (the file with raw NEON sample read counts - necessary for normalization), run this code in R (takes < 3 minutes):

```
library(neonUtilities)
library(dplyr)
metadata <- loadByProduct(dpID = 'DP1.10107.001', check.size = FALSE, package = 'expanded')
sequence_metadata <- metadata$mms_metagenomeSequencing %>%
	select(dnaSampleID, sampleTotalReadNumber, sampleFilteredReadNumber) %>%
	distinct() %>% rename(sampleID = dnaSampleID)
write.csv(sequence_metadata, "./metadata_metagenome.csv")
```
