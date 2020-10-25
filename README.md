# metagenomes_NEON

This repository contains code for downloading and processing shotgun metagenomics data from the National Ecological Observatory Network (NEON). This accompanies a tutorial for processing the raw reads.

Two BLAST analyses are included in the tutorial: one focuses on N-cycling genes, and the other on antibiotic resistance genes. There are scripts here for parsing the BLAST output (in Python), writing it to a CSV file, then visualizing the results in R. 

To recreate or update `metadata_metagenome.csv` (the file with raw NEON sample read counts), run this code in R (takes < 3 minutes):

```
library(neonUtilities)
library(dplyr)
metadata <- loadByProduct(dpID = 'DP1.10107.001', check.size = FALSE, package = 'expanded')
sequence_metadata <- metadata$mms_metagenomeSequencing %>%
	select(dnaSampleID, sampleTotalReadNumber, sampleFilteredReadNumber) %>%
	distinct() %>% rename(sampleID = dnaSampleID)
saveRDS(sequence_metadata, "./metadata_metagenome.csv")
```
