# load libraries
library(dplyr)          
library(tidyr)   
library(data.table)
library(ggplot2)

### Specify paths - must change for your own file system ###
# Path to NEON DP1.10107.001 (metagenome metadata) - must download for total read counts
sequence_metadata_path <- "/projectnb/talbot-lab-data/zrwerbin/metagenomes_raw/test_pipeline/metadata_metagenome.rds"
# Path to Antibiotic Resistance Ontology (ARO) - download by accessing https://card.mcmaster.ca/latest/ontology 
aro_path <- "/projectnb/talbot-lab-data/zrwerbin/metagenomes_NEON/ARO/aro.tsv"
# Path to BLASTp or BLASTx results
blast_path <- '/projectnb/talbot-lab-data/zrwerbin/metagenomes_raw/test_pipeline/blastp_card.csv'

### Read in the above files. ###
blastp_raw <- fread(blast_path, header = F, data.table = F)
sequence_metadata <- readRDS(sequence_metadata_path)
aro <- fread(aro_path, sep = "\t", data.table = F)

# Subset ARO to genes that mention "tetracycline." Approx. 103 genes.
tetra <- aro %>% filter(grepl("tetracycline", Description)) 

# Separate BLAST results into columns, subset to tetracycline gene hits, 
# and count total number of gene hits as well as number of distinct genes
# Note: does not handle the "No definiton line" results well.
blastp_counts <- blastp_raw %>% 
	separate(V1, sep=".x", into=c("sampleID",NA)) %>% 
	separate(V6, sep=c("]|,"), into = c("blast_hit",NA,"alignment_length","e_val")) %>% 
	rename(Accession = V5) %>% 
	filter(Accession %in% tetra$Accession) %>% 
	group_by(sampleID) %>% add_count(name = "gene_count") %>%
	mutate("distinct_genes" = n_distinct(Accession)) %>% 
	distinct(sampleID, gene_count, distinct_genes)

# Now merge with metadata to get a normalized gene count.
merged <- merge(blastp_counts, sequence_metadata) %>% 
	mutate(norm = (gene_count/sampleTotalReadNumber)*10000, 
				 siteID = substr(sampleID, 1, 4)) 

### Let's make some plots! ###

# Visualize normalized number of hits per site
ggplot(merged) + geom_boxplot(aes(x = siteID, y = norm, fill = siteID)) + 
	ggtitle("Tetracycline resistance gene abundance", "Normalized by total reads per sample") + 
	xlab("Sample") + ylab("Relative abundance") + theme_bw() + 
	theme(axis.text.x  = element_text(angle=290, hjust=-.05, size=12))

# Visualize number of distinct genes
ggplot(merged) + geom_point(aes(x = sampleTotalReadNumber, y = distinct_genes, color = siteID)) + 
	ggtitle("Unique tetracycline resistance genes", "As a function of total reads per sample") + 
	xlab("Sequencing depth") + ylab("Number of unique tetracycline-resistance genes") + theme_bw() + 
	theme(axis.text.x  = element_text(angle=290, hjust=-.05, size=12))
