# load libraries
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("/projectnb/talbot-lab-data/zrwerbin/metagenomes_NEON")

### Specify data path - must change for your own file system ###
# Path to BLASTp or BLASTx results
blast_path <- '/projectnb/talbot-lab-data/zrwerbin/metagenomes_raw/test_pipeline/sunbeam_output/BLAST/blastp_ncyc.csv'
blastp_raw <- data.table::fread(blast_path, data.table = F)

# Path to NCycDB key (included in Github repo)
N_pathway_key <- read.csv('https://raw.githubusercontent.com/zoey-rw/metagenomes_NEON/main/NCycDB_pathways.csv', 
													stringsAsFactors = F, strip.white = T)  
# Path to metagenome metadata, NEON DP1.10107.001 (included in Github repo) - must download for total read counts
sequence_metadata <- read.csv("https://raw.githubusercontent.com/zoey-rw/metagenomes_NEON/main/metadata_metagenome.csv")

# Separate BLAST results into columns
# Note: Does not handle the "No definiton line" results well.
blastp <- blastp_raw %>% 
	separate(V1, sep=".x", into=c("sampleID",NA)) %>% 
	mutate(V2 = gsub("\\[|\\]| homolog|", "", V2)) %>% 
  separate(V2, sep=c(" "), into = c("ID","gene","ontology","source")) %>% 
  separate(gene, "=", into=c(NA, "gene")) %>% 
  separate(ontology, "=", into=c(NA, "ontology")) %>% 
	separate(source, "=", into=c(NA, "source")) %>% 
	rename(alignment_length = V3, e_val = V4)

# Add pathway data and summarize gene hits by pathway
blastp$pathway <- N_pathway_key[match(blastp$gene, N_pathway_key$Gene_family),]$Pathway
blastp_counts <- blastp %>% filter(!grepl("Others", pathway) & !is.na(pathway)) %>% 
	group_by(sampleID, pathway) %>% add_count(name = "gene_count") %>%
	mutate("distinct_genes" = n_distinct(gene)) %>% 
	distinct(sampleID, gene_count, distinct_genes)

# Now merge with metadata to get a normalized gene count.
merged <- merge(blastp_counts, sequence_metadata) %>% 
	mutate(norm = (gene_count/sampleTotalReadNumber),
				 siteID = substr(sampleID, 1, 4)) 

#merged$sampleID <- gsub("-COMP-DNA[12]","",merged$sampleID)

### Let's make some plots! ###

# Visualize normalized number of hits per site
ggplot(merged) + geom_boxplot(aes(x = siteID, y = norm, fill = siteID)) + 
	facet_wrap(~pathway) + ggtitle("Nitrogen cycling genes, by pathway", 
																 "Normalized by total reads per sample") + 
	xlab("Site") + ylab("Relative abundance") + theme_bw() + 
	theme(axis.text.x  = element_text(angle=290, hjust=-.05, size=12))
	
# Visualize number of distinct genes
ggplot(merged) + geom_point(aes(x = sampleTotalReadNumber, y = distinct_genes, color = siteID)) + 
	facet_wrap(~pathway) + ggtitle("Unique N-cycling genes per pathway", 
																 "As a function of total reads per sample") + 
	xlab("Sequencing depth") + ylab("Number of unique N-cycling genes") + theme_bw() + 
	theme(axis.text.x  = element_text(angle=290, hjust=-.05, size=12))


# blastx_raw <- read.table('/projectnb/talbot-lab-data/zrwerbin/metagenomes_raw/test_pipeline/blastx.csv')
# blastx <- blastx_raw %>% 
#   separate(V1, sep=c(","), into = c("ID","gene","homolog","ontology","source")) %>% 
#   separate(gene, "=", into=c(NA, "gene")) %>% 
#   separate(ontology, "=", into=c(NA, "ontology"))
# blastx$pathway <- map[match(blastx$gene, map$Gene_family),]$Pathway