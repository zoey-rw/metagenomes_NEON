# load libraries
library(dplyr)          
library(tidyr)   
library(data.table)
library(ggplot2)
library(ggpubr)

### Specify paths - must change for your own file system ###
# Path to BLASTp or BLASTx results
blast_path <- '/projectnb/talbot-lab-data/zrwerbin/metagenomes_raw/test_pipeline/sunbeam_output/BLAST/blastp_card.csv'
blastp_raw <- data.table::fread(blast_path, header = F, data.table = F)

# Path to metagenome metadata, NEON DP1.10107.001 (included in Github repo) - must download for total read counts
sequence_metadata <- read.csv("https://raw.githubusercontent.com/zoey-rw/metagenomes_NEON/main/metadata_metagenome.csv")
# Path to Antibiotic Resistance Ontology (ARO) - download newest version at https://card.mcmaster.ca/latest/ontology 
aro <- read.csv("https://raw.githubusercontent.com/zoey-rw/metagenomes_NEON/main/ARO/aro.tsv", sep = "\t")

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
	mutate(norm = (gene_count/sampleTotalReadNumber), 
				 siteID = substr(sampleID, 1, 4)) 

### Let's make some plots! ###

# Visualize normalized number of hits per site
hits_plot <- ggplot(merged) + geom_dotplot(binaxis = "y", 
																					 aes(x = siteID, y = norm, fill = siteID, color = siteID), show.legend = F) + 
	ggtitle("Tetracycline resistance gene abundance", "Normalized by total reads per sample") + 
	xlab("NEON site") + ylab("Normalized gene abundance") + theme_bw() + 
	theme(axis.text.x  = element_text(angle=290, hjust=-.05, size=12))

# Visualize number of distinct genes
unique_genes_plot <- ggplot(merged) + geom_point(aes(x = sampleTotalReadNumber, y = distinct_genes, color = siteID)) + 
	ggtitle("Unique tetracycline resistance genes") + 
	xlab("Sequencing depth") + ylab("Number of unique tetracycline-resistance genes") + theme_bw() + 
	theme(axis.text.x  = element_text(angle=290, hjust=-.05, size=12)) + guides(color=guide_legend(title="NEON site"))

ggarrange(hits_plot, unique_genes_plot, labels = c("A", "B"))
