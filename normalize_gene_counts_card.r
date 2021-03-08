
library(DESeq2)
library(tidyverse)
library(dplyr)

### Specify data path - must change for your own file system ###\


# Path to BLASTp or BLASTx results
blast_path <- '/projectnb/talbot-lab-data/zrwerbin/metagenomes_raw/manuscript_test/metagenome_analysis//sunbeam_output/annotation/blastp/blastp_card_allsamples.csv'
blastp_raw <- data.table::fread(blast_path, data.table = F, header = F)

# Path to metagenome metadata, NEON DP1.10107.001 (included in Github repo) - must download for total read countss
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
	dplyr::rename(Accession = V5) 

gene_counts <- blastp_counts %>% 
	group_by(sampleID, blast_hit,Accession) %>% add_count(name = "gene_count") %>%
	distinct(sampleID, blast_hit, Accession, gene_count, .keep_all = T)

# This step is only because we're using legacy NEON data, with wonky sample names.
sequence_metadata$sampleID <- gsub("COMP-DNA1", "comp", sequence_metadata$sampleID)

# Prep count data
countdata <- merge(gene_counts, sequence_metadata, by = "sampleID") %>% 
	pivot_wider(id_cols = "Accession", names_from = "sampleID", values_from = "gene_count", values_fill=0) %>% 
	tibble::column_to_rownames("Accession")

# Prep sample data
coldata <- sequence_metadata %>% filter(sampleID %in% gene_counts$sampleID)
rownames(coldata) <- coldata$sampleID
coldata <- coldata[colnames(countdata),] 
coldata$X <- NULL

# Check that sample names match up
identical(colnames(countdata), rownames(coldata))

# Add column for "Other" read counts (those that did not match our genes of interest)
other <- coldata$sampleTotalReadNumber - colSums(countdata)
countdata <- rbind(countdata, other)
rownames(countdata)[nrow(countdata)] <- "other"

# create DE 'object'
ddsFullCountTable <- DESeqDataSetFromMatrix(
	countData = countdata,
	colData = coldata,
	design = ~ sampleID)
ddsFullCountTable
dds <- ddsFullCountTable

#as.data.frame( colData(ddsFullCountTable) )
sf <- estimateSizeFactors(dds)
normalized_counts <- counts(sf, normalized=TRUE) %>% 
	data.frame() %>%
	rownames_to_column(var="gene")
gathered_normalized_counts <- normalized_counts %>%
	gather(colnames(normalized_counts)[2:3], key = "samplename", value = "normalized_counts")

# Merge with ARO data
merged <- merge(gathered_normalized_counts, aro, all.y = F, by.x = "gene", by.y="Accession")
coldata$samplename <- gsub("-",".", coldata$sampleID)
merged <- merge(merged, coldata, all.y = F)#, by.x = "samplename", by.y="sampleID")

# subset to specific genes
#merged <- merged %>% filter(grepl("tet[[:alpha:]]|tet\\([[:alpha:]]\\)", Name))
merged <- merged %>% filter(grepl("tet[[:alpha:]]", Name))

## plot using ggplot2
card <- ggplot(merged) +
	geom_point(aes(x = Name, y = normalized_counts, color = samplename), size = 3) +
	scale_y_log10() +
	#facet_grid(~Pathway, scales = "free", space = "free") +
	xlab("") +
	ylab("") +
	ggtitle("Selected tetracycline-resistance genes") +
	theme_bw() +
	theme(text = element_text(size=20),
				axis.text.x = element_text(angle = 45, hjust = 1),
				plot.title = element_text(hjust = 0.5)) + labs(color = "Sample name")



library(grid)
library(gridExtra)
library(ggpubr)

figure <- ggarrange(card, ncyc, nrow=2, common.legend = TRUE, legend="top", labels = c("A", "B"))

annotate_figure(figure, 
								left = text_grob("log10 Normalized Counts", 
																 size = 20, rot = 90), 
								bottom = text_grob("Gene name", 
																	 size = 20))


bitmap("Figure5.tiff", height = 10, width = 14, units = 'in', res = 600)


annotate_figure(figure, 
								left = text_grob("log10 Normalized Counts", 
																 size = 20, rot = 90), 
								bottom = text_grob("Gene name", 
																	 size = 20))
# print(grid.arrange(arrangeGrob(ncyc, left = textGrob("a)", x = unit(1, "npc"), 
# 																									 y = unit(.95, "npc"))), 
# 									 arrangeGrob(card, left =textGrob("b)", x = unit(1, "npc"), 
# 									 															 y = unit(.95, "npc"))),
# 									 layout_matrix = rbind(c(1, 1),
# 									 											c(2,2))))

dev.off()

