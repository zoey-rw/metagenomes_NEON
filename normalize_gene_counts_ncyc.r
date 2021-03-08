
library(DESeq2)
library(tidyverse)
library(dplyr)

### Specify data path - must change for your own file system ###

# Path to BLASTp or BLASTx results
blast_path <- '/projectnb/talbot-lab-data/zrwerbin/metagenomes_raw/manuscript_test/metagenome_analysis//sunbeam_output/annotation/blastp/blastp_ncyc_allsamples.csv'
blastp_raw <- data.table::fread(blast_path, data.table = F)

# Path to NCycDB key (included in Github repo)
N_pathway_key <- read.csv('https://raw.githubusercontent.com/zoey-rw/metagenomes_NEON/main/NCycDB_pathways.csv', 
													stringsAsFactors = F, strip.white = T)  
# Path to metagenome metadata, NEON DP1.10107.001 (included in Github repo) - must download for total read countss
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
	dplyr::rename(alignment_length = V3, e_val = V4)

# Add N-cycling pathway data
blastp$pathway <- N_pathway_key[match(blastp$gene, N_pathway_key$Gene_family),]$Pathway

# Get counts per gene
gene_counts <- blastp %>% filter(!grepl("Others", pathway) & !is.na(pathway)) %>% 
	group_by(sampleID, gene) %>% add_count(name = "gene_count") %>%
	distinct(sampleID, gene, gene_count, .keep_all = T)

# Make sample names compatible
sequence_metadata$sampleID <- gsub("COMP-DNA1", "comp", sequence_metadata$sampleID)
# Prep count data for DeSeq
countdata <- merge(gene_counts, sequence_metadata, by = "sampleID") %>% 
	mutate(norm = (gene_count/sampleTotalReadNumber)) %>% 
	pivot_wider(id_cols = "gene", names_from = "sampleID", values_from = "gene_count", values_fill=0) %>% 
	tibble::column_to_rownames("gene") 

# prep sample data
coldata <- sequence_metadata %>% filter(sampleID %in% merged$sampleID)
rownames(coldata) <- coldata$sampleID
coldata <- coldata[colnames(countdata),] 
coldata$siteID = substr(coldata$sampleID, 1, 4)
coldata$horizon = ifelse(grepl("-M-", coldata$sampleID), "M", "O")
coldata$X <- NULL

# Check that sample names match up
identical(colnames(countdata), rownames(coldata))

# Add column for "Other" read counts (those that did not match our genes of interest)
other <- coldata$sampleTotalReadNumber - colSums(countdata)
countdata <- rbind(countdata, other)

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

# Merge with pathway and horizon data
merged <- merge(gathered_normalized_counts, N_pathway_key, all.y = F, by.x = "gene", by.y="Gene_family")
coldata$samplename <- gsub("-",".", coldata$sampleID)
merged <- merge(merged, coldata, all.y = F)

# subset to one pathway
merged <- merged %>% filter(grepl("Organic degradation and synthesis", Pathway))
## plot using ggplot2
ncyc <- ggplot(merged) +
	geom_point(aes(x = gene, y = normalized_counts, color = samplename), size = 3) +
	scale_y_log10() +
	#facet_grid(~Pathway, scales = "free", space = "free") +
	xlab("") +
	ylab("") +
	ggtitle("Organic degradation and synthesis genes") +
	theme_bw() +
	theme(text = element_text(size=20),
				axis.text.x = element_text(angle = 45, hjust = 1),
				plot.title = element_text(hjust = 0.5)) + labs(color = "Sample name")
