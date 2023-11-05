library(tidyverse)
library(phyloseq)
library(ape)
library(vegan)

# metadata
meta <- read_tsv(file = "metadata/combined_md.txt") %>%
  rename(SampleID = `#SampleID`)

samp_df <- as.data.frame(meta[, -1])
rownames(samp_df) <- meta$'SampleID'
SAMP <- sample_data(samp_df)

# OTU 
otua <- read_delim(file = "anemia/anemia_export/table_export/feature-table.txt", delim = "\t", skip = 1)
otui <- read_delim(file = "infant/infant_export/table_export/feature-tablei.txt", delim = "\t", skip = 1)

otu_mata <- as.matrix(otua[,-1])
rownames(otu_mata) <- otua$'#OTU ID'
OTUA <- otu_table(otu_mata, taxa_are_rows = TRUE)
otu_mati <- as.matrix(otui[,-1])
rownames(otu_mati) <- otui$'#OTU ID'
OTUI <- otu_table(otu_mati, taxa_are_rows = TRUE)

# load taxonomy
taxa <- read_delim(file = "anemia/anemia_export/taxonomy_export/taxonomy.tsv", delim = "\t")
taxi <- read_delim(file = "infant/infant_export/taxonomy_export/taxonomy.tsv", delim = "\t")

# taxonomy anemia
tax_mata <- taxa %>% select(-Confidence) %>%
  separate(col = Taxon, sep = ";"
           , into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  as.matrix()
tax_mata <- tax_mata[, -1]
rownames(tax_mata) <- taxa$'Feature ID'
TAXA <- tax_table(tax_mata)

# taxonomy infant
tax_mati <- taxi %>% select(-Confidence) %>%
  separate(col = Taxon, sep = ";"
           , into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  as.matrix()
tax_mati <- tax_mati[, -1]
rownames(tax_mati) <- taxi$'Feature ID'
TAXI <- tax_table(tax_mati)

# tree 
treea <- read.tree("anemia/anemia_export/rooted_tree_export/tree.nwk")
treei <- read.tree("infant/infant_export/rooted_tree_export/tree.nwk")

treef <- bind.tree(treea, treei, where = "root", position = 0, interactive = FALSE)
duplicate <- duplicated(treef$tip.label)
treei <- treef$tip.label
for (i in which(duplicate)) {
  count <- sum(duplicate[1:i])
  treei[i] <- paste0(treef$tip.label[i], "-", count)
}
treef$tip.label <- treei

# phyloseq objects 
phyloanemia <- phyloseq(SAMP, OTUA, TAXA, treef)
phyloinfant <- phyloseq(SAMP, OTUI, TAXI)

phylofinal <- merge_phyloseq(phyloanemia, phyloinfant)

saveRDS(phylofinal, file = "rdata/phyloseq_final.rds")
