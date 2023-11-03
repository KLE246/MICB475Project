library(phyloseq)
library(tidyverse)
library(ape)
library(indicspecies)

metaFP <- "metadata/combined_md.txt"
meta <- read.delim(file=metaFP)

# TODO, fill in with correct files after export, join files for both infant and anemia
otuFP <- "feature-table.txt"
otu <- read.delim(file=otuFP, skip=1, row.names = 1)
taxFP <- "taxonomy.tsv"
tax <- read.delim(file=taxFP)
treeFP <- "tree.nwk"
tree <- read.tree(treeFP)

otu_mat <- as.matrix(otu)
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$sample_name
SAMP <- sample_data(samp_df)

tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix()
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature.ID`
TAX <- tax_table(tax_mat)

mouse_pseq <- phyloseq(OTU, SAMP, TAX, tree)
