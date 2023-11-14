library(DESeq2)
library(phyloseq)
library(tidyverse)

# Read in phyloseq object and convert to deseq, where the control is female infants
combined_phyloseq <- readRDS("rdata/phyloseq_final.rds")

combined_plus1 <- transform_sample_counts(combined_phyloseq, function(x) x+1)
combined_deseq <- phyloseq_to_deseq2(combined_plus1, ~ sex)
DESeq_out <- DESeq(combined_deseq)
res <- results(DESeq_out, tidy=TRUE)

## Volcano plot
# Effect size VS significance
ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
vol_plot <- res %>%
  mutate(significant = padj<0.05 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant)) +
  labs(x = "male/female change", title="Significant abundance differences between male and female")

vol_plot
ggsave(filename="plots/vol_plot_deseq.png",vol_plot)

## Bar plot
# To get table of results
sigASVs <- res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)

# Get a vector of ASV names
sigASVs_vec <- sigASVs %>%
  pull(ASV)

# Prune phyloseq file
combined_DESeq <- prune_taxa(sigASVs_vec,combined_phyloseq)
sigASVs <- tax_table(combined_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))%>%
  filter(!grepl("NA",Genus), !is.na(Genus))

bar_plot <- ggplot(sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill=log2FoldChange>0), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(y = "male/female change", title="Change in Genus abundance between male and female")+
  theme(legend.position="none")

bar_plot
ggsave(filename="plots/bar_plot_deseq.png",bar_plot)
