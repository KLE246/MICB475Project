library(DESeq2)
library(phyloseq)
library(tidyverse)

# Read in phyloseq object and convert to deseq, where the control is female infants
combined_phyloseq <- readRDS("rdata/phyloseq_final.rds")

anemia_phyloseq <-  subset_samples(combined_phyloseq, cohort == "anemia")


anemia_phyloseq <- transform_sample_counts(anemia_phyloseq, function(x) x+1)
anemia_deseq <- phyloseq_to_deseq2(anemia_phyloseq, ~ sex)
anemia_DESeq_out <- DESeq(anemia_deseq)

anemia_res <- results(anemia_DESeq_out, tidy=TRUE)

anemia_sigASVs <- anemia_res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)%>%
  cbind()

# Get a vector of ASV names
anemia_sigASVs_vec <- anemia_sigASVs %>%
  pull(ASV)

# Prune phyloseq file
anemia_pruned <- prune_taxa(anemia_sigASVs_vec,anemia_phyloseq)
anemia_sigASVs_tax <- tax_table(anemia_pruned) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(anemia_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))%>%
  filter(!grepl("NA.",Genus), !is.na(Genus)) %>%
  mutate(cohort = "anemia")

### Repeat for infant cohort
infant_phyloseq = subset_samples(combined_phyloseq, cohort == "infant")

infant_phyloseq <- transform_sample_counts(infant_phyloseq, function(x) x+1)
infant_deseq <- phyloseq_to_deseq2(infant_phyloseq, ~ sex)
infant_DESeq_out <- DESeq(infant_deseq)

infant_res <- results(infant_DESeq_out, tidy=TRUE)

infant_sigASVs <- infant_res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>2) %>%
  dplyr::rename(ASV=row)%>%
  cbind()

# Get a vector of ASV names
infant_sigASVs_vec <- infant_sigASVs %>%
  pull(ASV)

# Prune phyloseq file
infant_pruned <- prune_taxa(infant_sigASVs_vec,infant_phyloseq)
infant_sigASVs_tax <- tax_table(infant_pruned) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(infant_sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))%>%
  filter(!grepl("NA.",Genus), !is.na(Genus)) %>%
  mutate(cohort = "infant")


combined_results <- rbind(anemia_sigASVs_tax, infant_sigASVs_tax)

bar_plot <- ggplot(combined_results) +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill=log2FoldChange>0), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(y = "male/female change", title="Change in Genus abundance between male and female")+
  theme(legend.position="none")+
  facet_wrap(~cohort)

bar_plot
ggsave(filename="plots/separated_cohort_bar_plot.png",bar_plot)

similar_bar_plot <- combined_results %>% 
  filter(Genus %in% (split(combined_results$Genus, combined_results$cohort) %>%
                       reduce(intersect))) %>%
  ggplot() +
  geom_bar(aes(x=Genus, y=log2FoldChange, fill=log2FoldChange>0), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  labs(y = "male/female change", title="Change in Genus abundance between male and female")+
  theme(legend.position="none")+
  facet_wrap(~cohort)

similar_bar_plot
ggsave(filename="plots/similar_cohort_bar_plot.png",similar_bar_plot)


