library(phyloseq)
library(tidyverse)
library(ape)
library(vegan)
library(picante)

# Read in phyloseq object and convert to deseq, where the control is female infants

combined_phyloseq <- readRDS("phyloseq_final.rds")
combined_phyloseq = subset_samples(combined_phyloseq, sex != "Missing: Not collected")

# Rarefraction at depth of 1600
mpt_rare <- rarefy_even_depth(combined_phyloseq, rngseed = 1, sample.size = 1600)
samp_dat_wdiv <- data.frame(sample_data(mpt_rare), estimate_richness(mpt_rare))


#### Beta diversity #####
bc_dm <- distance(mpt_rare, method="bray")
# check which methods you can specify
?distance

pcoa_bc <- ordinate(mpt_rare, method="PCoA", distance=bc_dm)

gg_pcoa <- plot_ordination(mpt_rare, pcoa_bc, color = "sex") +
  labs(pch="cohort", col = "sex")
gg_pcoa

gg_pcoa_cohort <- plot_ordination(mpt_rare, pcoa_bc, color = "cohort") +
  labs(pch="cohort", col = "sex")
gg_pcoa_cohort

# preform statistical testing using permanova comparing sex and location 
# as a factor
permanova <- adonis2(bc_dm ~ sex*cohort, data=samp_dat_wdiv)
permanova

# re-plot the above PCoA with ellipses to show a significant difference 

plot_ordination(mpt_rare, pcoa_bc, color = "sex") +
  stat_ellipse(type = "norm")

plot_ordination(mpt_rare, pcoa_bc, color = "cohort") +
  stat_ellipse(type = "norm")

#### Taxonomy bar plots ####

# Plot bar plot of taxonomy
plot_bar(mpt_rare, fill="Phylum") 

# Convert to relative abundance
mpt_RA <- transform_sample_counts(mpt_rare, function(x) x/sum(x))

mpt_phylum <- tax_glom(mpt_RA, taxrank = "Phylum", NArm=FALSE)

gg_taxa_cohort <- plot_bar(mpt_phylum, fill="Phylum") + 
  facet_wrap(.~cohort, scales = "free_x")
gg_taxa_cohort

gg_taxa_sex <- plot_bar(mpt_phylum, fill="Phylum") + 
  facet_wrap(.~sex, scales = "free_x")
gg_taxa_sex
