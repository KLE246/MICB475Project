library(tidyverse)
library(phyloseq)
library(ape)
library(ggsignif)
library(picante)

# Shannon alpha diversity plot
set.seed(1)
phylof <- readRDS("rdata/phyloseq_final.RDS")

gg_richness <- plot_richness(phylof, x = "sex", measures = "Shannon") +
  geom_boxplot()
gg_richness

ggsave(filename = "plot_richness.png", 
       gg_richness)

# combine data and richness
alphadiv <- estimate_richness(phylof, measures = "Shannon")
samp_data <- sample_data(phylof)
samp_data_wdiv <- data.frame(samp_data, alphadiv)


# wilcoxon rank sum 
wilcox.test(Shannon ~ sex, data = samp_data_wdiv, exact = FALSE)
# 0.7824
