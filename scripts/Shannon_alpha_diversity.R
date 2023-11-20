library(tidyverse)
library(phyloseq)
library(ape)
library(ggsignif)
library(picante)
library(microViz)
# BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"), update = FALSE)
# install.packages(
#   "microViz",
#   repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
# )

# Shannon alpha diversity plot
set.seed(1)
phylof <- readRDS("rdata/phyloseq_final.RDS") %>%
  ps_filter((sex == "female" | sex == "male"))

gg_richness <- plot_richness(phylof, x = "sex", measures = "Shannon") +
  geom_boxplot()
gg_richness

ggsave(filename = "plots/plot_richness.png", 
       gg_richness)

# combine data and richness
alphadiv <- estimate_richness(phylof, measures = "Shannon")
samp_data <- sample_data(phylof)
samp_data_wdiv <- data.frame(samp_data, alphadiv)


# wilcoxon rank sum 
wilcox.test(Shannon ~ sex, data = samp_data_wdiv, exact = FALSE)
# 0.7824
# 

# Aim 4 split
cohort_split <- gg_richness +
  facet_grid(~ cohort)
  
ggsave(filename = "plots/plot_richness_cohort_split.png", 
       cohort_split)
