# Load packages
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(dplyr)
library(DESeq2)
library(GGally)
library(ggplot2)
library(pheatmap)

#### Combined - Sex ####
# Load metadata
metadata_cmb <- read_delim("metadata/combined_md.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE) #combined
filtered_metadata_cmb <- metadata_cmb %>% filter(sex %in% c("female", "male")) #filtered to remove all except female/male in sex column

## KEGG Orthology Pathway 
# Input data
input_file_ko_a <- "picrust2_output_anemia/KO_metagenome_out/pred_metagenome_unstrat.tsv" #anemia
input_file_ko_i <- "infant/picrust2_output_infant/KO_metagenome_out/pred_metagenome_unstrat.tsv" #infant
ko_abundance_a <- read_delim(input_file_ko_a) #anemia KO abundance
ko_abundance_i <- read_delim(input_file_ko_i) #infant KO abundance
ko_abundance_merg <- merge(ko_abundance_a, ko_abundance_i) #merge 
                                                         
# Convert KO abundance to KEGG pathway abundance
kegg_abundance_cmb <- ko2kegg_abundance(data = ko_abundance_merg) #KO to KEGG abundance; combined

# Preparing for DAA using DESeq2 
common_samples <- filtered_metadata_cmb$`#SampleID` %>% intersect(colnames(kegg_abundance_cmb)) #select for common identifiers
kegg_abundance_filtered <- kegg_abundance_cmb[, common_samples] #filter kegg_abundance data

# Perform pathway differential abundance analysis (DAA) using DESeq2 method
daa_results_df_cmb <- pathway_daa(abundance = kegg_abundance_filtered, metadata = filtered_metadata_cmb, group = "sex", daa_method = "DESeq2", select = NULL, reference = NULL) #using filtered metadata (only female/male in sex column)

# Annotate pathway 
daa_annotated_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_df_cmb, ko_to_kegg = TRUE)
  # RESULTS: NO STATISTICAL SIGNIFICANCE

# Annotate pathway without KO to KEGG conversion
daa_annotated_sub_method_results_df <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_df_cmb, ko_to_kegg = FALSE)
  # RESULTS: NO STATISTICAL SIGNIFICANCE

# Generate PCA plot
ko_pca_plot <- pathway_pca(abundance = kegg_abundance_filtered, metadata = filtered_metadata_cmb, group = "sex")


## Metacyc Pathway 
# Input data
input_file_metacyc_a <- "picrust2_output_anemia/EC_metagenome_out/pred_metagenome_unstrat.tsv" #anemia
input_file_metacyc_i <- "infant/picrust2_output_infant/EC_metagenome_out/pred_metagenome_unstrat.tsv" #infant
metacyc_abundance_a <- read_delim(input_file_metacyc_a, delim = "\t") #anemia MetaCyc abundance
metacyc_abundance_i <- read_delim(input_file_metacyc_i, delim = "\t") #infant MetaCyc abundance
metacyc_abundance_merge <- merge(metacyc_abundance_a, metacyc_abundance_i) #merge!

# Preparing for DAA using DESeq2 
common_samples_meta <- intersect(colnames(metacyc_abundance_merge), filtered_metadata_cmb$`#SampleID`) #select for common SampleID's
metacyc_abundance_filtered <- metacyc_abundance_merge[, c("function", common_samples_meta)] #filter metacyc_abundance data

# Perform pathway DAA using DESeq2
metacyc_daa_results_df <- pathway_daa(abundance = metacyc_abundance_filtered %>% column_to_rownames("function"), metadata = filtered_metadata_cmb, group = "sex", daa_method = "DESeq2")

# Annotate MetaCyc Pathway results
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_daa_results_df, ko_to_kegg = FALSE)
daa_results_filtered_sub_df <- subset(metacyc_daa_annotated_results_df, p_values < 0.05)


# Adding annotation to MetaCyc Pathway results
input_file_metacyc_descrip_a <- "picrust2_output/picrust2_output_anemia/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv" #inputting anemia description file
metacyc_abundance_descrip_a <- read_delim(input_file_metacyc_descrip_a, delim = "\t")
input_file_metacyc_descrip_i <- "picrust2_output/picrust2_output_infant/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv" #inputting infant description file
metacyc_abundance_descrip_i <- read_delim(input_file_metacyc_descrip_i, delim = "\t")
metacyc_abundance_descrip_merg <- merge(metacyc_abundance_descrip_a, metacyc_abundance_descrip_i) #merge -- WORKED!
names(metacyc_abundance_descrip_merg)[names(metacyc_abundance_descrip_merg) == "function"] <- "feature" #renaming function into feature to match DESeq2 data
metacyc_daa_annotated_results_df <- merge(metacyc_daa_results_df, metacyc_abundance_descrip_merg[, c("feature", "description")], by = "feature", all.x = TRUE) #annotate

# Generate PCA plot
metacyc_pca_plot <- pathway_pca(abundance = metacyc_abundance_filtered %>% column_to_rownames("function"), metadata = filtered_metadata_cmb, group = "sex")

# Generate pathway error bar plot
metacyc_pathway <- pathway_errorbar(
  abundance = metacyc_abundance_filtered %>% column_to_rownames("function"),
  daa_results_df = daa_results_filtered_sub_df,
  Group = filtered_metadata_cmb$sex,
  ko_to_kegg = FALSE,
  p_values_threshold = 0.06,
  order = "group",
  select = NULL,
  p_value_bar = TRUE,
  colors = NULL,
  x_lab = "description"
)
  # RESULTS: NO STATISTICAL SIGNIFICANCE


## Pathway (EC) Abundance
# Input data
input_file_path_a <- "picrust2_output/picrust2_output_anemia/pathways_out/path_abun_unstrat.tsv" #anemia
input_file_path_i <- "picrust2_output/picrust2_output_infant/pathways_out/path_abun_unstrat.tsv" #infant
path_abundance_a <- read_delim(input_file_path_a, delim = "\t") #anemia path abundance
path_abundance_i <- read_delim(input_file_path_i, delim = "\t") #infant path abundance
path_abundance_merge <- merge(path_abundance_a, path_abundance_i) #merge!

# Preparing for DAA using DESeq2 
common_samples_path <- intersect(colnames(path_abundance_merge), filtered_metadata_cmb$`#SampleID`) #select for common SampleID's
path_abundance_filtered <- path_abundance_merge[, c("pathway", common_samples_path)] #filter metacyc_abundance data

# Perform pathway DAA using DESeq2
path_daa_results_df <- pathway_daa(abundance = path_abundance_filtered %>% column_to_rownames("pathway"), metadata = filtered_metadata_cmb, group = "sex", daa_method = "DESeq2")

# Adding annotation to MetaCyc Pathway results
input_file_path_descrip_a <- "picrust2_output/picrust2_output_anemia/pathways_out/path_abun_unstrat_descrip.tsv" #inputting anemia description file
path_abundance_descrip_a <- read_delim(input_file_path_descrip_a, delim = "\t")
input_file_path_descrip_i <- "picrust2_output/picrust2_output_infant/pathways_out/path_abun_unstrat_descrip.tsv" #inputting infant description file
path_abundance_descrip_i <- read_delim(input_file_path_descrip_i, delim = "\t")
path_abundance_descrip_merg <- merge(path_abundance_descrip_a, path_abundance_descrip_i) #merge -- WORKED!
names(path_abundance_descrip_merg)[names(path_abundance_descrip_merg) == "pathway"] <- "feature" #renaming function into feature to match DESeq2 data
path_daa_annotated_results_df <- merge(path_daa_results_df, path_abundance_descrip_merg[, c("feature", "description")], by = "feature", all.x = TRUE) #annotate

# Generate PCA plot
path_pca_plot <- pathway_pca(abundance = path_abundance_filtered %>% column_to_rownames("pathway"), metadata = filtered_metadata_cmb, group = "sex")

# Generate pathway error bar plot
path_pathway <- pathway_errorbar(
  abundance = path_abundance_filtered %>% column_to_rownames("pathway"),
  daa_results_df = path_daa_annotated_results_df,
  Group = filtered_metadata_cmb$sex,
  ko_to_kegg = FALSE,
  p_values_threshold = 0.06,
  order = "group",
  select = NULL,
  p_value_bar = TRUE,
  colors = NULL,
  x_lab = "description"
)


#### Combined - Location ####

## KEGG Orthology Pathway 
# Input data
ko_abundance_merg <- merge(ko_abundance_a, ko_abundance_i) #merge 

# Convert KO abundance to KEGG pathway abundance
kegg_abundance_cmb_location <- ko2kegg_abundance(data = ko_abundance_merg) #KO to KEGG abundance; combined

# Preparing for DAA using DESeq2 
common_samples_location <- metadata_cmb$`#SampleID` %>% intersect(colnames(kegg_abundance_cmb_location)) #select for common identifiers
kegg_abundance_filtered_location <- kegg_abundance_cmb_location[, common_samples_location] #filter kegg_abundance data

# Perform pathway differential abundance analysis (DAA) using DESeq2 method
daa_results_df_cmb_location <- pathway_daa(abundance = kegg_abundance_filtered_location, metadata = metadata_cmb, group = "cohort", daa_method = "DESeq2", select = NULL, reference = NULL) #using combined metadata (region)

# Annotate pathway 
subset_features_ko_location <- head(daa_results_df_cmb_location$feature, 30) #subsetting top 175 as dataset too large
daa_results_subset_ko_location <- daa_results_df_cmb_location[daa_results_df_cmb_location$feature %in% subset_features_ko_location, ]
daa_annotated_results_df_location <- pathway_annotation(pathway = "KO", daa_results_df = daa_results_subset_ko_location, ko_to_kegg = TRUE) #annotating

# Generate pathway error bar plot
p_ko_location <- pathway_errorbar(abundance = kegg_abundance_filtered_location, 
                          daa_results_df = daa_annotated_results_df_location,
                          Group = metadata_cmb$cohort,
                          ko_to_kegg = TRUE,
                          p_values_threshold = 0.05,
                          order = "pathway_class",
                          select = NULL, 
                          p_value_bar = TRUE, 
                          colors = NULL, 
                          x_lab = "pathway_name")

# Generate PCA plot
ko_pca_plot_location <- pathway_pca(abundance = kegg_abundance_filtered_location, metadata = metadata_cmb, group = "cohort")


## Metacyc Pathway 
metacyc_abundance_merge_location <- merge(metacyc_abundance_a, metacyc_abundance_i) #merge!

# Preparing for DAA using DESeq2 
common_samples_meta_location <- intersect(colnames(metacyc_abundance_merge_location), metadata_cmb$`#SampleID`) #select for common SampleID's
metacyc_abundance_filtered_location <- metacyc_abundance_merge_location[, c("function", common_samples_meta_location)] #filter metacyc_abundance data

# Perform pathway DAA using DESeq2
metacyc_daa_results_df_location <- pathway_daa(abundance = metacyc_abundance_filtered_location %>% column_to_rownames("function"), metadata = metadata_cmb, group = "cohort", daa_method = "DESeq2")

# Adding annotation to MetaCyc Pathway results
metacyc_daa_annotated_results_df_location <- merge(metacyc_daa_results_df_location, metacyc_abundance_descrip_merg[, c("feature", "description")], by = "feature", all.x = TRUE) #annotate

# Generate PCA plot
metacyc_pca_plot_location <- pathway_pca(abundance = metacyc_abundance_filtered_location %>% column_to_rownames("function"), metadata = metadata_cmb, group = "cohort")

# Generate pathway error bar plot
selected_p_0.05_metacyc_location <- feature_with_p_0.05_metacyc_location %>% slice_head(n = 30)
meta_pathway_location <- pathway_errorbar(
  abundance = metacyc_abundance_filtered_location %>% column_to_rownames("function"),
  daa_results_df = selected_p_0.05_metacyc_location,
  Group = metadata_cmb$cohort,
  ko_to_kegg = FALSE,
  p_values_threshold = 0.05,
  order = "group",
  select = NULL,
  p_value_bar = TRUE,
  colors = NULL,
  x_lab = "description"
)


## Pathway (EC) Abundance
path_abundance_merge_location <- merge(path_abundance_a, path_abundance_i) #merge!

# Preparing for DAA using DESeq2 
common_samples_path_location <- intersect(colnames(path_abundance_merge_location), metadata_cmb$`#SampleID`) #select for common SampleID's
path_abundance_filtered_location <- path_abundance_merge_location[, c("pathway", common_samples_path_location)] #filter metacyc_abundance data

# Perform pathway DAA using DESeq2
path_daa_results_df_location <- pathway_daa(abundance = path_abundance_filtered_location %>% column_to_rownames("pathway"), metadata = metadata_cmb, group = "cohort", daa_method = "DESeq2")

# Adding annotation to MetaCyc Pathway results
path_daa_annotated_results_df_location <- merge(path_daa_results_df_location, path_abundance_descrip_merg[, c("feature", "description")], by = "feature", all.x = TRUE) #annotate

# Generate PCA plot
path_pca_plot_location <- pathway_pca(abundance = path_abundance_filtered_location %>% column_to_rownames("pathway"), metadata = metadata_cmb, group = "cohort")

# Generate pathway error bar plot
feature_with_p_0.05_path_location <- path_daa_annotated_results_df_location %>% 
  filter(p_adjust < 0.05) # Filter features with p < 0.05
selected_p_0.05_path_location <- feature_with_p_0.05_path_location %>% slice_head(n = 30)

path_pathway_location <- pathway_errorbar(
  abundance = path_abundance_filtered_location %>% column_to_rownames("pathway"),
  daa_results_df = selected_p_0.05_path_location,
  Group = metadata_cmb$cohort,
  ko_to_kegg = FALSE,
  p_values_threshold = 0.05,
  order = "group",
  select = NULL,
  p_value_bar = TRUE,
  colors = NULL,
  x_lab = "description"
)
