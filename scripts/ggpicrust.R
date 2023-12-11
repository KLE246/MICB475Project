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

#### Combined - Sex ####
# Load metadata
metadata_cmb <- read_delim("metadata/combined_md.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE) #combined
filtered_metadata_cmb <- metadata_cmb %>% filter(sex %in% c("female", "male")) #filtered to remove all except female/male in sex column

## KEGG Orthology Pathway 
# Input data
input_file_ko_a <- "anemia/picrust2_output_anemia/KO_metagenome_out/pred_metagenome_unstrat.tsv" #anemia
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
ko_pca_plot <- pathway_pca(
  abundance = kegg_abundance_filtered, 
  metadata = filtered_metadata_cmb, 
  group = "sex")

/// 

## Enzyme Commission Pathway 
# Input data
input_file_EC_a <- "picrust2_output/picrust2_output_anemia/EC_metagenome_out/pred_metagenome_unstrat.tsv" #anemia
input_file_EC_i <- "picrust2_output/picrust2_output_infant/EC_metagenome_out/pred_metagenome_unstrat.tsv" #infant
EC_abundance_a <- read_delim(input_file_EC_a, delim = "\t") #anemia MetaCyc abundance
EC_abundance_i <- read_delim(input_file_EC_i, delim = "\t") #infant MetaCyc abundance
EC_abundance_merge <- merge(EC_abundance_a, EC_abundance_i) #merge!

# Preparing for DAA using DESeq2 
common_samples_EC <- intersect(colnames(EC_abundance_merge), filtered_metadata_cmb$`#SampleID`) #select for common SampleID's
EC_abundance_filtered <- EC_abundance_merge[, c("function", common_samples_EC)] #filter metacyc_abundance data

# Perform pathway DAA using DESeq2
EC_daa_results_df <- pathway_daa(abundance = EC_abundance_filtered %>% column_to_rownames("function"), metadata = filtered_metadata_cmb, group = "sex", daa_method = "DESeq2")

# Annotate EC Pathway results
EC_daa_annotated_results_df <- pathway_annotation(pathway = "EC", daa_results_df = EC_daa_results_df, ko_to_kegg = FALSE)
daa_results_filtered_sub_df <- subset(EC_daa_annotated_results_df, p_values < 0.05)

# Adding annotation to EC Pathway results
input_file_EC_descrip_a <- "picrust2_output/picrust2_output_anemia/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv" #inputting anemia description file
EC_abundance_descrip_a <- read_delim(input_file_EC_descrip_a, delim = "\t")
input_file_EC_descrip_i <- "picrust2_output/picrust2_output_infant/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv" #inputting infant description file
EC_abundance_descrip_i <- read_delim(input_file_EC_descrip_i, delim = "\t")
EC_abundance_descrip_merg <- merge(EC_abundance_descrip_a, EC_abundance_descrip_i) #merge -- WORKED!
names(EC_abundance_descrip_merg)[names(EC_abundance_descrip_merg) == "function"] <- "feature" #renaming function into feature to match DESeq2 data
EC_daa_annotated_results_df <- merge(EC_daa_results_df, EC_abundance_descrip_merg[, c("feature", "description")], by = "feature", all.x = TRUE) #annotate

# Generate PCA plot
EC_pca_plot <- pathway_pca(abundance = EC_abundance_filtered %>% column_to_rownames("function"), metadata = filtered_metadata_cmb, group = "sex")

# Generate pathway error bar plot
EC_pathway <- pathway_errorbar(
  abundance = EC_abundance_filtered %>% column_to_rownames("function"),
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
# NO STATISTICAL SIGNIFICANCE

///

## Metabolic Pathway Abundance
# Input data
input_file_metacyc_a <- "picrust2_output/picrust2_output_anemia/pathways_out/path_abun_unstrat.tsv" #anemia
input_file_metacyc_i <- "picrust2_output/picrust2_output_infant/pathways_out/path_abun_unstrat.tsv" #infant
metacyc_abundance_a <- read_delim(input_file_metacyc_a, delim = "\t") #anemia path abundance
metacyc_abundance_i <- read_delim(input_file_metacyc_i, delim = "\t") #infant path abundance
metacyc_abundance_merge <- merge(metacyc_abundance_a, metacyc_abundance_i) #merge!

# Preparing for DAA using DESeq2 
common_samples_metacyc <- intersect(colnames(metacyc_abundance_merge), filtered_metadata_cmb$`#SampleID`) #select for common SampleID's
metacyc_abundance_filtered <- metacyc_abundance_merge[, c("pathway", common_samples_metacyc)] #filter metacyc_abundance data

# Perform pathway DAA using DESeq2
metacyc_daa_results_df <- pathway_daa(abundance = metacyc_abundance_filtered %>% column_to_rownames("pathway"), metadata = filtered_metadata_cmb, group = "sex", daa_method = "DESeq2")

# Adding annotation to EC Pathway results
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_daa_results_df, ko_to_kegg = FALSE)


# Generate PCA plot
metacyc_pca_plot <- pathway_pca(abundance = metacyc_abundance_filtered %>% column_to_rownames("pathway"), metadata = filtered_metadata_cmb, group = "sex")

# Generate pathway error bar plot
metacyc_pathway <- pathway_errorbar(
  abundance = metacyc_abundance_filtered %>% column_to_rownames("pathway"),
  daa_results_df = metacyc_daa_annotated_results_df,
  Group = filtered_metadata_cmb$sex,
  ko_to_kegg = FALSE,
  p_values_threshold = 0.05,
  order = "group",
  select = NULL,
  p_value_bar = TRUE,
  colors = NULL,
  x_lab = "description"
)

///

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

///

## Enzyme Commission Pathway 
EC_abundance_merge_location <- merge(EC_abundance_a, EC_abundance_i) #merge!

# Preparing for DAA using DESeq2 
common_samples_EC_location <- intersect(colnames(EC_abundance_merge_location), metadata_cmb$`#SampleID`) #select for common SampleID's
EC_abundance_filtered_location <- EC_abundance_merge_location[, c("function", common_samples_EC_location)] #filter metacyc_abundance data

# Perform pathway DAA using DESeq2
EC_daa_results_df_location <- pathway_daa(abundance = EC_abundance_filtered_location %>% column_to_rownames("function"), metadata = metadata_cmb, group = "cohort", daa_method = "DESeq2")

# Adding annotation to MetaCyc Pathway results
EC_daa_annotated_results_df_location <- merge(EC_daa_results_df_location, EC_abundance_descrip_merg[, c("feature", "description")], by = "feature", all.x = TRUE) #annotate

# Generate PCA plot
EC_pca_plot_location <- pathway_pca(abundance = EC_abundance_filtered_location %>% column_to_rownames("function"), metadata = metadata_cmb, group = "cohort")


# Generate pathway error bar plot
selected_p_0.05_EC_location <- EC_daa_annotated_results_df_location %>% slice_head(n = 30)
EC_pathway_location <- pathway_errorbar(
  abundance = EC_abundance_filtered_location %>% column_to_rownames("function"),
  daa_results_df = selected_p_0.05_EC_location,
  Group = metadata_cmb$cohort,
  ko_to_kegg = FALSE,
  p_values_threshold = 0.05,
  order = "group",
  select = NULL,
  p_value_bar = TRUE,
  colors = NULL,
  x_lab = "description"
)

///

## Metabolic Pathway Abundance
metacyc_abundance_merge_location <- merge(metacyc_abundance_a, metacyc_abundance_i) #merge!

# Preparing for DAA using DESeq2 
common_samples_metacyc_location <- intersect(colnames(metacyc_abundance_merge_location), metadata_cmb$`#SampleID`) #select for common SampleID's
metacyc_abundance_filtered_location <- metacyc_abundance_merge_location[, c("pathway", common_samples_metacyc_location)] #filter metacyc_abundance data

# Perform pathway DAA using DESeq2
metacyc_daa_results_df_location <- pathway_daa(abundance = metacyc_abundance_filtered_location %>% column_to_rownames("pathway"), metadata = metadata_cmb, group = "cohort", daa_method = "DESeq2")

# Adding annotation to MetaCyc Pathway results
metacyc_daa_annotated_results_df_location <- pathway_annotation(pathway = "MetaCyc", daa_results_df = metacyc_daa_results_df_location, ko_to_kegg = FALSE)

# Generate PCA plot
metacyc_pca_plot_location <- pathway_pca(abundance = metacyc_abundance_filtered_location %>% column_to_rownames("pathway"), metadata = metadata_cmb, group = "cohort")
metacyc_pca_plot_location_new <- pathway_pca(abundance = metacyc_abundance_filtered_location %>% column_to_rownames("pathway"), metadata = metadata_name_change, group = "cohort") #Uses new cohort names

# Generate pathway error bar plot
feature_with_p_0.05_metacyc_location <- metacyc_daa_annotated_results_df_location %>% 
  filter(p_adjust < 0.05) # Filter features with p < 0.05
selected_p_0.05_metacyc_location <- feature_with_p_0.05_metacyc_location %>% slice_head(n = 30)

metacyc_pathway_location <- pathway_errorbar(
  abundance = metacyc_abundance_filtered_location %>% column_to_rownames("pathway"),
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

# Organizing - big shout out to Avril for providing the code to organize by log2 Fold change
pathway_errorbar_onesigresult = function (abundance, daa_results_df, Group, ko_to_kegg = FALSE, 
                                          wrap_label = F, wraplength=20, fc_cutoff = 0, order_by_log = F,
                                          p_values_threshold = 0.05, order = "group", select = NULL, 
                                          p_value_bar = TRUE, colors = NULL, x_lab = NULL){
  
  # abundance = abundance_data_first_high %>% column_to_rownames("pathway");
  # daa_results_df = metacyc_daa_annotated_results_df; order_by_log = T;
  # Group = metadata_nicole_high$feed; wrap_label = F;wraplength=20; fc_cutoff = 0;
  # p_values_threshold = 0.001; order = "group"; select = NULL; ko_to_kegg = F; p_value_bar = TRUE;
  # colors = NULL; x_lab = "description"
  
  require(stringr)
  
  ##################### Calculate log 2 fold changes first so we can filter them
  errorbar_abundance_mat <- as.matrix(abundance)
  pseudo = min(abundance[abundance>0])/2
  relative_abundance_mat <- apply(t(errorbar_abundance_mat), 
                                  1, function(x) (x+pseudo)/sum(x)) # Applying a pseudocount of 1
  
  sub_relative_abundance_mat <- relative_abundance_mat[rownames(relative_abundance_mat) %in% 
                                                         (daa_results_df %>% filter(!is.na(p_adjust),p_adjust<p_values_threshold) %>% pull(feature)), ,drop=F]
  
  error_bar_matrix <- cbind(sample = colnames(sub_relative_abundance_mat), 
                            group = Group, 
                            as.data.frame(t(sub_relative_abundance_mat)))
  error_bar_df <- as.data.frame(error_bar_matrix)
  error_bar_df$group <- factor(Group, levels = levels(as.factor(Group)))
  error_bar_pivot_longer_df <- tidyr::pivot_longer(error_bar_df, 
                                                   -c(sample, group)) %>% as_tibble()
  error_bar_pivot_longer_tibble_summarised <- error_bar_pivot_longer_df %>% 
    group_by(name, group) %>% summarise(mean = mean(value)) %>% 
    ungroup() %>% 
    group_by(name) %>% arrange(group) %>% 
    summarise(log2FC = log2(mean[1L]/mean[2L]))
  
  passes_fc = error_bar_pivot_longer_tibble_summarised %>% 
    filter(abs(log2FC)>fc_cutoff) %>% pull(name)
  
  log_order = error_bar_pivot_longer_tibble_summarised %>% 
    arrange(log2FC) %>% pull(name)
  
  #####################
  
  missing_pathways <- daa_results_df[is.na(daa_results_df$pathway_name), 
                                     "feature"]
  if (length(missing_pathways) > 0) {
    message("The following pathways are missing annotations and have been excluded: ", 
            paste(missing_pathways, collapse = ", "))
    message("You can use the 'pathway_annotation' function to add annotations for these pathways.")
  }
  column_names <- colnames(daa_results_df)
  group_columns <- grepl("^group", column_names)
  if (sum(group_columns) > 2) {
    p_value_bar <- FALSE
    message("There are more than two 'group' columns in the 'daa_results_df' data frame. As a result, it is not possible to compute the log2 fold values. The 'p_value_bar' has been automatically set to FALSE.")
  }
  daa_results_df <- daa_results_df[!is.na(daa_results_df[, x_lab]), ]
  if (is.null(x_lab)) {
    if (ko_to_kegg == TRUE) {
      x_lab <- "pathway_name"
    } else {
      x_lab <- "description"
    }
    if (is.null(daa_results_df$pathway_name) & is.null(daa_results_df$description)) {
      message("Please utilize the 'pathway_annotation' function to annotate the 'daa_results_df' data frame.")
    }
  }
  if (!(x_lab %in% colnames(daa_results_df))) {
    message("The 'x_lab' you defined does not exist as a column in the 'daa_results_df' data frame.")
  }
  if (nlevels(factor(daa_results_df$method)) != 1) {
    message("The 'method' column in the 'daa_results_df' data frame contains more than one method. Please filter it to contain only one method.")
  }
  if (nlevels(factor(daa_results_df$group1)) != 1 || nlevels(factor(daa_results_df$group2)) != 
      1) {
    message("The 'group1' or 'group2' column in the 'daa_results_df' data frame contains more than one group. Please filter each to contain only one group.")
  }
  if (is.null(colors)) {
    colors <- c("#d93c3e", "#3685bc", "#6faa3e", "#e8a825", 
                "#c973e6", "#ee6b3d", "#2db0a7", "#f25292")[1:nlevels(as.factor(Group))]
  }
  errorbar_abundance_mat <- as.matrix(abundance)
  daa_results_filtered_df <- daa_results_df[daa_results_df$p_adjust < p_values_threshold & 
                                              daa_results_df$feature %in% passes_fc, ]
  if (!is.null(select)) {
    daa_results_filtered_sub_df <- daa_results_filtered_df[daa_results_filtered_df$feature %in% 
                                                             select, ]
  } else {
    daa_results_filtered_sub_df <- daa_results_filtered_df
  }
  if (nrow(daa_results_filtered_sub_df) > 30) {
    message(paste0("The number of features with statistical significance exceeds 30, leading to suboptimal visualization. ", 
                   "Please use 'select' to reduce the number of features.\n", 
                   "Currently, you have these features: ", paste(paste0("\"", 
                                                                        daa_results_filtered_sub_df$feature, "\""), collapse = ", "), 
                   ".\n", "You can find the statistically significant features with the following command:\n", 
                   "daa_results_df %>% filter(p_adjust < 0.05) %>% select(c(\"feature\",\"p_adjust\"))"))
    stop()
  }
  if (nrow(daa_results_filtered_sub_df) == 0) {
    stop("Visualization with 'pathway_errorbar' cannot be performed because there are no features with statistical significance. ", 
         "For possible solutions, please check the FAQ section of the tutorial.")
  }
  relative_abundance_mat <- apply(t(errorbar_abundance_mat), 
                                  1, function(x) x/sum(x))
  
  sub_relative_abundance_mat <- relative_abundance_mat[rownames(relative_abundance_mat) %in% 
                                                         daa_results_filtered_sub_df$feature, ,drop=F]
  
  error_bar_matrix <- cbind(sample = colnames(sub_relative_abundance_mat), 
                            group = Group, t(sub_relative_abundance_mat))
  error_bar_df <- as.data.frame(error_bar_matrix)
  error_bar_df$group <- factor(Group, levels = levels(as.factor(Group)))
  error_bar_pivot_longer_df <- tidyr::pivot_longer(error_bar_df, 
                                                   -c(sample, group))
  error_bar_pivot_longer_tibble <- mutate(error_bar_pivot_longer_df, 
                                          group = as.factor(group))
  error_bar_pivot_longer_tibble$sample <- factor(error_bar_pivot_longer_tibble$sample)
  error_bar_pivot_longer_tibble$name <- factor(error_bar_pivot_longer_tibble$name)
  error_bar_pivot_longer_tibble$value <- as.numeric(error_bar_pivot_longer_tibble$value)
  error_bar_pivot_longer_tibble_summarised <- error_bar_pivot_longer_tibble %>% 
    group_by(name, group) %>% summarise(mean = mean(value), 
                                        sd = stats::sd(value))
  error_bar_pivot_longer_tibble_summarised <- error_bar_pivot_longer_tibble_summarised %>% 
    mutate(group2 = "nonsense")
  switch(order, p_values = {
    order <- order(daa_results_filtered_sub_df$p_adjust)
  }, name = {
    order <- order(daa_results_filtered_sub_df$feature)
  }, group = {
    daa_results_filtered_sub_df$pro <- 1
    for (i in levels(error_bar_pivot_longer_tibble_summarised$name)) {
      error_bar_pivot_longer_tibble_summarised_sub <- error_bar_pivot_longer_tibble_summarised[error_bar_pivot_longer_tibble_summarised$name == i, ]
      pro_group <- error_bar_pivot_longer_tibble_summarised_sub[error_bar_pivot_longer_tibble_summarised_sub$mean == 
                                                                  max(error_bar_pivot_longer_tibble_summarised_sub$mean),]$group
      pro_group <- as.vector(pro_group)
      daa_results_filtered_sub_df[daa_results_filtered_sub_df$feature == i, ]$pro <- pro_group
    }
    order <- order(daa_results_filtered_sub_df$pro, daa_results_filtered_sub_df$p_adjust)
  }, pathway_class = {
    if (!"pathway_class" %in% colnames(daa_results_filtered_sub_df)) {
      stop("The 'pathway_class' column is missing in the 'daa_results_filtered_sub_df' data frame. ", 
           "Please use the 'pathway_annotation' function to annotate the 'pathway_daa' results.")
    }
    order <- order(daa_results_filtered_sub_df$pathway_class, 
                   daa_results_filtered_sub_df$p_adjust)
  }, {
    order <- order
  })
  daa_results_filtered_sub_df <- daa_results_filtered_sub_df[order, ]
  error_bar_pivot_longer_tibble_summarised_ordered <- data.frame(name = NULL, 
                                                                 group = NULL, mean = NULL, sd = NULL)
  for (i in daa_results_filtered_sub_df$feature) {
    error_bar_pivot_longer_tibble_summarised_ordered <- rbind(error_bar_pivot_longer_tibble_summarised_ordered, 
                                                              error_bar_pivot_longer_tibble_summarised[error_bar_pivot_longer_tibble_summarised$name == i, ])
  }
  if (ko_to_kegg == FALSE) {
    error_bar_pivot_longer_tibble_summarised_ordered[, x_lab] <- rep(daa_results_filtered_sub_df[, x_lab], each = length(levels(factor(error_bar_pivot_longer_tibble_summarised_ordered$group))))
  }
  if (ko_to_kegg == TRUE) {
    error_bar_pivot_longer_tibble_summarised_ordered$pathway_class <- rep(daa_results_filtered_sub_df$pathway_class, 
                                                                          each = length(levels(factor(error_bar_pivot_longer_tibble_summarised_ordered$group))))
  }
  
  if(order_by_log == F){
    error_bar_pivot_longer_tibble_summarised_ordered$name <- factor(error_bar_pivot_longer_tibble_summarised_ordered$name, 
                                                                    levels = rev(daa_results_filtered_sub_df$feature))
  } else {
    error_bar_pivot_longer_tibble_summarised_ordered$name <- factor(error_bar_pivot_longer_tibble_summarised_ordered$name, 
                                                                    levels = rev(log_order))
  }
  
  
  bar_errorbar <- ggplot2::ggplot(error_bar_pivot_longer_tibble_summarised_ordered, 
                                  ggplot2::aes(mean, name, fill = group)) + 
    ggplot2::geom_errorbar(ggplot2::aes(xmax = mean + sd, xmin = 0), 
                           position = ggplot2::position_dodge(width = 0.8), 
                           width = 0.5, size = 0.5, color = "black") + 
    ggplot2::geom_bar(stat = "identity", 
                      position = ggplot2::position_dodge(width = 0.8), width = 0.8) + 
    GGally::geom_stripped_cols(width = 10) + ggplot2::scale_fill_manual(values = colors) + 
    ggplot2::scale_color_manual(values = colors) + ggprism::theme_prism() + 
    ggplot2::scale_x_continuous(expand = c(0, 0))+ #, guide = "prism_offset_minor") + 
    ggplot2::scale_y_discrete(labels = rev(daa_results_filtered_sub_df[, x_lab])) + 
    ggplot2::labs(x = "Relative Abundance", y = NULL) + 
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(), 
                   axis.line.y = ggplot2::element_blank(), axis.line.x = ggplot2::element_line(size = 0.5), 
                   axis.ticks.x = ggplot2::element_line(size = 0.5), 
                   panel.grid.major.y = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), 
                   axis.text = ggplot2::element_text(size = 10, color = "black"), 
                   axis.text.x = ggplot2::element_text(margin = ggplot2::margin(r = 0)), 
                   axis.text.y = ggplot2::element_text(size = 10, color = "black", 
                                                       margin = ggplot2::margin(b = 6)), 
                   axis.title.x = ggplot2::element_text(size = 10, color = "black", hjust = 0.5), legend.position = "top", 
                   legend.key.size = ggplot2::unit(0.1, "cm"), legend.direction = "vertical", 
                   legend.justification = "left", legend.text = ggplot2::element_text(size = 8, face = "bold"), 
                   legend.box.just = "right", plot.margin = ggplot2::margin(0, 0.5, 0.5, 0, unit = "cm")) + 
    ggplot2::coord_cartesian(clip = "off")
  if (ko_to_kegg == TRUE) {
    pathway_class_group_mat <- daa_results_filtered_sub_df$pathway_class %>% 
      table() %>% data.frame() %>% column_to_rownames(".")
    pathway_class_group <- data.frame(. = unique(daa_results_filtered_sub_df$pathway_class), 
                                      Freq = pathway_class_group_mat[unique(daa_results_filtered_sub_df$pathway_class), 
                                      ])
    start <- c(1, rev(pathway_class_group$Freq)[1:(length(pathway_class_group$Freq) - 
                                                     1)]) %>% cumsum()
    end <- cumsum(rev(pathway_class_group$Freq))
    ymin <- start - 1/2
    ymax <- end + 1/2
    nPoints <- length(start)
    pCol <- c("#D51F26", "#272E6A", "#208A42", "#89288F", 
              "#F47D2B", "#FEE500", "#8A9FD1", "#C06CAB", "#E6C2DC", 
              "#90D5E4", "#89C75F", "#F37B7D", "#9983BD", "#D24B27", 
              "#3BBCA8", "#6E4B9E", "#0C727C", "#7E1416", "#D8A767", 
              "#3D3D3D")[1:nPoints]
    pFill <- pCol
    for (i in 1:nPoints) {
      bar_errorbar <- bar_errorbar + ggplot2::annotation_custom(
        grob = grid::rectGrob(gp = grid::gpar(col = pCol[i],fill = pFill[i], lty = NULL, lwd = NULL, alpha = 0.2)),
        xmin = ggplot2::unit(-2, "native"), xmax = ggplot2::unit(0, "native"), ymin = ggplot2::unit(ymin[i], "native"), 
        ymax = ggplot2::unit(ymax[i], "native"))
    }
  }
  daa_results_filtered_sub_df <- cbind(daa_results_filtered_sub_df, 
                                       negative_log10_p = -log10(daa_results_filtered_sub_df$p_adjust), 
                                       group_nonsense = "nonsense", log_2_fold_change = NA)
  for (i in daa_results_filtered_sub_df$feature) {
    mean <- error_bar_pivot_longer_tibble_summarised_ordered[error_bar_pivot_longer_tibble_summarised_ordered$name %in% 
                                                               i, ]$mean
    daa_results_filtered_sub_df[daa_results_filtered_sub_df$feature == 
                                  i, ]$log_2_fold_change <- log2(mean[1]/mean[2])
  }
  if(order_by_log == F) {
    daa_results_filtered_sub_df$feature <- factor(daa_results_filtered_sub_df$feature, 
                                                  levels = rev(daa_results_filtered_sub_df$feature))
  } else {
    daa_results_filtered_sub_df$feature <- factor(daa_results_filtered_sub_df$feature, 
                                                  levels = rev(log_order))
  }
  
  p_values_bar <- daa_results_filtered_sub_df %>% 
    ggplot2::ggplot(ggplot2::aes(feature, log_2_fold_change, fill = group_nonsense)) + 
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.8), width = 0.8) + 
    ggplot2::labs(y = "log2 fold change", x = NULL) + GGally::geom_stripped_cols() + 
    ggplot2::scale_fill_manual(values = "#87ceeb") + ggplot2::scale_color_manual(values = "#87ceeb") + 
    ggplot2::geom_hline(ggplot2::aes(yintercept = 0), linetype = "dashed", 
                        color = "black") + 
    ggprism::theme_prism() +
    # ggplot2::scale_y_continuous(expand = c(0, 0), guide = "prism_offset_minor") + 
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(), 
                   axis.line.y = ggplot2::element_blank(), axis.line.x = ggplot2::element_line(size = 0.5), 
                   axis.ticks.x = ggplot2::element_line(size = 0.5), panel.grid.major.y = ggplot2::element_blank(), 
                   panel.grid.major.x = ggplot2::element_blank(),                                                                                                         axis.text = ggplot2::element_text(size = 10,color = "black"), axis.text.y = ggplot2::element_blank(), 
                   axis.text.x = ggplot2::element_text(size = 10, color = "black", margin = ggplot2::margin(b = 6)), 
                   axis.title.x = ggplot2::element_text(size = 11, color = "black", hjust = 0.5), legend.position = "non") + 
    ggplot2::coord_flip()
  if (ko_to_kegg == TRUE) {
    pathway_class_y <- (ymax + ymin)/2 - 0.5
    pathway_class_plot_df <- as.data.frame(cbind(nonsense = "nonsense", 
                                                 pathway_class_y = pathway_class_y, 
                                                 pathway_class = rev(unique(daa_results_filtered_sub_df$pathway_class))))
    if(wrap_label==T){
      pathway_class_plot_df$pathway_class = stringr::str_wrap(pathway_class_plot_df$pathway_class,wraplength)
    }
    pathway_class_plot_df$pathway_class_y <- as.numeric(pathway_class_plot_df$pathway_class_y)
    if(nrow(pathway_class_plot_df)==2){ 
      pathway_class_plot_df$pathway_class_y = 1; pathway_class_plot_df = pathway_class_plot_df[1,,drop=F]
    }
    pathway_class_annotation <- pathway_class_plot_df %>% 
      ggplot2::ggplot(ggplot2::aes(nonsense, pathway_class_y)) + 
      ggplot2::geom_text(ggplot2::aes(nonsense, pathway_class_y, 
                                      label = pathway_class), size = 3.5, color = "black", 
                         fontface = "bold", family = "sans") + ggplot2::scale_y_discrete(position = "right") + 
      ggprism::theme_prism() + ggplot2::theme(axis.ticks = ggplot2::element_blank(), 
                                              axis.line = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank(), 
                                              panel.grid.major.x = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                                              axis.text = ggplot2::element_blank(), plot.margin = ggplot2::unit(c(0, 0.2, 0, 0), "cm"), 
                                              axis.title.y = ggplot2::element_blank(), 
                                              axis.title.x = ggplot2::element_blank(), legend.position = "non")
  }
  daa_results_filtered_sub_df$p_adjust <- as.character(daa_results_filtered_sub_df$p_adjust)
  daa_results_filtered_sub_df$unique <- nrow(daa_results_filtered_sub_df) - 
    seq_len(nrow(daa_results_filtered_sub_df)) + 1
  daa_results_filtered_sub_df$p_adjust <- substr(daa_results_filtered_sub_df$p_adjust, 
                                                 1, 5)
  p_annotation <- daa_results_filtered_sub_df %>% ggplot2::ggplot(ggplot2::aes(group_nonsense,p_adjust)) + 
    ggplot2::geom_text(ggplot2::aes(group_nonsense, unique, label = p_adjust), size = 3.5, color = "black",fontface = "bold", family = "sans") +
    ggplot2::labs(y = "p-value (adjusted)") + 
    ggplot2::scale_y_discrete(position = "right") + ggprism::theme_prism() + 
    ggplot2::theme(axis.ticks = ggplot2::element_blank(), 
                   axis.line = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank(), 
                   panel.grid.major.x = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                   axis.text = ggplot2::element_blank(), 
                   plot.margin = ggplot2::unit(c(0,0.2, 0, 0), "cm"), 
                   axis.title.y = ggplot2::element_text(size = 11, color = "black", vjust = 0), 
                   axis.title.x = ggplot2::element_blank(), 
                   legend.position = "non")
  if (p_value_bar == TRUE) {
    if (ko_to_kegg == TRUE) {
      combination_bar_plot <- pathway_class_annotation + 
        bar_errorbar + p_values_bar + p_annotation + 
        patchwork::plot_layout(ncol = 4, widths = c(1, 
                                                    1.2, 0.5, 0.1))
    } else {
      combination_bar_plot <- bar_errorbar + p_values_bar + 
        p_annotation + patchwork::plot_layout(ncol = 3, 
                                              widths = c(2.3, 0.7, 0.3))
    }
  } else {
    if (ko_to_kegg == TRUE) {
      combination_bar_plot <- pathway_class_annotation + 
        bar_errorbar + p_annotation + patchwork::plot_layout(ncol = 3, 
                                                             widths = c(1, 1.2, 0.1))
    } else {
      combination_bar_plot <- bar_errorbar + p_annotation + 
        patchwork::plot_layout(ncol = 2, widths = c(2.5, 
                                                    0.2))
    }
  }
  return(combination_bar_plot)
}

metacyc_abundance_filtered_location <- metacyc_abundance_filtered_location %>% ungroup()

metacyc_pathway_location_organized <- pathway_errorbar_onesigresult(
  abundance = metacyc_abundance_filtered_location %>% column_to_rownames("pathway"),
  daa_results_df = metacyc_daa_annotated_results_df_location,
  wrap_label = F, wraplength = 20, fc_cutoff = 2,
  Group = metadata_cmb$cohort,
  ko_to_kegg = FALSE,
  p_values_threshold = 0.05,
  order = "group",
  select = NULL,
  p_value_bar = TRUE,
  colors = NULL,
  x_lab = "description"
)
