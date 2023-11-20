# MICB475 Project Lab Notebook

## Team 10  
Risa Fox  
Stella Lin  
Kevin Le  
Michael Qiu  
Dennis Xie

# Table of Contents
|Experiment # | Title | Dates | Link |
| --- | ------------------- | ----- | ---------|
| W5-TM | Week 5 Team Meeting | 2023-10-04 - 2023-10-10 | [jump](#w5-tm---week-5-team-meeting) |
| W6-TM | Week 6 Team Meeting | 2023-10-11 - 2023-10-17 | [jump](#w6-tm---week-6-team-meeting) |
| W7-TM | Week 7 Team Meeting | 2023-10-18 - 2023-10-24 | [jump](#w7-tm---week-7-team-meeting) |
| EXP-1| Summary of Dataset Processing | 2023-10-21 | [jump](#exp-1-infant-dataset-processing-methods) |
| EXP-2| Combining of Both Metadata Sets | 2023-10-27 | [jump](#exp-2-combining-metadata) |
| W8-TM | Week 8 Team Meeting | 2023-10-25 - 2023-10-31 | [jump](#w8-tm---week-8-team-meeting) |
| W9-TM | Week 9 Team Meeting | 2023-11-01 - 2023-11-07 | [jump](#w9-tm---week-9-team-meeting) |
| W10-TM | Week 10 Team Meeting | 2023-11-08 - 2023-11-15 | [jump](#w10-tm---week-10-team-meeting) |
| EXP-3| Deseq Analysis | 2023-11-07 | [jump](#exp-3-deseq-analysis) |
| EXP-4| Deseq Analysis (cohort investigation) | 2023-11-10 | [jump](#exp-4-deseq-analysis-cohort-investigation) |
| W11-TM | Week 11 Team Meeting | 2023-11-16 - 2023-11-22 | [jump](#w11-tm---week-11-team-meeting) |
| | | |


## W5-TM - Week 5 Team Meeting

### Agenda
1. Discuss potential for combining the anemia and feeding datasets. Found that things like sex, weight, height, and age would be shared between the datasets. [EDA findings](https://github.com/KLE246/MICB475Project/blob/main/references/initial_metadata_analysis.xlsx)
2. Discuss some potential questions
3. Creating novel question given our limited overlap between the two sets. Deciding on using one dataset or both

Potential questions:
- Is there a gender-related disparity in the prevalence of anemia in infants, and how might this be related to other factors like dietary habits and infections?
- Are there associations between an infant's dietary habits (appetite, food responsiveness) and anemia status (e.g. growth delays or underweight status) (from the anemia dataset)?
- How do various health parameters (e.g., ferritin levels, hemoglobin) change over time in infants, and can these changes be linked to specific dietary habits or medical interventions?
- At how many months should iron supplements be introduced to an infant’s diet, and how does this alter their dietary habits?
- Are the microbiome changes in the intestines from anemia and diet westernisation comparable? Do these changes affect rates of wasting or stunting?
- Is there a significant weight difference between the two cohorts? If so what is causing this difference? (eg. Iron supplements, feeding type, microbiome)

### Meeting notes:
- Is there a gender-related disparity in the prevalence of anemia in infants, and how might this be related to other factors like dietary habits and infections?
  - This seems novel (according to TA)
  - Combine: controls of anemia dataset, compare dietary habits – assume that kids don’t have anemia?
    - Sex differences between the healthy kids in the anemia dataset (cohort 1), infant dataset (cohort 2)
    - Filter for 6 months, breastmilk feed (?)
    - Look for sex differences with respect to microbiome
    - Compare the two groups: which correlations impacted by what dataset (e.g. which microbes affected by sex)
      - Include dataset as variable
  - For both, approach would be similar

### Tasks for the week
- Finalize dataset and research question
- Start proposal outline
    - Split responsibilities
- \(stretch goal\) Load all data and begin QIIME work


## W6-TM - Week 6 Team Meeting

### Agenda
1. Go through assigned responsibilities for the proposal. Original proposal:
   - Risa - Introduction and Background
   - Kevin, Stella - Research Objective, Experimental Aim, Proposed Approach, Overview Flowchart, Weekly Timeframe
   - Dennis, Michael - Dataset Overview
2. Discussion Areas:
   - Discuss experimental aims - amount of aims we should have and how we should approach this section
   - Starting points to be included in the weekly timeframe
   - Confirm collection sites for the anemia and infant datasets and potential methods of comparing them for the research question

### Meeting notes:
- Reassess assigned responsibilities to split roles evenly:
    - Risa -  Introduction and Background, Weekly Timeframe
    - Kevin, Stella - Research Objective, Experimental Aim, Overview Flowchart
    - Dennis, Michael - Proposed Approach, Dataset Overview
- Breaking down outline:
    - Research objective and experimental aims should be done together
    - Background should consist of a broad introduction giving context; no direct reference to dataset
    - Research objective should consist of a topical overview
    - Experimental aims should consist of specific variables of focus and how we plan to tackle them, providing context on why aims were chosen
        -  Want minimum 3 aims: Diversity, Bacterial taxa, comparison of cohorts
        -  Diversity recommendation: Wilcoxon test (Alpha) and Permanova (Beta)
        -  Differential abundance recommendation: DEseq2
        -  Other aims to consider: Functional analysis (using PICRUSt2), Generate Taxonomic Barplots
    - Weekly timeframe - copy gantt chart in examples; include as much information as possible
    - Overview flowchart - add as many parameters as possible, mentioning specific tools (QIIME2 VS R), diversity metrics, statistical calculation method, and differential abundance measurement used
 - Differentiating and incorporating cohorts:
    - Provide brief overview for each cohort, specifically focusing on geographical differences
    - Can refer to the different living standards and their impact on the microbiome
    - Would need to find and reference papers to explain experimental conclusions
      
### Tasks for the week
- Begin proposal content
    - Split responsibilities internally between general teams


## W7-TM - Week 7 Team Meeting

### Agenda
1. Go through current proposal draft. Current progress:

  | Section | Member(s) | Status |
  | :--- | :---: | :---: |
  | Introduction and Background | Risa | Drafted |
  | Research Objective | Kevin | Drafted |
  | Experimental Aims | Stella | In Progress |
  | Overview Flowchart | Kevin | Not Started |
  | Data Overview | Michael, Dennis | Drafted |
  | Proposed Approach | Michael, Dennis | Not Started |
  | Weekly Timeframe | Risa | Not Started |
 
  WIP for research proposal: https://docs.google.com/document/d/1UzibL7BBn7bLV8ccPwqZRif5IkvXrVEKDcAQhmFeIg0/edit?usp=sharing

2. Discuss how hypothesis should be stated in the Research Objectives
3. Is there proper split in the sections regarding information? Is there overlapping information that has to be removed, specifically between the Introduction and Research Objectives?
4. How do we find what primers were used in the experiments?
5. Discuss whether we should include an Experimental Aim to investigate Functional Analysis?

### Meeting notes:
- Avril highlighted that the goal with the proposals is to clearly outline gaps in the field, and why main outcomes are important and relevant
- Avril gave several notes on the proposal draft
  - Introduction and Background done well; Avril hinted to go through assignment rubric to ensure that every point is fully addressed
    - Missing: Reference to previous UJEMI papers
      - Evelyn would like to see evidence that we've looked into UJEMI papers
      - Can give a brief overview of what was researched previously
    - Risa verified that the Infant Feeding dataset does not have an associated paper
    - Can also mention ideas that are not directly related to our research question
    - Last paragraph: Avril advised to include additional benefits that the outcomes of the research can provide (ie. future research areas, interventions)
  - Research objectives: Avril highlighted specific parts of the draft
    - Paragraph 3: try to give more specific examples of effects of different geological location and connect to the real world
    - Paragraph 2 Sentence 4: try to specify that the factors are from the metadata
      - External factors can be mentioned in the Discussion
  - Experimental aims: need to be more specific and clarify each statistical test outcome
    - Avril recommended for us to not focus on geological locations for aims 1 and 2, but keep for aim 3
    - For alpha and beta diversity, specify outcomes from specific analyses and how they differ from other tests
    - Can look at specific studies and include expected results from previous papers
    - Use words "Bacterial functional phenotype" for aim 3 to clarify obtained functions from measuring DNA
  - Proposed Approach: looks good so far
    - Avril outlined that we should mention the type of Beta Diversity analysis we will be using
    - May want to control which cohort the infants are from by inputting cohort into the y=mx+b formula (y=mx + cohort)
      - Will allow for further correlation between outcomes 
- When stating the hypothesis, try to create a logical flow from reference papers to justify our predicted outcomes
- There is a proper split in the sections currently without too much overlapping information; overlapping is fine as long as we elaborate in further sections
  - The Objective will include the general research, while the Aims will have specific research examples and outcomes from statistical tests
- Infant feeding paper used 515-806 primers (from previous UJEMI paper)
- Anemia paper used 515-806 primers as well
- Avril thinks that including a functional analysis will be doable and an important addition to the other statistical tests
- Other Discussion
  - Will not need to justify type of statistical test chosen
  - Will need to justify reasoning for truncating reads
  - There is no major significance in infants who are exclusively fed breast milk versus infants fed with combination breast milk and other
    - Key: breast milk will offer most difference in gut microbiota
  - Expected workload will be primarily for troubleshooting the codes
  - Plan to create ~3 figures and additional supplementary figures
 

 ### Tasks for the week
- Complete proposal with given edits
- Run proposal through spell and grammar check


## EXP-1 Infant dataset processing methods
1.  Connect to UBC’s VPN and log into MICB 475 class server through secure shell ssh on Windows Terminal.
2.  Create a directory in the data directory called infant: /data/infant/
3.  Import and demultiplex infant data using a manifest. Outputted files: infant_seqs.qza
    - Infant dataset located in /mnt/datasets/project_2/infant/
4.  A visualization of demultiplexed samples was created and scp command was used to transfer the file to local directory. Outputted files: infant_seqs.qzv
5.  Appropriate truncation length was determined through visualization on Qiime 2 View by looking at the interactive quality plot.
6.  Denoise with truncation length 150, left 0. Outputted files: rep-seqsi.qza, tablei.qza, statsi.qza
7.  Visualization for DADA 2 stats were generated and scp command was used to transfer files to the local directory. Outputted files: statsi.qzv
8.  Table (tablei.qza) file was filtered using Qiime for 6 month old infants and infants that were at least partially breastfed. Outputted files: final-filtered-table.qza, final-filtered-table.qzv
    - 6 month old infants were filtered using “[age_category]='6 months '”
    - Partially breastfed infants were filtered using “[feed] IN ('breast', 'combined')”
9.  Taxonomic analysis done by training Silva 138-99 database using 515F-806R primers. Outputted files: ref-seqs-trimmed.qza, classifier.qza, taxonomyi.qza
10.  Tree for phylogenetic diversity analysis was generated using the rep-seqsi.qza file. Outputted files: rooted-treei.qza
11.  Alpha rarefaction generates a rarefaction curve which was used to determine sampling depth. 23,647 was set as maximum sampling depth to retain all samples. Scp command was used to transfer visualization to local directory. Outputted files: alpha-rarefactioni.qzv
12.  Diversity metrics for alpha and beta diversity are done at a sampling depth of 20,000. Outputted directory: core-metrics-results

## EXP-1 Anemia dataset processing methods 
qiime2-2023.7 was used for data processing
1.  Connect to UBC’s VPN and log into MICB 475 class server through secure shell ssh on Windows Terminal.
2.  Create a directory in the data directory called infant: /data/anemia/
3.  Import and demultiplex infant data using a manifest. Outputted files: demux_seqs_anemia.qza
    - Infant dataset located in /mnt/datasets/project_2/anemia/
4.  A visualization of demultiplexed samples was created and scp command was used to transfer the file to local directory. Outputted files: demux_seqs_anemia.qzv
5.  Appropriate truncation length was determined through visualization on Qiime 2 View by looking at the interactive quality plot.
6.  Denoise using DADA 2 package with truncation length 150, left 0. Outputted files: rep-seqs_anemia.qza, table_anemia.qza, stats_anemia.qza
7.  Visualization for DADA 2 stats were generated and scp command was used to transfer files to the local directory. Outputted files: stats_anemia.qzv
8.  Table (tablei.qza) file was filtered using the metadata for 6 month old infants and infants that were at least partially breastfed, as well as not anemic. Outputted files: filtered_table_infant.qza
9.  Alpha rarefaction generates a rarefaction curve which was used to determine sampling depth. 23,647 was set as maximum sampling depth to retain all samples. Scp command was used to transfer visualization to local directory. Outputted files: alpha-rarefaction_anemia.qzv
    
For documentation of code please see [anemia_oct21.23](https://github.com/KLE246/MICB475Project/blob/bb10776ef5fbdb8577f50310681158929169e47d/documentation/anemia_oct21.23.sh)

Oct21.23 DX


## EXP-2 Combining Metadata
The metadata for both cohorts needs to be combined on similar columns that are relevant for the future analysis.  
Columns were determined through initial [EDA findings](https://github.com/KLE246/MICB475Project/blob/main/references/initial_metadata_analysis.xlsx)  
Script used for combining found [HERE](https://github.com/KLE246/MICB475Project/blob/main/scripts/metadata_combining.R)  
Output found [HERE](https://github.com/KLE246/MICB475Project/blob/main/metadata/combined_md.txt)


## W8-TM - Week 8 Team Meeting

### Agenda
1. Discuss new experimental aim added in the proposal: "Determine if geographical location is a factor that influences the gut microbiota diversity in infants of different sexes."
2. Discuss how to split up the responsibilities for the project

### Meeting notes:
- Avril approved of the experimental aim:
  - Analyses will be similar to other experimental aims
  - Focus will be on comparing how microbiota diversity between regions rather than replicating results
- Specifics regarding project + analyses:
  - QIIME2 processing has been completed and needs to be exported
    - For taxonomy,  ensure to justify how classifer was trained
  - Avril mentioned that Alpha diversity will be the easiest to complete, but Beta diversity and differential abundance will be more complicated
  - Avril highlighted that we should create a combined metadata and begin formatting codes for each analysis
    - Metadata should be formated as a R file and put into github for everyone to use
    - Can begin formatting codes for DESeq2, Alpha-, Beta-Diversity
  - Details about each analysis:
    - Alpha- and beta-diversity will use vegan R package; visualized using ggplot; let Avril know of issues when running PERMANOVA
    - DESeq2 code similar to one in module and will generate calcualted statistic
      - May need 2 people to troubleshoot generated Volcano plot to ensure plot is clear; made through ggplot
    - PICRUSt package may have to be installed onto QIIME, otherwise code consists of 1 big code
      - Uses feature table and reads to generate a taxonomy table
  - Truncation length and other parameters needs to be the same length for both datasets
    - Option: combine all reads together and process everything
  - Avril highlighted that we should prioritize compiling metadata, and generating diversity and taxonomic composition reads; phyloseq object will be used in further analyses
    - We should aim to complete a rough analysis within the next 3 weeks to begin working on the rough manuscript
- Role Splitting: *subject to change 

   | Role | # People | Member |
  | :--- | :---: | :---: |
  | Compiling metadata | 1 | Kevin |
  | Alpha Diversity | 1 | Michael |
  | Beta Diversity | 1 | Dennis |
  | Differntial Analysis (DESEq2) | 2 | Risa, Kevin |
  | PICRUSt (QIIME2) | 1-2 | Dennis, Stella |


 ### Tasks for the week
- Dennis and Michael collaborate to align parameters used
- Complete metadata compilation
- Write skeleton codes for each analysis
- Create phyloseq object


## W9-TM - Week 9 Team Meeting

### Agenda
1. Discuss feedback from proposal
2. Update on analysis progress

### Meeting notes:
- Proposal feedback:
  - Experimental aims and research objectives felt unbalanced; more general information in aims can be redirected into research objectives
    - Research objectives goal: provide more evidence to why a specific aspect (ie. geographical regions) is significant and how it would impact the research
  - For resubmission, focus on comments added to rubric
  - When performing correlation, Avril proposed us to either manually compare results of analyses or conduct an interaction effect when conducting analyses
    - Manually: directly compare cohorts and establish a cutoff to determine if results are significantly different
    - Interaction Effect: will test if differences between 2 cohorts is significant
      - Code: Microbe~sex*cohort
      - When conducting: could create box plots to correlate results from male to female in each cohort, then compare slopes between 2 cohorts
  - Avril proposed that we should remove Experimental Aim 5 and merge analyses with aims 3 and 4 
    - Resulting aims: Diversity, Differential Abundance, Bacterial Functional Phenotype, Geographical Comparison
- Update on analysis progress
  - Metadata has been compiled as a .txt file in the repo
    - Phyloseq object will be saved as an R data
  - Roles have been assigned


 ### Tasks for the week
- Complete rough analysis for each aim (skeleton code)

## W10-TM - Week 10 Team Meeting

### Agenda
1. Discuss code progress so far
2. Debug picrust2 error

### Meeting notes:
- Picrust debugging:
  - Need to install the tool in a directory that you have write permissions to
  - Use the same version of picrust as Avril, use the command picrust2_pipeline.py 
- DESEq2: looking good!
  - Lots of upregulation as opposed to downregulation in volcano plot
  - Female was used as control; more upregulation in males
  - Label male/female for log(2) fold change
  - Can also filter out NA 
- Alpha + Beta:
  - Issue with tree
    - Post-combining – Error: duplicate taxa
    - Send error to Avril 
    - Tree fails for both datasets
  - Need to attach tree A to phyloseq A, then tree I to phyloseq I
    - Make each phyloseq object separately then combine them
    - Will get duplicated taxa with this: need to combine trees first then create phyloseq object
  - Proposed: merge phylo first 
    - Send Avril the zipped phyloseq object creation R file
    - Leave tree out of phyloseq for now using types that don’t use tree data
- Will be using ggpicrust2 in R for PICRUSt
  - Will allow us to run DESEq2 – what differential abundance to run, etc
  - Stick to DESEq2 for this
- Results discussion:
  - Figure out significant microbes
  - All numbers at end of each name refers to ASV 
  - Start with microbes with greater significance (option 1)
    - Can focus on microbes with greater fold changes
  - Combine ASVs and rerun (option 2)
  - Start with general overall discussions, then cohort specific findings
  - Can also discuss results different from literature
  - Think about how we would like to portray that
    - Ex. boxplot side-by-side comparison
    - Cohort vs cohort comparison 
    - Internal sex difference of cohort 1 vs internal sex difference of cohort 2
  - Think about where we want each plot + which ones to put into supplementary
    - Ex. volcano plot – likely will go into supplementary
    - Plots not discussed in depth will go into supplementary


### Tasks for the week
- Complete all coding analysis for aim 1-3, aim 4 if possible
- Start writing methods + introductions for the final manuscript
- Begin making presentation
- For deseq analysis, implement Avril's suggestions: make figures presentation ready (remove NAs, male vs female in y-axis, add colours)
- Confirm with Avril whether we’d like a meeting next week + preferred time

## EXP-3 Deseq Analysis
Deseq analysis has been conducted to evaluate how gut microbial composition differs in infants of different sexes. A volcano plot has been created to evaluate differences in microbial taxa abundance between male and female infants. A bar plot has been created to show relative abundance of specific microbial taza and provide a clear representation of microbial prevalence in each sex.

The following pull request shows the R script and resulting plots: https://github.com/KLE246/MICB475Project/pull/30/files.


## EXP-4 Deseq Analysis (cohort investigation)
As stated in Aim 4, analysis of previous aims are to be done separately between cohorts.

[Original Deseq bar plot with separated cohorts](https://github.com/KLE246/MICB475Project/blob/main/plots/separated_cohort_bar_plot.png)

[Looking at only similar Genus](https://github.com/KLE246/MICB475Project/blob/main/plots/similar_cohort_bar_plot.png)

Results indicate only one of the genus in the Deseq analysis were similar between the cohorts.  
Analyis of these factors could be used for discussion  
- role in gut, oral cavity  
- role in diet  
- role in environment  

Since the level of sex difference varies between the two locations (comparatively higher in male in infant cohort), conclusions could be made that one location has some factors that push for a more unequal presence. Further discussion required.

Update:

Additional changes were made to look at only the similar Genus after grouping to allow observation of trends, not just relying on statistical significant results. 

[Plot with similar Genus, highest changes](https://github.com/KLE246/MICB475Project/blob/main/plots/separated_similar_genus_bar_plot.png)

## W11-TM - Week 11 Team Meeting

- Bottom of lab notebook – 2 links to 2 plots
  - Change in genus abundance between male + female plot
    - Only 1 genus in common between both cohorts – based off DESeq2 analysis plotting significant ones
    - Could be due to filtering leading to a smaller dataset; anemia may be more powerful as there’s more datasets
    - plot values for ALL of microbes so every column has a value for anemia + infant; make it a bit transparent if not significant
      - Could just not be significant due to difference in power
    - Genus shared between cohorts would be good; if there are differences could be due to geographical factors
    - Potential workaround: May see more overlap at genus level rather than asv level
      - ASV’s more likely to be more specific to the population studying
- Bar plot 2 – male/female change
    - Looking at the one genus similar to both
    - Would be the plot we want for aim 4 – looking at just 1 microbe isn’t as robust as we’d like
      - Can just look at genus level to see for we can get more results
        - Avril suggests making box plots if we plan on subsetting a specific microbe
          - Can convert to relative abundance, then show amount for male/females in box plot and put them side by side
          - Result – 4 boxes in 1 boxplot 
            - To see what the actual abundances are and give more information 
        - Would be typical to have DESeq as main results; box plots and subsetted plots can be supplementary
          - Could make a page of box plots and put in supplementary and reference when needed
        - Avril recommends adding in values for all microbes even if not significant
          - Need a lot more power for infant one – could explain insignificance
          - Aim: have plot like deseq2 one, show significance/non-signifiance using  transparency, then discuss results – can refer to results visually rather than statistically
  - Command for aggregating genus: tax_glom('Species')
    - Will be done directly to phyloseq object before subsetting 
- Shannon Alpha Diversity
  - Similar story to other boxplot
  - Basically none of alpha diversity had any significance
    - Comment: pretty typical for alpha diversity to have much significance
      - Should be fine for it to not be significant – can just say “Alpha diversity is not significant”
  - Asked if it’s ok if we just briefly talk about it 
    - Worth it to talk about DESeq2 but with something like Alpha diversity, it’s fine to just leave it this way
- ggpicrust2
  - Calculate statistics for pca for sex – worth running
  - Show – overall pathway + function = subtle
  - PCA = good
  - Metacyc – include pathways; could say that these are significant 
    - Take one and make a box plot of pathway with each group of interest
  - Kegg – if there aren’t significant results overall, don’t worry too much
    - Focus on the ones that ARE significant
    - Since infant isn’t as big, it might result in difference in power
  - Plotting – plot all 4 groups and see significance
    - Then see trends
  - Don’t HAVE to analysis everything – can just focus on one analysis pathway
- Beta diversity
  - No significant when looking at sex, but there is significance when looking at cohort
  - Can we put more focus on cohort?
    - May lose some of the work already done with sex
  - Avril – potentially could focus more on geography if sex is not working out well (would need to run with Evelyn)
  - Could also flip it – focus on geography and then use sex differences as a supplementary
- For next week – need to start building up manuscript
  - Make a big google slides with all figures to start piecing together what the overall story is
  - Don’t need to include ALL data but need to include all necessary data to tell the story
- GOALS:
  - Aggregate everything into document
  - Should have an understanding of what the paper should look like, what figures will be in the paper
    - Manuscript skeleton/outline OR manuscript draft (stretch)
- Make sure to send draft manuscript to Avril ASAP to give him enough time to work through it

