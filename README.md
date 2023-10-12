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
| | | |
| | | |
| | | |
| | | |
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
        -  Other aims to consider: Functional analysis (using PyCrust), Generate Taxonomic Barplots
    - Weekly timeframe - copy gantt chart in examples; include as much information as possible
    - Overview flowchart - add as many parameters as possible, mentioning specific tools (QIIME2 VS R), diversity metrics, statistical calculation method, and differential abundance measurement used
 - Differentiating and incorporating cohorts:
    - Provide brief overview for each cohort, specifically focusing on geographical differences
    - Can refer to the different living standards and their impact on the microbiome
    - Would need to find and reference papers to explain experimental conclusions
      
### Tasks for the week
- Begin proposal content
    - Split responsibilities internally between general teams
  
