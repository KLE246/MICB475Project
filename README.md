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
