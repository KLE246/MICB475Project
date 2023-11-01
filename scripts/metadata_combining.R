library(tidyverse)
library(readxl)
anemia_metadata <- read_excel("metadata/anemia_metadata.xlsx")
infant_metadata <- read_excel("metadata/infant_metadata.xlsx")

anemia_cols <- anemia_metadata %>%
  filter(age_months == 6, grepl("BM", diet, fixed=TRUE))%>%
  select(`#SampleID`, "age_months", "sex", "diet") %>%
  mutate(cohort = "anemia") %>%
  mutate(sex = recode(sex, `M` = "male", `F` = "female")) %>%
  rename(age = age_months)

infant_cols <- infant_metadata %>%
  filter(age_category == "6 months") %>%
  filter(feed == "combined" | feed == "breast") %>%
  rename(diet = feed) %>%
  select(`#SampleID`, "age_category", "sex", "diet") %>%
  mutate(cohort = "infant") %>%
  rename(age = age_category) %>%
  mutate(age = word(age , 1  , -2))

combined_df <- rbind(anemia_cols, infant_cols)

write.table(combined_df, file="metadata/combined_md.txt", sep="\t", quote=FALSE)

