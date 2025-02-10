### Script to merge all dataframes of all samples / BLASTn results

## run script_unmapped_all.R first

library(tidyverse)
library(reshape2)
library(dplyr)

#---- merge all samples together 

df_list = list(Aca, Ass, Hga, Hls, Hti, Nps, Ona, Osa, Ppr1, Ppr2, Sms)

sample_name = c("Aca", "Ass", "Hga", "Hls", "Hti", "Nps", "Ona", "Osa",
                   "Ppr1", "Ppr2", "Sms")

merged = Reduce(function(x, y) {
  merge(x, y, by = c("Category", "Value"), all = TRUE)
}, df_list)


#---- merges all df_xxx

# list of dataframes for all samples

df_list_df = list(df_aca, df_ass, df_hga, df_hls, df_hti, df_nps, df_ona, df_osa, 
                  df_ppr1, df_ppr2, df_sms)

sample_names_df = c("aca", "ass", "hga", "hls", "hti", "nps",
                     "ona", "osa", "ppr1", "ppr2", "sms")

# add sample-specific columns and count occurrences
for (i in seq_along(df_list_df)) {
  df_list_df[[i]] <- df_list_df[[i]] %>%
    select(Domain, Phylum, Class, Order, Family, Genus, Species) %>%
    group_by(Domain, Phylum, Class, Order, Family, Genus, Species) %>%
    summarise(Sample_Count = n(), .groups = "drop") %>%  # Count occurrences of each unique taxonomy
    mutate(Sample = sample_names_df[i])  # Add sample name
}
# combine all dataframes into one
combined_df = do.call(rbind, df_list_df)

# summarize across all samples
merged_df = combined_df %>%
  pivot_wider(
    names_from = Sample,
    values_from = Sample_Count,
    values_fill = 0  # fill missing counts with 0
  ) %>%
  group_by(Domain, Phylum, Class, Order, Family, Genus, Species) %>%
  summarise(across(everything(), sum), .groups = "drop")  # sum counts for identical entries



