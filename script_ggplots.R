### Script to viszualise your data with ggplot###

library(tidyverse)
library(reshape2)
library(dplyr)
library(ggplot2)


#---- create num dataframe ----
# num = number of bacterial hits for the given Category across all samples

# set filter; you can change "Genus" to every other Category

category_num = merged %>% filter(Category == "Genus")

# identify the columns corresponding to the sample names
columns_num = grep("^(aca|ass|hga|hls|hti|nps|ona|osa|ppr1|ppr2|sms)$", 
                    colnames(category_num), value = TRUE)

# select the 'Category', 'Value', and the sample columns
category_num_data = category_num %>%
  select(Category, Value, all_of(columns_num))

# dataframe as txt file (if you want)
write_delim(category_num_data, "/home/tete/Documents/Project/RStudio/Unmapped/Infos/category_num_data_genus.txt", delim = "\t")

# create a long dataframe 
category_num_long <- melt(category_num_data, id.vars = c("Category", "Value"))

#---- Visualize ur filtered dataframe ----

ggplot(category_num_long, aes(x = variable, y = value, fill = Value)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() + 
  scale_x_discrete(labels = sample_name) +
  labs(x = "Samples", y = "Bacteria Match Counts", 
       title = paste("Bacteria Match Counts for each Sample on", category_num_long$Category[1],
                     "Level"),
       fill = "Bacteria") 
ggsave("/home/tete/Documents/Project/RStudio/Unmapped/Figures/ggplot_num_genus.png", 
       width = 14, height = 8, dpi = 300)


#---- create datframe for filtered per ----

# set filter for ("Category") 
category_per = merged %>% filter(Category == "Genus")

# select only the columns that end with '_per' and the Category and Value columns
columns_per = grep("_per$", colnames(category_per), value = TRUE)
category_per_data <- category_per %>%
  select(Category, Value, all_of(columns_per))

# create a long data frame 

category_per_long = melt(category_per_data, id.vars = c("Category", "Value"))

#---- Visualize ur filtered dataframe ----

# ggplot for per value = proportion of each bacteria hit
ggplot(category_per_long, aes(x = variable, y = value, fill = Value)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_bw() + 
  scale_x_discrete(labels = sample_name) +
  scale_y_continuous(breaks = seq( 0, 1, by = 0.2 )) +
  labs(x = "Samples", y = "Portion", 
       fill = "Bacteria",
       title = paste("Portion of Bacteria Matches on", category_per_long$Category[1], 
                     "Level"))
  ggsave("Documents/Project/RStudio/Unmapped/Figures/ggplot_per_genus.png",
        width = 14, height = 8, dpi = 300)


