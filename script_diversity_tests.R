### The script for alpha and beta diversity tests ###

library(tidyverse)
library(tibble)
library(reshape2)
library(ggplot2)
library(vegan)
library(dplyr)
library(writexl)

#---- *** alpha diversity *** ----
#---- create long per dataframe + filter ----
# set filter 
category_per = merged %>% filter(Category == "Genus")

# select only the columns that end with '_per' and the Category and Value columns
columns_per = grep("_per$", colnames(category_per), value = TRUE)
category_per_data = category_per %>%
  select(Category, Value, all_of(columns_per))

# create a long data frame 

category_per_long = melt(category_per_data, id.vars = c("Category", "Value"))

#---- calculate relative abundance (p_i) for each taxon in each sample ----

# this is the "value" = xxx_per!!!

#---- calculate Shannon-Index for each sample (variable) ----
shannon_index = category_per_long %>%
  group_by(variable) %>% 
  summarise(Shannon = -sum(value * log(value), na.rm = TRUE))  # this calculates the Shannon Index

# print the Shannon Index for each sample
print(shannon_index)s
columns_per = grep("_per$", colnames(category_per), value = TRUE) 
category_per_data = category_per %>%
  select(Category, Value, all_of(columns_per))
show(shannon_index)
# save shannon_index table
write.csv(shannon_index, "/home/tete/Documents/Project/RStudio/Unmapped/Figures/shannon_index_table_genus.csv", row.names = FALSE)

#---- *** beta diversity *** ----
#---- create long dataframe for num + filter ----

# set filter 
category_num = merged %>% filter(Category == "Genus")

# identify the columns corresponding to the sample names 
columns_num = grep("^(aca|ass|hga|hls|hti|nps|ona|osa|ppr1|ppr2|sms)$", 
                    colnames(category_num), value = TRUE)

# select the 'Category', 'Value', and the sample columns
category_num_data = category_num %>%
  select(Category, Value, all_of(columns_num))


# create a long dataframe 
category_num_long = melt(category_num_data, id.vars = c("Category", "Value"))

community_matrix_transposed = category_num_long %>%
  filter(!is.na(value)) %>%                             
  pivot_wider(names_from = variable,                    
              values_from = value, 
              values_fill = 0) %>%                      
  column_to_rownames(var = "Value") %>%                 
  select(-Category) %>%                                 
  t() %>%  # transpose the matrix so that samples become rows
  as.data.frame() # this converts matrix to a regular dataframe

#---- Filter: exclude samples min hit count ----

# samples to exclude
excluded_samples = c("hti", "nps", "ppr1", "ppr2")

# remove excluded samples from the transposed community matrix
community_matrix_filtered = community_matrix_transposed[!(rownames(community_matrix_transposed) %in% excluded_samples), ]


#---- Bray-Curtis on filtered Samples ----
# calculate Bray-Curtis dissimilarity on the filtered matrix

bray_curtis_filtered = vegdist(community_matrix_filtered, method = "bray")

# Convert the resulting dissimilarity into a matrix
bray_curtis_matrix_filtered = as.matrix(bray_curtis_filtered)

# Update sample names for the filtered matrix
filtered_sample_names = rownames(community_matrix_filtered)
rownames(bray_curtis_matrix_filtered) = filtered_sample_names
colnames(bray_curtis_matrix_filtered) = filtered_sample_names

### Visualiztion ###
# perform hierarchical clustering
hc_filtered = hclust(as.dist(bray_curtis_matrix_filtered), method = "average")

# plot dendrogram for filtered samples
png("/home/tete/Documents/Project/RStudio/Unmapped/Figures/bray_curtis_dendrogram_filtered.png", 
    width = 800, height = 600)

# plot the dendrogram
plot(hc_filtered, 
     main = "Dendrogram of Bray-Curtis Dissimilarity (Filtered)",
     xlab = "Samples", 
     ylab = "Height (Dissimilarity)")

dev.off()


# Generate a heatmap
png("/home/tete/Documents/Project/RStudio/Unmapped/Figures/bray_curtis_heatmap_genus_filtered.png", width = 800, height = 800)
heatmap.2(
  bray_curtis_matrix_filtered,
  trace = "none",
  dendrogram = "none",  # Remove hierarchical clustering
  Rowv = NA,
  Colv = NA,
  col = topo.colors(100),
  main = "Bray-Curtis Dissimilarity Heatmap (Filtered)",
  xlab = "Samples",
  ylab = "Samples",
  key = TRUE,
  key.title = "Dissimilarity",
  key.xlab = "Dissimilarity",
  key.par = list(mar = c(4, 1, 4, 1)),
  keysize = 1.5
)
dev.off()

#---- PERMANOVA - on Filtered ----
# samples to exclude
excluded_samples = c("hti", "nps", "ppr1", "ppr2")

# filter reproductive_mode to exclude the specified samples
reproductive_mode_filtered = reproductive_mode[!(reproductive_mode$sample %in% excluded_samples), ]

# calculate Bray-Curtis distance on the filtered community matrix
bray_curtis_filtered = vegdist(community_matrix_filtered, method = "bray")

# perform PERMANOVA using the filtered data
adonis_result_bray_filtered = adonis2(
  as.dist(bray_curtis_filtered) ~ mode,
  data = reproductive_mode_filtered,
  permutations = 999, 
  method = "bray"
)
print(adonis_result_bray_filtered)

# for the record: the permutation is set to 999, but the function detected a low data size and automatically changed it to 5039 

#---- NMDS for FILTERED ----

# samples to exclude
excluded_samples = c("hti", "nps", "ppr1", "ppr2")

# filter reproductive_mode to exclude the specified samples
reproductive_mode_filtered = reproductive_mode[!(reproductive_mode$sample %in% excluded_samples), ]

# filter the community matrix to match the filtered reproductive_mode
community_matrix_filtered = community_matrix_transposed[!(rownames(community_matrix_transposed) %in% excluded_samples), ]

# calculate Bray-Curtis distance on the filtered community matrix
bray_curtis_filtered = vegdist(community_matrix_filtered, method = "bray")

# Perform NMDS on the filtered Bray-Curtis dissimilarity matrix
nmds_bray_filtered = metaMDS(as.dist(bray_curtis_filtered), k = 2, trymax = 100)

# extract NMDS coordinates
nmds_coordinates_bray_filtered = as.data.frame(nmds_bray_filtered$points)

# add sample names and reproductive mode to the NMDS coordinates
nmds_coordinates_bray_filtered$sample = rownames(nmds_coordinates_bray_filtered)
nmds_coordinates_bray_filtered = merge(nmds_coordinates_bray_filtered, reproductive_mode_filtered, by = "sample")

#---- Visualize Filtered NMDS Bray Results ----

# plot NMDS with sample names above the points
ggplot(nmds_coordinates_bray_filtered, aes(x = MDS1, y = MDS2, color = mode)) +
  geom_point(size = 4) +
  geom_text(aes(label = sample), size = 3, vjust = -1, hjust = 0.5, color = "black") +  
  labs(title = "Filtered NMDS Plot of Bray-Curtis Dissimilarity",
       x = "NMDS Dimension 1",
       y = "NMDS Dimension 2",
       color = "Reproductive Mode") +
  scale_color_manual(values = c("asex" = "blue", "sex" = "red")) +  
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("/home/tete/Documents/Project/RStudio/Unmapped/Figures/nmds_bray_curtis_genus_filtered.png", width = 8, height = 5)


