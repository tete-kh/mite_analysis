### Import BLASTn Results to RStudio and Create Dataframe ###

library(tidyverse)
library(reshape2)
library(dplyr)

# the following shows only the code for the Aca BLAStn file
# but this is done for every BLASTn file (I didnt use a loop) 

#---- Aca ----

# load data

blastout = read_delim(
  "/home/tete/Documents/Project/Blast/Unmapped/Formatted/Aca_formatted_unmapped.txt",
  delim = "\t", col_name = F) 

# add column names

colnames(blastout) = c("qseqid", "stitle", "evalue", "bitscore", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send")

# split stitle column

df = separate(blastout, col = stitle, 
              into = c("Accession","Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "X", "Y", "Z"), 
              sep =";")

# merge species colum

df_aca = unite(df, col = "Species", c("Species","X", "Y", "Z"), sep = " ")
df_aca$Species = gsub("NA", "", df_aca$Species) # replace the NA in Species col with void / blank 

# get count data and filter 

Aca_all = df_aca %>%
  select(Domain, Order, Phylum, Class, Order, Family, Genus, Species) %>% 
  pivot_longer(cols = everything(), names_to = "Category", values_to = "Value") %>%
  group_by(Category, Value) %>%
  dplyr::summarise(aca = n(), .groups = "drop") %>%
  arrange(Category, desc(aca))

# filtered data 

df_aca_bacteria = df_aca %>%
  filter(Domain == "Bacteria")  # keep only rows where Domain is Bacteria

Aca = df_aca_bacteria %>%
  select(Domain, Order, Phylum, Class, Family, Genus, Species) %>%  # select relevant columns
  pivot_longer(cols = everything(), names_to = "Category", values_to = "Value") %>%
  filter(!(Category == "Genus" & Value %in% c("marine", "17", "metagenome", "Verruc-01", "Incertae"))) %>% # filter out ambiguous labels  
  group_by(Category, Value) %>%
  dplyr::summarise(aca = n(), .groups = "drop") %>%
  arrange(Category, desc(aca))

# group by the Category column and calculate / add proportion value of each genus (aca_per)  

Aca = Aca %>%
  group_by(Category) %>%
  mutate(aca_per = aca / sum(aca))

