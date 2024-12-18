library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(stringr)
library(ggpubr)
library(ggupset)

# Read in the overall aggregated dRep information table
aggregated <- read.csv("/Users/dylanmaghini/scg4/projects/awigen2/response/external_compare/compare_all/01.drep/out/data_tables/Cdb.csv", header = TRUE)

#### OVERALL NOVELTY - count by study #####
aggregated <- aggregated %>% mutate(Study = gsub("_.*", "", genome))
aggregated <- aggregated %>% mutate(Study = ifelse(substr(Study, 1, 3) == "REF", "Carter", 
                                                   ifelse(substr(Study, 1, 3) == "MGY", "UHGG", 
                                                          ifelse(substr(Study, 1, 2) %in% c("AG", "BF", "GH", "DM", "KY", "SW"), "AWIGen", Study))))

# pull list of secondary clusters with UHGG 
uhgg <- aggregated %>% filter(Study == "UHGG") %>% select(secondary_cluster) %>% mutate(Type = "InUHGG")
aggregated <- aggregated %>% mutate(Type = ifelse(secondary_cluster %in% uhgg$secondary_cluster, "InUHGG", "Novel"))

# aggregate by study
aggregated_novelty <- aggregated %>% group_by(Study, Type) %>% 
  summarise(Count = n())

aggregated_novelty$Study <- str_to_title(aggregated_novelty$Study)
aggregated_novelty <- aggregated_novelty %>% mutate(Study = ifelse(Study == "Awigen", "AWIGen", Study))
aggregated_novelty$Type <- as.factor(aggregated_novelty$Type)
aggregated_novelty <- aggregated_novelty %>% mutate(Type = fct_relevel(Type, "Novel", "InUHGG"))

aggregated_novelty$Study <- as.factor(aggregated_novelty$Study)
aggregated_novelty <- aggregated_novelty %>% mutate(Study = fct_relevel(Study, "AWIGen", "Uhgg", "Carter", "Yachida", "Franzosa", "Schirmer", "Lochlainn", ))

#### UPSET PLOT ####
upset_prep <- aggregated %>% select(secondary_cluster, Study) %>% mutate(Value = 1)
upset_prep <- unique(upset_prep)
wide_upset <- pivot_wider(upset_prep, names_from = Study, values_from = Value, values_fill = 0)
names(wide_upset) <- c("secondary_cluster", "Franzosa", "Lochlainn", "Schirmer", "Yachida", "AWIGen", "Carter", "UHGG")

wide_upset_matrix <- t(as.matrix(wide_upset))
colnames(wide_upset_matrix) <- wide_upset_matrix[1, ]
wide_upset_matrix <- wide_upset_matrix[-1, ]

tidy_matrix <- wide_upset_matrix %>% 
  as_tibble(rownames = "Study") %>% 
  gather(Cluster, Member, -Study) %>% 
  filter(Member == 1) %>% select(- Member)

tidy_matrix_grouped <- tidy_matrix %>% group_by(Cluster) %>% summarize(Studies = list(Study))
temp_mat <- tidy_matrix_grouped %>% group_by(Studies) %>% summarise(count = n())

filter_list <- c("UHGG", "Carter", "AWIGen", "Franzosa", "Lochlainn", 
                 "Schirmer", "Yachida", 
                 'c("AWIGen", "UHGG")', 'c("AWIGen", "Carter")', 'c("Carter", "UHGG")', 
                 'c("Yachida", "UHGG")', 'c("Franzosa", "UHGG")', 'c("Lochlainn", "UHGG")', 
                 'c("Schirmer", "UHGG")', 
                 'c("AWIGen", "Carter", "UHGG")', 
                 'c("Franzosa", "Lochlainn", "Schirmer", "Yachida", "AWIGen", "Carter", "UHGG")')


upsetplot <- ggplot(tidy_matrix_grouped %>% filter(Studies %in% filter_list), aes(x = Studies)) + 
  geom_bar() + 
  scale_x_upset(order_by = "degree") + 
  theme_bw() + 
  ylab("Genome Clusters") + 
  theme(panel.grid.major.x = element_blank())
upsetplot

ggsave("/Users/dylanmaghini/scg4/projects/awigen2/00.plots/mag_upset.pdf", plot = upsetplot, dpi = 300, w = 6.9, h = 3)

tidy_matrix_grouped$Studies <- as.character(tidy_matrix_grouped$Studies)
write.table(tidy_matrix_grouped, "/Users/dylanmaghini/Downloads/sourcedata/ed7d.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

temp_tab <- tidy_matrix_grouped %>% group_by(Studies) %>% 
  summarise(n = n())
