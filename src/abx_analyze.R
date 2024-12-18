library(ggplot2)
library(dplyr)
library(forcats)
library(ggbeeswarm)
library(data.table)
library(pheatmap)
library(RColorBrewer)

infile <- read.csv("/Users/dylanmaghini/scg4/projects/awigen2/binning/05.treponema/04.trep_response/07.abx/out/abx_summary_minimal.tsv", header = TRUE, sep = "\t")
infile <- infile %>% mutate(Site=substr(Sample, 1, 2))

checkm <- read.csv("/Users/dylanmaghini/scg4/projects/awigen2/binning/05.treponema/04.trep_response/02.checkm/checkm_pass.csv", header = TRUE, sep = ",")
names(checkm) <- c("Sample")

sitedf <- data.frame(Site = c("AG", "BF", "DM", "GC", "GH", "GU", "KY", "RE", "SW"), 
                     SiteName = c("Agincourt", "Nanoro", "DIMAMO", "Ref", "Navrongo", "UHGG", 
                                  "Nairobi", "Carter", "Soweto"))

infile <- merge(infile, sitedf, by = c("Site"), all.x = TRUE)
infile <- mutate(infile, Study = ifelse(Site %in% c("AG", "DM", "KY", "SW", "GH", "BF"), "AWIGen", SiteName))
infile_filtered <- infile %>% filter(Sample %in% checkm$Sample)
infile_filtered_quality <- infile_filtered %>% filter(Cut_Off != "Loose")

write.table(infile_filtered_quality, "/Users/dylanmaghini/Downloads/abx_results.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

infile_filtered_quality_genomecount <- infile_filtered_quality %>% 
  group_by(Sample, SiteName, Best_Hit_ARO) %>% 
  summarise(count = n())
infile_filtered_quality_genomecount <- infile_filtered_quality_genomecount %>% 
  mutate(Best_Hit_ARO = gsub(" gene in ", " (", gsub("cluster", "cluster)", Best_Hit_ARO)))
infile_filtered_quality_genomecount <- infile_filtered_quality_genomecount %>% 
  mutate(presence = ifelse(count == 0, 0, 1))
infile_filtered_quality_genomecount <- infile_filtered_quality_genomecount %>% 
  select(!count)

### calculate overall prevalence
pres_abs_long <- infile_filtered_quality_genomecount %>% 
  filter(substr(Sample, 1, 2) %in% c("AG", "BF", "GH", "DM", "KY", "SW"))
pres_abs_abx <- pres_abs_long %>% group_by(Best_Hit_ARO) %>% 
  summarise(count = sum(presence))
total_genome_count = 244 

write.table(pres_abs_long, "/Users/dylanmaghini/Downloads/sourcedata/fig2g.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

pres_abs_abx <- pres_abs_abx %>% mutate(prevalence = count / total_genome_count)

neworder <- c("nimJ","OXA-347",
              "aadS","ErmF",
              "CfxA3", "CfxA5","CfxA6",
              "tet(O)","tet(W)","tet(X)","tet(32)",
              "vanR (vanG cluster)","vanY (vanB cluster)","vanY (vanG cluster)",
              "vanY (vanF cluster)","vanH (vanD cluster)","vanT (vanG cluster)")

pres_abs_abx$Best_Hit_ARO <- fct_relevel(pres_abs_abx$Best_Hit_ARO, "nimJ","OXA-347",
                                         "aadS","ErmF",
                                         "CfxA3", "CfxA5","CfxA6",
                                         "tet(O)","tet(W)","tet(X)","tet(32)",
                                         "vanR (vanG cluster)","vanY (vanB cluster)","vanY (vanG cluster)",
                                         "vanY (vanF cluster)","vanH (vanD cluster)","vanT (vanG cluster)")

ggplot(pres_abs_abx, aes(x = Best_Hit_ARO, y = prevalence)) + 
  geom_col() + 
  theme_bw() + 
  ylab("Prevalence") + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black"), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        panel.background = element_rect(color="black", linewidth = 1), 
        axis.text.y = element_text(color="black"), 
        axis.title.y = element_text(color = "black"))

temp_df <- infile_filtered_quality %>% select(Best_Hit_ARO, Drug.Class, Resistance.Mechanism, AMR.Gene.Family, Antibiotic) %>% 
  unique() 
ggsave("/Users/dylanmaghini/scg4/projects/awigen2/binning/05.treponema/04.trep_response/07.abx/prevalence_bar_ordered.pdf", w = 5, h = 2, dpi = 300)

  
