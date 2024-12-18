#library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(here)
library(cowplot)
library(scales)
library(ggpubr)
library(forcats)


#### GENOME NOVELTY ####
cdb <- read.csv("/Users/dylanmaghini/scg4/projects/awigen2/binning/06.uhgg_compare/out/data_tables/Cdb.csv", header=TRUE, sep=",")
cdb <- cdb %>% mutate(Study = ifelse(substr(genome, 1, 3) == "MGY", "UHGG", "AWIGen2"))

cluster_ids <- cdb %>%
  group_by(secondary_cluster, Study) %>%
  summarise(entries_count = n()) %>%
  pivot_wider(names_from = Study, values_from = entries_count, values_fill = 0)

cluster_ids <- cluster_ids %>% mutate(ClusterType = ifelse(AWIGen2 == 0 & UHGG > 0, "UHGG Only",
                                                           ifelse(AWIGen2 > 0 & UHGG == 0, "AWI-Gen 2", "Shared")))


sitedf <- data.frame(SiteCode = c("AG", "BF", "GH", "KY", "SW", "DM"),
                     Site = c("Bushbuckridge", "Nanoro", "Navrongo", "Nairobi", "Soweto", "DIMAMO"))

cdb_awi <- cdb %>% filter(Study == "AWIGen2")
cdb_awi <- merge(cdb_awi, cluster_ids, by="secondary_cluster", all.x = TRUE)
cdb_awi <- cdb_awi %>% mutate(SiteCode = substr(genome, 1, 2))
cdb_awi <- merge(cdb_awi, sitedf, by="SiteCode")
cdb_awi <- cdb_awi %>% mutate(ClusterType = gsub("Shared", "In UHGG", gsub("AWI-Gen 2", "Novel", ClusterType)))
cdb_awi$ClusterType <- as.factor(cdb_awi$ClusterType)
cdb_awi$ClusterType <- relevel(cdb_awi$ClusterType, 'Novel')

awi_new <- cdb_awi %>% filter(ClusterType == "Novel") %>% nrow()
awi_old <- cdb_awi %>% filter(ClusterType == "In UHGG") %>% nrow()
bacteria_compare <- ggplot(cdb_awi, aes(x=Site, fill=ClusterType)) + 
  geom_bar() + 
  ylab("Bacterial Genomes") +
  #scale_fill_manual(values = c("#a3c3c9", "#45818e"))+ 
  scale_fill_manual(values = c("#cfcfcf", "#969696"))+ 
  scale_x_discrete(limits = c("Nanoro", "Navrongo", "DIMAMO", "Bushbuckridge", "Soweto", "Nairobi")) + 
  theme_bw()

bacteria_compare


write.table(cdb_awi %>% mutate(site = ifelse(SiteCode == "AG", "Agincourt", 
                                             ifelse(SiteCode == "DM", "DIMAMO", 
                                                    ifelse(SiteCode == "GH", "Navrongo", 
                                                           ifelse(SiteCode == "BF", "Nanoro", 
                                                                  ifelse(SiteCode == "KY", "Nairobi", 
                                                                         "Soweto")))))), "/Users/dylanmaghini/Downloads/sourcedata/fig2c.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#### PROTEIN NOVELTY ####
proteins <- read.csv("/Users/dylanmaghini/scg4/projects/awigen2/proteins/231202_uhgp_compare/uhgp_compare_counts.tsv", sep="\t", header=TRUE)

names(proteins) <- c("SiteCode", "In UHGP", "Novel")
proteins <- merge(proteins, sitedf, by="SiteCode")

proteins_long <- melt(proteins, id.vars=c("SiteCode", "Site"))
names(proteins_long) <- c("SiteCode", "Site", "ClusterType", "Count")

proteins_long$ClusterType <- as.factor(proteins_long$ClusterType)
proteins_long$ClusterType <- relevel(proteins_long$ClusterType, 'Novel')

# this is where to get the total value for the manuscript
proteins_long %>% group_by(ClusterType) %>% summarize(total = sum(Count))

write.table(proteins_long, "/Users/dylanmaghini/Downloads/sourcedata/ed7a.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)
proteins_compare <- ggplot(proteins_long, aes(x=Site, y=Count, fill=ClusterType)) + 
  geom_col() + 
  ylab("Bacterial Proteins") +
  #scale_fill_manual(values = c("#a3c3c9", "#45818e"))+ 
  scale_fill_manual(values = c("#cfcfcf", "#969696"))+ 
  scale_x_discrete(limits = c("Nanoro", "Navrongo", "DIMAMO", "Bushbuckridge", "Soweto", "Nairobi"), 
                   labels = c("Nanoro", "Navrongo", "DIMAMO", "Agincourt", "Soweto", "Nairobi")) + 
  theme_bw() 
proteins_compare

#### VIRAL NOVELTY ####

site_df <- data.frame(Site=c("AG", "BF", "DM", "GH", "KY", "SW"), 
                      SiteName = c("Bushbuckridge", "Nanoro", "DIMAMO", "Navrongo", "Nairobi", "Soweto"))

class_df <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/phages/class_counts.tsv", header=TRUE, sep="\t")
class_df_long <- melt(class_df, id.vars=c("Site"))
names(class_df_long) <- c("Site", "Class", "Count")
class_df_long <- class_df_long %>% mutate(Class = gsub("NULL.", "NoClassification", Class)) %>% 
  mutate(Class = gsub("NoMGVRep", "NoMGVMatch", Class)) %>% filter(Site != "Bu" & Site != "Zy")

class_df_long <- class_df_long %>% mutate(MGVRep = ifelse(Class == "NoMGVMatch", "Novel", "In MGV"))
class_df_long$MGVRep <- as.factor(class_df_long$MGVRep)
class_df_long$MGVRep <- relevel(class_df_long$MGVRep, 'Novel')
class_df_long <- merge(class_df_long, site_df, by="Site")

class_df_long_summary <- class_df_long %>% group_by(SiteName, MGVRep, Site) %>%
  summarise(n = sum(Count))

# this is where to get the value for the manuscript
class_df_long_summary %>% group_by(MGVRep) %>% summarise(total = sum(n))

viruses_compare <- ggplot(class_df_long_summary, aes(x=SiteName, y=n, fill=MGVRep)) + 
  geom_col() + 
  ylab("Viral Genomes") +
  #scale_fill_manual(values = c("#a3c3c9", "#45818e"))+ 
  scale_fill_manual(values = c("#cfcfcf", "#969696"))+ 
  scale_x_discrete(limits = c("Nanoro", "Navrongo", "DIMAMO", "Bushbuckridge", "Soweto", "Nairobi")) + 
  theme_bw()

viruses_compare

write.table(class_df_long_summary, "/Users/dylanmaghini/Downloads/sourcedata/fig3b.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

plottheme <- theme(
  panel.grid.major.x = element_blank(), 
  panel.grid.minor.x = element_blank(), 
  panel.grid.minor.y = element_blank(),
  legend.title = element_blank(), 
  axis.title.x = element_blank(), 
  axis.text.x = element_text(angle=45, hjust = 1),
  legend.position = c(1, 1),
  legend.justification = c("right", "top"),
  legend.box.just = "right",
  legend.margin = margin(1, 1, 1, 1),
  legend.background = element_rect(color = NA, fill='transparent'),
  legend.box.background = element_rect(color = NA, fill='transparent'),
  legend.key.size = unit(0.2, 'cm'),
  text = element_text(size=8), 
  title = element_text(size=6), 
  legend.text = element_text(size=4, margin = margin(l = -4, unit = "pt"))
)


plot_grid(bacteria_compare+plottheme, 
          viruses_compare + plottheme,
          proteins_compare+plottheme, 
          nrow=1, ncol=3, rel_widths=c(1,1.05,1))

##### RAREFACTIONS #####


metadata <- data.frame(site_abbreviation = c("AG", "BF", "KY", "SW", "DM", "GH"),
                       site_name = c("Bushbuckridge", "Nanoro", "Nairobi", "Soweto", "Dikgale", "Navrongo"), 
                       site_color = c("#889A4D","#BE9856","#4A4C59","#82354E","#DBBF5A","#99ADD0"))

pal <- c("#889A4D","#BE9856","#4A4C59","#82354E","#DBBF5A","#99ADD0", "#8A8A8A")
names(pal) <- c("AG", "BF", "KY", "SW", "DM", "GH", "All")

sitelabels <- c("Bushbuckridge", "Nanoro", "Nairobi", "Soweto", "DIMAMO", "Navrongo", "All")
names(sitelabels) <- c("AG", "BF", "KY", "SW", "DM", "GH", "All")

plot_curves <- function(dataframe, ytitle, xscale, yscale, title) {
  ggplot(dataframe, aes(x=NumParticipants, y=meanClusters, fill = Category, color=Category)) + 
    geom_ribbon(aes(ymin=meanClusters-sdClusters, ymax=meanClusters+sdClusters), alpha=0.3, color=NA) + 
    geom_line(size=1) + 
    theme_bw() + 
    scale_color_manual(values = pal, labels=sitelabels, breaks=c("All", "BF", "GH", "DM", "AG", "SW", "KY")) + 
    scale_fill_manual(values=pal, labels=sitelabels, breaks=c("All", "BF", "GH", "DM", "AG", "SW", "KY")) + 
    ylab(ytitle) + 
    xlab("Number of Participants") + 
    scale_x_continuous(limits=xscale) + 
    scale_y_continuous(limits=yscale)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          text = element_text(size=8), legend.title=element_blank(), 
          title = element_text(size=6), 
          legend.position = "bottom", 
          legend.key.size = unit(0.5, 'cm')) + 
    guides(colour = guide_legend(nrow = 1))
  
}


#########  PROTEIN RAREFACTION #####

infile_protein95 <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/proteins/231202_uhgp_compare/mmseqs_sample_clusters_95_concat.tsv", sep="\t", header=TRUE)

names(infile_protein95) <- c("Category", "NumParticipants", "Repetition", "NumClusters")

curve_protein95 <- infile_protein95 %>%
  group_by(NumParticipants, Category) %>%
  summarise(
    meanClusters = mean(NumClusters),
    sdClusters = sd(NumClusters))

proteins95_all <- plot_curves(curve_protein95, "Bacterial Proteins", c(0,1600), c(0,1e7), "Protein Clusters (95% AAI)")
proteins95_all
proteins95_zoomed <- plot_curves(curve_protein95, "Bacterial Proteins", c(0,500), c(0,6e6), "Protein Clusters (95% AAI)")
proteins95_zoomed

write.table(distinct(infile_protein95) %>% mutate(Site = ifelse(Category == "AG", "Agincourt", 
                                                     ifelse(Category == "BF", "Nanoro", 
                                                            ifelse(Category == "GH", "Navrongo", 
                                                                   ifelse(Category == "KY", "Nairobi", 
                                                                          ifelse(Category == "DM", "DIMAMO", ifelse(Category == "SW", "Soweto", "All"))))))), "/Users/dylanmaghini/Downloads/sourcedata/ed7f.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

#########  VIRUS RAREFACTION #####
infile_virus <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/phages/cluster_subsampling.tsv", sep="\t", header=TRUE)

curve_virus <- infile_virus %>%
  group_by(NumParticipants, Category) %>%
  summarise(
    meanClusters = mean(NumClusters),
    sdClusters = sd(NumClusters))

viruses_all <- plot_curves(curve_virus, "Viral Genomes", c(0,1600), c(0,40100), "Viral Clusters (95% ANI)")
viruses_zoomed <- plot_curves(curve_virus, "Viral Genomes", c(0,500), c(0,17000), "Viral Clusters (95% ANI)")

viruses_zoomed + scale_x_continuous(limits = c(0, 500), breaks = c(0, 250, 500)) + theme(legend.position = "none", 
                       axis.title = element_text(size = 12), 
                       axis.text = element_text(size = 12))
ggsave("/Users/dylanmaghini/scg4/projects/awigen2/00.plots/virus_rarefaction.pdf", dpi = 300, w = 2.5, h = 2)

write.table(infile_virus%>% mutate(Site = ifelse(Category == "AG", "Agincourt", 
                                                 ifelse(Category == "BF", "Nanoro", 
                                                        ifelse(Category == "GH", "Navrongo", 
                                                               ifelse(Category == "KY", "Nairobi", 
                                                                      ifelse(Category == "DM", "DIMAMO", ifelse(Category == "SW", "Soweto", "All"))))))), "/Users/dylanmaghini/Downloads/sourcedata/ed9c.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)
#########  BIN RAREFACTION #####
infile_bins <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/binning/01.dereplication/bin_sampling.tsv", sep="\t", header=TRUE)

curve_bins <- infile_bins %>%
  group_by(NumParticipants, Category) %>%
  summarise(
    meanClusters = mean(NumClusters),
    sdClusters = sd(NumClusters))

write.table(infile_bins %>% mutate(Site = ifelse(Category == "AG", "Agincourt", 
                                                 ifelse(Category == "BF", "Nanoro", 
                                                        ifelse(Category == "GH", "Navrongo", 
                                                               ifelse(Category == "KY", "Nairobi", 
                                                                      ifelse(Category == "DM", "DIMAMO", ifelse(Category == "SW", "Soweto", "All"))))))), 
            "/Users/dylanmaghini/Downloads/sourcedata/ed7e.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

bins_all <- plot_curves(curve_bins, "Bacterial Genomes", c(0,1600), c(0,2500), "Genome Clusters (95% ANI)")
bins_zoomed <- plot_curves(curve_bins, "Bacterial Genomes", c(0,500), c(0,1800), "Genome Clusters (95% ANI)")
bins_zoomed


###### NEW FEATURES PER SAMPLE #####

# BACTERIA
cdb_original <- read.csv("/Users/dylanmaghini/scg4/projects/awigen2/binning/01.dereplication/out/data_tables/Cdb.csv", header=TRUE, sep=",")
cdb_original <- cdb_original %>% select(genome, secondary_cluster)

wdb_original <- read.csv("/Users/dylanmaghini/scg4/projects/awigen2/binning/01.dereplication/out/data_tables/Wdb.csv", header=TRUE, sep=",")
wdb_original <- wdb_original %>% select(genome, cluster)
names(wdb_original) <- c("representative_genome", "secondary_cluster")

original_results <- merge(cdb_original, wdb_original, by=c("secondary_cluster"), all.x = TRUE)

cdb_uhgg <- cdb_awi %>% select(SiteCode, genome, ClusterType, Site)
names(cdb_uhgg) <- c("Sitecode", "representative_genome", "ClusterType", "Site")

bacteria <- merge(original_results, cdb_uhgg, by=c("representative_genome"), all.x = TRUE)
bacteria <- bacteria %>% mutate(Sample = gsub("_.*", "", genome))

bacteria_persample <- bacteria %>% group_by(Sample) %>% 
  summarise(n =sum(ClusterType == "Novel")) 
bacteria_persample <- bacteria_persample %>% mutate(Site = substr(Sample, 1, 2))


sitedf <- data.frame(sitename = c("Agincourt", "Nanoro", "Nairobi", "Soweto", "DIMAMO", "Navrongo"), 
                     Site = c("AG", "BF", "KY", "SW", "DM", "GH"))
sitepal <- c("#889A4D","#BE9856","#4A4C59","#82354E","#DBBF5A","#99ADD0")
names(sitepal) <- c("Agincourt", "Nanoro", "Nairobi", "Soweto", "DIMAMO", "Navrongo")

bacteria_persample <- merge(bacteria_persample, sitedf, by=c("Site"), all.x = TRUE)
bacteria_persample <- bacteria_persample %>% mutate(sitename = fct_relevel(sitename, "Nanoro", "Navrongo", "DIMAMO", "Agincourt", "Soweto", "Nairobi"))
bac_per <- ggplot(bacteria_persample, aes(x=sitename, y=n, fill=sitename)) + 
  geom_jitter(color="grey", width = 0.1, size = 0.3) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.9) + 
  ylab("Novel Bacterial Genomes") + 
  scale_fill_manual(values = sitepal) + 
  theme_bw() + 
  stat_compare_means(method = "kruskal.test", label.y = 15, label.x = 1.5, size = 2.5) + 
  theme(axis.title.x = element_blank(), 
        legend.position = "None", 
        axis.text.x = element_text(angle = 45, hjust =1), 
        panel.grid.major.x = element_blank())
bac_per

write.table(bacteria_persample, "/Users/dylanmaghini/Downloads/sourcedata/ed7b.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)
bacteria_persample %>% group_by(sitename) %>% 
  summarise(med = mean(n))


# VIRUSES
virus_awi_mgv <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/phages/per_phage_mgv_novelty.tsv", sep="\t", header = FALSE)
names(virus_awi_mgv) <- c("Virus", "InMGV")
virus_awi_mgv <- virus_awi_mgv %>% mutate(Sample = gsub("_.*", "", Virus))
virus_persample <- virus_awi_mgv %>% group_by(Sample) %>% 
  summarise(n = sum(InMGV == "False")) %>% filter(Sample != "Buffer4" & Sample != "Buffer8")
virus_persample <- virus_persample %>% mutate(Site = substr(Sample, 1, 2))

virus_persample <- merge(virus_persample, sitedf, by=c("Site"), all.x = TRUE)
virus_persample <- virus_persample %>% mutate(sitename = fct_relevel(sitename, "Nanoro", "Navrongo", "DIMAMO", "Agincourt", "Soweto", "Nairobi"))

vir_per <- ggplot(virus_persample, aes(x=sitename, y=n, fill=sitename)) + 
  geom_jitter(color="grey", width = 0.2, size = 0.5) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.9) + 
  ylab("Novel Viral Genomes") + 
  scale_fill_manual(values = sitepal) + 
  theme_bw() + 
  stat_compare_means(method = "kruskal.test", label.y = 160, label.x = 1.5, size = 2.5) + 
  theme(axis.title.x = element_blank(), 
        legend.position = "None", 
        axis.text.x = element_text(angle = 45, hjust =1), 
        panel.grid.major.x = element_blank())
vir_per
ggsave("/Users/dylanmaghini/scg4/projects/awigen2/00.plots/vir_per_indiv.pdf", dpi = 300, h = 2, w = 2)
virus_persample %>% group_by(sitename) %>% 
  summarise(med = mean(n))

write.table(virus_persample, "/Users/dylanmaghini/Downloads/sourcedata/ed9a.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# PROTEINS
protein_sample <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/proteins/231202_uhgp_compare/per_sample_counts.tsv", header=FALSE, sep="\t")
names(protein_sample) <- c("Sample", "Matching", "Nonmatching", "NotInDB")
protein_sample <- protein_sample %>% mutate(n = Nonmatching + NotInDB)
proteins_persample <-protein_sample %>% select(Sample, n) %>% mutate(Site = substr(Sample, 1, 2))

proteins_persample <- merge(proteins_persample, sitedf, by=c("Site"), all.x = TRUE)
proteins_persample <- proteins_persample %>% mutate(sitename = fct_relevel(sitename, "Nanoro", "Navrongo", "DIMAMO", "Agincourt", "Soweto", "Nairobi"))


proteins_per <- ggplot(proteins_persample, aes(x=sitename, y=n, fill=sitename)) + 
  geom_jitter(color="grey", width = 0.2, size = 0.5) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.9) + 
  ylab("Novel Proteins") + 
  ylim(0, 85000) + 
  scale_fill_manual(values = sitepal) + 
  theme_bw() + 
  stat_compare_means(method = "kruskal.test", label.y = 80000, label.x = 1.5, size = 2.5) + 
  theme(axis.title.x = element_blank(), 
        legend.position = "None", 
        axis.text.x = element_text(angle = 45, hjust =1), 
        panel.grid.major.x = element_blank(), 
        )
proteins_per

write.table(proteins_persample, "/Users/dylanmaghini/Downloads/sourcedata/ed7c.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

proteins_persample %>% group_by(sitename) %>% 
  summarise(med = mean(n))

plot_grid(bac_per, vir_per, proteins_per, nrow = 1, ncol = 3, labels = c("a", "b", "c"))
ggsave(here("00.plots/novelty_figure_per_sample.pdf"), dpi=300, h=2.7, w=6.9)
ggsave(here("00.plots/novelty_figure_per_sample.jpeg"), dpi=300, h=2.7, w=6.9)

# calculate n
bacteria_persample %>% group_by(Site) %>% summarise(count = n())


## redo for new supplementary figure ##
plot_grid(proteins_compare + theme(legend.position = "none", 
                                   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
                                   axis.title.x = element_blank(), 
                                   text = element_text(size = 12), 
                                   panel.grid.major.x = element_blank()) ,
          bac_per + theme(text = element_text(size = 12)), 
          proteins_per + theme(text = element_text(size = 12)),  
          bins_zoomed + theme(legend.position = "none", text = element_text(size = 12), axis.title = element_text(size = 12)), 
          proteins95_zoomed + theme(legend.position = "none", text = element_text(size = 12), axis.title = element_text(size = 12)),   
          align = "hv", axis = "lr")

ggsave("/Users/dylanmaghini/scg4/projects/awigen2/00.plots/novelty_supplement_new.pdf", dpi = 300, h = 5, w = 7)


##### ZOLFO NOVELTY #####

site_df <- data.frame(Site=c("AG", "BF", "DM", "GH", "KY", "SW"),
                      SiteName = c("Agincourt", "Nanoro", "DIMAMO", "Navrongo", "Nairobi", "Soweto"))

awi_all <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/phages/response/zolfo_compare/awi_zolfo_novelty_allphages.tsv", header = FALSE, sep = "\t")
names(awi_all) <- c("Phage", "Novelty")
awi_all <- awi_all %>% mutate(individual = substr(Phage, 1, 5))


virus_persample <- awi_all %>% group_by(individual) %>% 
  summarise(n = sum(Novelty == "Novel")) 
virus_persample <- virus_persample %>% mutate(Site = substr(individual, 1, 2))
virus_persample <- merge(virus_persample, site_df, by = "Site")
virus_persample <- virus_persample %>% mutate(SiteName = fct_relevel(SiteName, "Nanoro", "Navrongo", "DIMAMO", "Agincourt", "Soweto", "Nairobi"))


virus_persample %>% group_by(SiteName) %>% 
  summarise(median_new = median(n))
median(virus_persample$n)
awi_all_plot <- ggplot(virus_persample, aes(x=SiteName, y=n, fill=SiteName)) + 
  geom_jitter(color="grey", width = 0.2, size = 0.5) +  
  geom_boxplot(outlier.shape = NA, alpha = 0.9) + 
  ylab("Novel Viral Genomes") + 
  scale_fill_manual(values = sitepal) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        legend.position = "None", 
        axis.text.x = element_text(angle = 45, hjust =1, color = "black"), 
        axis.text.y = element_text(color = "black"), 
        panel.grid.major.x = element_blank(), 
        text = element_text(size=10, color = "black"), 
        panel.border = element_rect(color = "black", linewidth = 1)
  )
awi_all_plot

write.table(virus_persample, "/Users/dylanmaghini/Downloads/sourcedata/ed9b.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

ggsave(here("/Users/dylanmaghini/scg4/projects/awigen2/00.plots/novelty_figure_zolfo_individual.jpg"), dpi=300, h=2, w=2)
ggsave(here("/Users/dylanmaghini/scg4/projects/awigen2/00.plots/novelty_figure_zolfo_individual.pdf"), dpi=300, h=2, w=2)

