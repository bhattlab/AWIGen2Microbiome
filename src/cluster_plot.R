library(ggplot2)
library(dplyr)
library(data.table)
library(here)
library(reshape2)
library(factoextra)
library(ComplexHeatmap)
library(circlize)
library(ggpmisc)
library(ggbeeswarm)
library(forcats)
library(cowplot)
library(tidyr)


cluster_file <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/phages/response/cluster_dist.tsv", header = FALSE, sep = "\t")

names(cluster_file) <- c("Representative", "Agincourt", "Nanoro", "DIMAMO", "Navrongo", "Nairobi", "Soweto", "Total")

cluster_file <- cluster_file %>% mutate(sites = rowSums(select(., c("Agincourt", "Agincourt", "Nanoro", "DIMAMO", "Navrongo", "Nairobi", "Soweto")) > 0))
cluster_file$sites <- as.factor(cluster_file$sites)

cluster_file <- cluster_file %>% mutate(percent = Total / 1820 * 100)

# calculate statistics
cluster_file %>% filter(percent >= 1) %>% nrow()
343/44506*100

cluster_file %>% filter(percent >= 5) %>% nrow()
38/44506*100

cluster_file %>% filter(Total >= 9) %>% nrow()

# Prevalence heatmap

# df of phages with global prevalence of ~1% or greater
prevalence_df <- cluster_file %>% select(!(sites)) %>% select(!(percent))

# n
# Agincort: 533
# DIMAMO: 203
# Nanoro: 384
# Navrongo: 235
# Nairobi: 239
# Soweto: 226

prevalence_df <- prevalence_df %>% mutate(Agincourt = Agincourt / 533)
prevalence_df <- prevalence_df %>% mutate(DIMAMO = DIMAMO / 203)
prevalence_df <- prevalence_df %>% mutate(Nanoro = Nanoro / 384)
prevalence_df <- prevalence_df %>% mutate(Navrongo = Navrongo / 235)
prevalence_df <- prevalence_df %>% mutate(Nairobi = Nairobi / 239)
prevalence_df <- prevalence_df %>% mutate(Soweto = Soweto / 226)
prevalence_df <- prevalence_df %>% mutate(Total = Total / 1820)

# site level prevalence statistics
prevalence_df <- prevalence_df %>% rowwise() %>% 
  mutate(maxsitevalue = max(Agincourt, DIMAMO, Nanoro, Navrongo, Nairobi, Soweto))

# calculate statistics
prevalence_df %>% filter(maxsitevalue >= 0.01) %>% nrow()
2071/44506*100

prevalence_df %>% filter(maxsitevalue >= 0.05) %>% nrow()
142/44506*100

##### HEATMAP ###

prevalence_heatmap_df <- prevalence_df %>% filter(Total >= 0.0098) %>% select(!(Total))
rownames(prevalence_heatmap_df) <- prevalence_heatmap_df$Representative
phage_beforecluster <- rownames(prevalence_heatmap_df)
prevalence_heatmap_df <- prevalence_heatmap_df[,-1]
prevalence_heatmap_df <- prevalence_heatmap_df %>% select(Nanoro, Navrongo, DIMAMO, Agincourt, Soweto, Nairobi)

## Phage novelty
novelty_df <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/phages/per_phage_mgv_novelty.tsv", header = FALSE, sep = "\t")
names(novelty_df) <- c("Representative", "inMGV")
prevalence_df <- merge(prevalence_df, novelty_df, by=c("Representative"), all.x = TRUE)

novel_prev <- prevalence_df %>% filter(inMGV == "False") 
novel_highprev <- novel_prev %>% filter(Total >= 0.01) %>% select(Representative)
write.table(novel_highprev, "/Users/dylanmaghini/scg4/projects/awigen2/phages/response/novel_list.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

novelty_heatmap <- novelty_df
rownames(novelty_heatmap) <- novelty_heatmap$Representative
novelty_heatmap <- novelty_heatmap %>% filter(Representative %in% phage_beforecluster)
novelty_heatmap <- novelty_heatmap %>% mutate(novelty_shade = ifelse(inMGV == "True", 0, 1))
# reorder novelty heatmap according to phage order in the prevalence heatmap
novelty_heatmap$Representative <- factor(novelty_heatmap$Representative, levels = phage_beforecluster)
novelty_heatmap <- novelty_heatmap[order(novelty_heatmap$Representative), ]

col_fun = colorRamp2(c(0, 0.31), c("white", "#404040"))
ht_opt(legend_border = "black")

# crassphage markers

# read in representatives that are crass phages
awi_crass <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/phages/response/crassphage/awi_crass_phages_list.tsv", header = FALSE, sep = "\t")
names(awi_crass) <- c("Representative")
awi_crass <- awi_crass %>% mutate(crassphage = "True")

prevalence_inc_crass <- merge(prevalence_df, awi_crass, by=c("Representative"), all.x = TRUE)
prevalence_inc_crass <- prevalence_inc_crass %>% mutate(crassphage = ifelse(is.na(crassphage), "False", "True"))

prevalence_inc_crass <- prevalence_inc_crass %>% filter(Representative %in% phage_beforecluster)
prevalence_inc_crass$Representative <- factor(prevalence_inc_crass$Representative, levels = phage_beforecluster)
prevalence_inc_crass <- prevalence_inc_crass[order(prevalence_inc_crass$Representative), ]
rownames(prevalence_inc_crass) <- prevalence_inc_crass$Representative

write.table(prevalence_inc_crass, "/Users/dylanmaghini/Downloads/sourcedata/fig3a.tsv", sep = "\t", 
            row.names = FALSE, quote = FALSE)
column_ha = HeatmapAnnotation(novelty = as.character(novelty_heatmap$novelty_shade), 
                              crassphage = as.character(prevalence_inc_crass$crassphage),
                              col = list(novelty = c("0" = "#ffffff", "1" = "#001663"), 
                                         crassphage = c("False" = "#ffffff", "True" = "#800800")),
                              simple_anno_size = unit(0.3, "cm"))

ht <- Heatmap(t(as.matrix(prevalence_heatmap_df)), show_row_names = TRUE, 
              show_column_names = FALSE, col = col_fun, border_gp = gpar(col = "black"), 
              cluster_rows = FALSE, column_dend_reorder = FALSE, column_split = 3, 
              row_names_side = "left", column_title = NULL, column_title_side = "bottom", 
              clustering_method_columns = "ward.D2", name = "Prevalence", 
              bottom_annotation = column_ha, column_names_gp = grid::gpar(fontsize = 4), 
              height = unit(2, "in"))

pdf(file="/Users/dylanmaghini/scg4/projects/awigen2/00.plots/phage_prevalence.pdf", width = 10, height = 3)
draw(ht, column_title = "Phage genomes", column_title_side = "bottom")
dev.off()

#### COVERM BIG PHAGES ####
big_phage_coverage <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/phages/response/bigphages/coverm/phage_abundances.tsv", header = FALSE, sep = "\t")
names(big_phage_coverage) <- c("Sample", "Genome", "Coverage")

newest_phages <- c("SW076_k141_174947__flag=0__multi=20.4459__len=368144_fragment_2", "GH039_k141_160891__flag=1__multi=37.0000__len=437741",
                   "BF036_k141_34448__flag=1__multi=25.3869__len=278266", "AG414_k141_427358__flag=3__multi=14.0001__len=249539", 
                   "AG493_k141_291015__flag=3__multi=56.0002__len=274949", "BF152_k141_381896__flag=1__multi=15.9945__len=351181", 
                   "BF122_k141_112032__flag=0__multi=18.0000__len=276716", "BF270_k141_594336__flag=3__multi=21.7851__len=252089", 
                   "BF087_k141_217560__flag=1__multi=23.6720__len=241270", "BF107_k141_63636__flag=3__multi=53.0000__len=237438")

big_phage_coverage <- big_phage_coverage %>% filter(Genome %in% newest_phages)

big_phage_coverage <- big_phage_coverage %>% mutate(Genome = gsub("__flag.*", "", Genome))
big_phage_coverage <- big_phage_coverage %>% mutate(site = substr(Sample, 1, 2))
big_phage_coverage <- big_phage_coverage %>% filter(site != "Bu" & site != "Zy")

big_phage_prevalence <- big_phage_coverage %>% 
  group_by(site, Genome) %>% 
  summarise(percentage = mean(Coverage > 0) * 100)

pal <- c("#889A4D","#BE9856","#4A4C59","#82354E","#DBBF5A","#99ADD0")
names(pal) <- c("AG", "BF", "KY", "SW", "DM", "GH")

big_phage_prevalence$site <- as.factor(big_phage_prevalence$site)
big_phage_prevalence <- big_phage_prevalence %>% mutate(site = fct_relevel(site, "BF", "GH", "DM", "AG", "SW", "KY"))

big_phage_prevalence_higherpercent <- big_phage_coverage %>% 
  group_by(site, Genome) %>% 
  summarise(percentage = mean(Coverage > 0.1) * 100)

phageids <- data.frame(Phage = c("SW076", "BF107", "BF087", "AG414", "BF270", "AG493", "BF122", "BF036", "BF152"), 
                           Order = c("Phage A", "Phage B", "Phage C", "Phage D", "Phage E", "Phage F", "Phage G", "Phage H", "Phage I"))

big_phage_prevalence_higherpercent <- big_phage_prevalence_higherpercent %>% 
  mutate(Genome = gsub("_.*", "", Genome))

big_phage_prevalence_higherpercent <- merge(big_phage_prevalence_higherpercent, phageids, 
                                            by.x = "Genome", by.y = "Phage")
big_phage_prevalence_higherpercent <- big_phage_prevalence_higherpercent %>% 
  mutate(site = fct_relevel(site, "BF", "GH", "DM", "AG", "SW", "KY"))
ggplot(big_phage_prevalence_higherpercent, aes(x = Order, y = percentage, fill = site)) + 
  geom_col(position = "dodge") + 
  theme_bw() + 
  ylab("Prevalence (%)") +
  scale_fill_manual(values = pal) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.title.x = element_blank(), 
        legend.position = "none")
write.table(big_phage_prevalence_higherpercent %>% mutate(Site = ifelse(site == "AG", "Agincourt", 
                                                                        ifelse(site == "SW", "Soweto", 
                                                                               ifelse(site == "GH", "Navrongo", 
                                                                                      ifelse(site == "BF", "Nanoro", 
                                                                                             ifelse(site == "KY", "Nairobi", 
                                                                                                  "DIMAMO")))))), "/Users/dylanmaghini/Downloads/sourcedata/ed9f.tsv", sep = "\t", 
            row.names = FALSE, quote = FALSE)
ggsave("/Users/dylanmaghini/scg4/projects/awigen2/phages/response/bigphages/jumbophage_prevalence.pdf", dpi = 300, w = 4, h = 2)

#### CRASS ####

# note, these are fraction abundances (out of 1, not out of 100)
phanta <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/phages/response/crassphage/phanta_crass_relative_abundance.txt", header = TRUE, sep = "\t")

phanta_crasslike  <- phanta %>% 
  filter(Taxon_Lineage_with_Names == "root_root|superkingdom_Viruses|clade_Duplodnaviria|kingdom_Heunggongvirae|phylum_Uroviricota|class_Caudoviricetes|order_Caudovirales|family_Podoviridae|no rank_unclassified Podoviridae|clade_crAss-like viruses" | 
           Taxon_Lineage_with_Names == "root_root|superkingdom_Viruses|clade_Duplodnaviria|kingdom_Heunggongvirae|phylum_Uroviricota|class_Caudoviricetes|order_Caudovirales|family_Podoviridae|no rank_unclassified Podoviridae|clade_crAss-like viruses|no rank_environmental samples|species_uncultured crAssphage")

phanta_crasslike <- melt(setDT(phanta_crasslike), id.vars = c("Taxon_Lineage_with_Names","Taxon_Lineage_with_IDs"), variable.name = "Sample")
phanta_crasslike <- phanta_crasslike %>% mutate(SiteCode = substr(Sample, 1, 2)) %>% filter(SiteCode != "Bu", SiteCode != "Zy")
phanta_crasslike <- phanta_crasslike %>% mutate(value = value * 100) %>%
  mutate(SiteCode = fct_relevel(SiteCode, "BF", "GH", "DM", "AG", "SW", "KY"))

pal <- c("#889A4D","#BE9856","#4A4C59","#82354E","#DBBF5A","#99ADD0")
names(pal) <- c("Agincourt", "Nanoro", "Nairobi", "Soweto", "DIMAMO", "Navrongo")

sitedf <- data.frame(SiteCode = c("BF", "GH", "DM", "AG", "SW", "KY"), 
                     Site = c("Nanoro", "Navrongo", "DIMAMO", "Agincourt", "Soweto", "Nairobi"))

phanta_crasslike <- phanta_crasslike %>% mutate(Taxon_Lineage_with_Names = gsub(".*\\|", "", Taxon_Lineage_with_Names)) %>% 
  mutate(Taxon_Lineage_with_IDs = gsub(".*\\|", "", Taxon_Lineage_with_IDs))

phanta_crasslike <- merge(phanta_crasslike, sitedf, by = "SiteCode")
phanta_crasslike <- phanta_crasslike %>% mutate(Taxon_Lineage_with_Names = ifelse(Taxon_Lineage_with_Names == "clade_crAss-like viruses", "crAss-like viruses", "crAssphage"))

# calculate prevalence
phanta_crasslike_prevalence <- phanta_crasslike %>% 
  group_by(Site, Taxon_Lineage_with_Names) %>% 
  summarise(percent_present = mean(value > 0))

phanta_crasslike_prevalence$Site <- as.factor(phanta_crasslike_prevalence$Site)
phanta_crasslike_prevalence <- phanta_crasslike_prevalence %>% mutate(Site = fct_relevel(Site, "Nanoro", "Navrongo", "DIMAMO", "Agincourt", "Soweto", "Nairobi"))

ggplot(phanta_crasslike_prevalence, aes(x = Site, y = percent_present, fill = Site)) + 
  geom_col() + 
  ylim(0, 1) + 
  ylab("Prevalence") + 
  facet_grid(~Taxon_Lineage_with_Names) + 
  scale_fill_manual(values = pal) + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), 
        axis.text = element_text(angle = 45, hjust = 1),  
        panel.grid.major.x = element_blank(), 
        legend.position = "none")

write.table(phanta_crasslike_prevalence, "/Users/dylanmaghini/Downloads/sourcedata/fig3d.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)
