library(ggplot2)
library(dplyr)
library(ade4)
library(vegan)
library(forcats)
library(ggridges)
library(cowplot)
library(tidyr)
library(readr)
library(data.table)


#### PHYLOGEO STATS ####
# read in country-level and site-level geographic coordinates
coords <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/binning/05.treponema/04.trep_response/05.plots_and_stats/coordinates.txt", header = TRUE, sep="\t")
names(coords) <- c("country", "lat", "long")

# read in genome information
uhgg_information <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/binning/05.treponema/01.publicdata/uhgg_genomes-all_metadata_treponemaD.tsv", header = FALSE, sep = "\t")
uhgg_information <- uhgg_information %>% select(V1, V22)
names(uhgg_information) <- c("genome", "country")

#genome_list <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/binning/05.treponema/04.trep_response/genomes_list_tmp.txt", header = FALSE, sep = "\t")
genome_list <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/binning/05.treponema/04.trep_response/02.checkm/checkm_pass.csv", header = TRUE, sep = "\t")
genome_list <- genome_list %>% mutate(Bin.Id = paste(Bin.Id, ".fa", sep = ""))
names(genome_list) <- c("genome")

genome_list <- genome_list %>% mutate(genome = gsub(".fa", "", genome))
genome_list <- merge(genome_list, uhgg_information, by=c("genome"), all.x = TRUE)

genome_list <- genome_list %>% mutate(country = ifelse(substr(genome, 1, 2) == "AG", "Agincourt", 
                                                       ifelse(substr(genome, 1,2) == "BF", "Nanoro", 
                                                              ifelse(substr(genome, 1, 2) == "GH", "Navrongo", 
                                                                     ifelse(substr(genome, 1, 2) == "SW", "Soweto",
                                                                            ifelse(substr(genome, 1, 2) == "KY", "Nairobi", 
                                                                                   ifelse(substr(genome, 1, 2) == "DM", "DIMAMO", 
                                                                                          ifelse(substr(genome, 1, 3) == "REF", "Tanzania", country))))))))

genome_list <- genome_list %>% filter(!is.na(country)) %>% mutate(country = ifelse(country == "United Republic of Tanzania", "Tanzania", country))
genome_list <- merge(genome_list, coords, by = c("country"), all.x = TRUE)                                      
genome_list <- genome_list %>% mutate(genome = gsub("-", "_", gsub("\\.", "_", genome)))

# generate genome's country distance matrix
country.dists <- as.matrix(dist(cbind(genome_list$lon, genome_list$lat), method = "euclidean"))

# read in genetic distance matrix
gen.dist <- read.csv("/Users/dylanmaghini/scg4/projects/awigen2/binning/05.treponema/04.trep_response/05.plots_and_stats/tree_dist_all.csv", header = TRUE, row.names = 1, sep = ",")
rownames(gen.dist) <- gsub("-", "_", gsub("\\.", "_", gsub(" ", "_", rownames(gen.dist))))
colnames(gen.dist) <- gsub("\\.", "_", colnames(gen.dist))


#  Remove GCF column
column_index <- which(colnames(gen.dist) == "GCF_000195275_1_ASM19527v1_genomic")
modified_matrix <- gen.dist[, -column_index, drop = FALSE]
# Remove GCF row
row_index <- which(rownames(modified_matrix) == "GCF_000195275_1_ASM19527v1_genomic")
modified_matrix <- modified_matrix[-row_index, , drop = FALSE]

# will need to figure out how to reorder this 
desired_order <- genome_list$genome

# Get the indices for the desired row and column order
row_indices <- match(desired_order, rownames(modified_matrix))
col_indices <- match(desired_order, colnames(modified_matrix)) 

# Reorder rows and columns of the distance matrix
gen.dist.mat <- as.matrix(modified_matrix) 
gen.dist.ordered <- gen.dist.mat[row_indices, col_indices]
gen.dist.ordered <- as.dist(gen.dist.ordered)
country.dists <- as.dist(country.dists)

#mantel.rtest(country.dists, gen.dist.ordered, nrepet = 9999)

sort(rownames(as.matrix(gen.dist.ordered)))

#### PERMANOVA
# permanova in vegan for within-country distance 
#  versus between country distance using adonis2

# load geographic groupings
regions <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/binning/05.treponema/04.trep_response/05.plots_and_stats/region_info.txt", header = TRUE, sep="\t")
names(genome_list) <- c("Site", "genome", "lat", "long")
genome_list <- merge(genome_list, regions, by = c("Site"))

n_permutations <- 1000
adonis_site <- adonis(gen.dist.ordered ~ Site, 
                        data = genome_list, 
                        permutations = n_permutations)
adonis_country <- adonis(gen.dist.ordered ~ Country, 
                        data = genome_list, 
                        permutations = n_permutations)
adonis_region <- adonis(gen.dist.ordered ~ Region, 
                         data = genome_list, 
                         permutations = n_permutations)
adonis_continent <- adonis(gen.dist.ordered ~ Continent, 
                         data = genome_list, 
                         permutations = n_permutations)

adonis_site[["aov.tab"]][["Pr(>F)"]][1]
adonis_country[["aov.tab"]][["Pr(>F)"]][1]
adonis_region[["aov.tab"]][["Pr(>F)"]][1]
adonis_continent[["aov.tab"]][["Pr(>F)"]][1]


##### PERMANOVA, AWI-Gen only #####
awi_genomes <- genome_list %>% filter(Site %in% c("Nanoro", "Navrongo", "Soweto", "Nairobi", "Agincourt", "DIMAMO"))
awi_genomes <- awi_genomes %>% mutate(genome = gsub("\\.", "_", genome))
awi_gen_dist <- data.frame(as.matrix(gen.dist.ordered))
awi_gen_dist <- awi_gen_dist %>% select(all_of(awi_genomes$genome)) # filtering down to 75? 
awi_gen_dist <- awi_gen_dist[rownames(awi_gen_dist) %in% awi_genomes$genome, ]
awi_gen_dist <- as.dist(as.matrix(awi_gen_dist))

# genome_list_genomes <- awi_genomes$genome # BF148_maxbin.013_sub
# distmat_genomes <- colnames(awi_gen_dist) 

adonis_site <- adonis(awi_gen_dist ~ Region, 
                      data = awi_genomes, 
                      permutations = 100000)
adonis_site[["aov.tab"]][["Pr(>F)"]][1]


#### PLOTTING METADATA #### 

genome_list <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/binning/05.treponema/04.trep_response/02.checkm/checkm_pass.csv", header = TRUE, sep = "\t")
genome_list <- genome_list %>% mutate(Bin.Id = paste(Bin.Id, ".fa", sep = ""))
names(genome_list) <- c("genome")

genome_list <- genome_list %>% mutate(genome = gsub(".fa", "", genome))
genome_list <- merge(genome_list, uhgg_information, by=c("genome"), all.x = TRUE)

genome_list <- genome_list %>% mutate(country = ifelse(substr(genome, 1, 2) == "AG", "Agincourt", 
                                                       ifelse(substr(genome, 1,2) == "BF", "Nanoro", 
                                                              ifelse(substr(genome, 1, 2) == "GH", "Navrongo", 
                                                                     ifelse(substr(genome, 1, 2) == "SW", "Soweto",
                                                                            ifelse(substr(genome, 1, 2) == "KY", "Nairobi", 
                                                                                   ifelse(substr(genome, 1, 2) == "DM", "DIMAMO", 
                                                                                          ifelse(substr(genome, 1, 3) == "REF", "Tanzania", country))))))))

genome_list <- genome_list %>% filter(!is.na(country)) %>% mutate(country = ifelse(country == "United Republic of Tanzania", "Tanzania", country))
genome_list <- merge(genome_list, coords, by = c("country"), all.x = TRUE)                                      

regions <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/binning/05.treponema/04.trep_response/05.plots_and_stats/region_info.txt", header = TRUE, sep="\t")
names(genome_list) <- c("Site", "genome", "lat", "long")
genome_list <- merge(genome_list, regions, by = c("Site"))

### 11 countries
### 6 regions
### 4 continents

continent_pal <- data.frame(ColorContinent = c("#C03221", "#E0CA3C", "#4e597a", "#3F8170"), 
                            Continent = c("Americas", "Europe", "Africa", "Oceania"))

#continent_pal <- data.frame(ColorContinent = c("#3771ff", "#3ddc98", "#df2a35", "#fdca41",  "#ffffff"), 
#                            Continent = c("Africa", "Oceania", "South America", "Europe", NA))

region_pal <- data.frame(ColorRegion = c("#962820", "#C0AC27", "#3C4253", "#959CB2", "#67718E", "#326659"), 
                         Region = c("Americas", "Europe", "East Africa", "West Africa", "Southern Africa", "Oceania"))

#country_pal <- data.frame(ColorCountry = c("#C2D3FF", "#7099FF", "#1F5EFF","#003ACC", "#00237A", "#000C29", "#3ddc98", "#EA7179", "#C31D28", "#FEDE86", "#F2B202"), 
#                          Country = c("Tanzania", "Madagascar", "Ghana", "Burkina Faso", "Kenya", "South Africa", "Fiji", "Peru", "El Salvador", "Denmark", "Netherlands"))

# need to finish this - question is whether to do country by region colors or AWI colors, or 
# to not show country at all 
country_pal <- data.frame(ColorCountry = c("#889A4D","#BE9856","#4A4C59","#99ADD0", "#771911", "#6A376B", "#CF6A5D", "#B3442B", "#E29721", "#1D4080", "#074736"), 
                          Country = c("South Africa", "Burkina Faso", "Kenya", "Ghana", "Denmark", "El Salvador", "Fiji", "Netherlands", "Peru", "Tanzania", "Madagascar"))


genome_list <- merge(genome_list, continent_pal, by=c("Continent"), all.x = TRUE)
genome_list <- merge(genome_list, region_pal, by=c("Region"), all.x = TRUE)
genome_list <- merge(genome_list, country_pal, by=c("Country"), all.x = TRUE)

genome_list <- genome_list %>% mutate(Study = ifelse(Site %in% c("Agincourt", "Navrongo", "Nanoro", "DIMAMO", "Soweto", "Nanoro"), "AWI-Gen 2", "Public Genome"))
genome_list <- genome_list %>% mutate(StudyColor = ifelse(Study == "AWI-Gen 2", "#e68d85", "#ffffff"))


write.table(genome_list %>% select(genome, ColorContinent, Continent), "/Users/dylanmaghini/scg4/projects/awigen2/binning/05.treponema/04.trep_response/04.fasttree/1.fasttree_all/continent_metadata.tsv", sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(genome_list %>% select(genome, ColorCountry, Country), "/Users/dylanmaghini/scg4/projects/awigen2/binning/05.treponema/04.trep_response/04.fasttree/1.fasttree_all/country_metadata.tsv", sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(genome_list %>% select(genome, ColorRegion, Region), "/Users/dylanmaghini/scg4/projects/awigen2/binning/05.treponema/04.trep_response/04.fasttree/1.fasttree_all/region_metadata.tsv", sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
write.table(genome_list %>% select(genome, StudyColor, Study), "/Users/dylanmaghini/scg4/projects/awigen2/binning/05.treponema/04.trep_response/04.fasttree/1.fasttree_all/study_metadata.tsv", sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)

cat(paste(continent_pal$ColorContinent, collapse="\t"))
cat(paste(continent_pal$Continent, collapse="\t"))

cat(paste(country_pal$ColorCountry, collapse="\t"))
cat(paste(country_pal$Country, collapse="\t"))

cat(paste(region_pal$ColorRegion, collapse="\t"))
cat(paste(region_pal$Region, collapse="\t"))


##### GENOME STATS #####
genome_length <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/binning/05.treponema/02.trep_redo/05.genome_stats/genome_sizes.txt", header=FALSE, sep="\t")
names(genome_length) <- c("genome", "GenomeLength")
genome_length <- genome_length %>% mutate(genome = gsub(".fa", "", genome))

genome_info <- merge(genome_list, genome_length, by=c("genome"))

site_stats <- genome_info %>% filter(substr(genome, 1, 2) %in% c("AG", "GH", "BF", "KY", "DM", "SW")) %>% 
  mutate(Site = substr(genome, 1, 2)) %>%
  group_by(Site) %>% 
  summarise(count = n())

site_labels<- data.frame(Site = c("AG", "BF", "KY", "SW", "DM", "GH"), 
                         SiteName = c("Bushbuckridge", "Nanoro", "Nairobi", "Soweto", "DIMAMO", "Navrongo"))

site_stats <- merge(site_stats, site_labels, by=c("Site"), all.x = TRUE)
pal <- c("#889A4D","#BE9856","#4A4C59","#82354E","#DBBF5A","#99ADD0")
names(pal) <- c("AG", "BF", "KY", "SW", "DM", "GH")

site_stats <- site_stats %>%
  mutate(SiteName = fct_relevel(SiteName, 
                                "Nanoro", "Navrongo", "DIMAMO", 
                                "Bushbuckridge", "Soweto", "Nairobi"))

panela <- ggplot(site_stats, aes(x=SiteName, y=count, fill = Site)) + 
  geom_col() + 
  labs(y = expression(~italic("T. succinifaciens")~"MAGs")) + 
  scale_fill_manual(values=pal) + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(), 
        legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust=1), 
        axis.title.x = element_blank()) 
panela

write.table(site_stats, "/Users/dylanmaghini/Downloads/sourcedata/ed8a.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)
genome_info_awi <- genome_info %>% filter(substr(genome, 1, 2) %in% c("AG", "GH", "BF", "KY", "DM", "SW")) %>%
  mutate(Site = substr(genome, 1, 2))
genome_info_awi <- merge(genome_info_awi, site_labels, by=c("Site"), all.x = TRUE)

genome_info_awi <- genome_info_awi %>%
  mutate(SiteName = fct_relevel(SiteName, 
                                "Nanoro", "Navrongo", "DIMAMO", 
                                "Bushbuckridge", "Soweto", "Nairobi"))

genome_info_awi <- genome_info_awi %>% mutate(GenomeLength = GenomeLength / 1000000)

panelb <- ggplot(genome_info_awi %>% filter(SiteName != "Soweto"), aes(x = GenomeLength, y = SiteName, fill=Site, color=Site)) + 
  stat_density_ridges(alpha=0.8, scale=2, quantile_lines = TRUE, quantiles = 0.5) + 
  scale_fill_manual(values = pal) + 
  scale_color_manual(values=pal) +
  scale_y_discrete(expand = expansion(add = c(.3, 2.1))) + 
  xlab("MAG Length (Mbp)") +
  theme_bw() + 
  theme(legend.position = "none", 
        axis.title.y = element_blank(), 
        panel.grid.major.x = element_blank())
panelb

write.table(genome_info_awi %>% mutate(Site2 = ifelse(SiteName == "Bushbuckridge", "Agincourt", SiteName)), 
            "/Users/dylanmaghini/Downloads/sourcedata/ed8b.tsv", sep = "\t", quote = FALSE, row.names = FALSE )
# details for paper
mean(genome_info_awi$GenomeLength)
sd(genome_info_awi$GenomeLength)

#gene_pres <- read.csv("/Users/dylanmaghini/Downloads/240423_offline/gene_pres_abs_trimmed.csv", header = TRUE)
gene_counts <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/binning/05.treponema/04.trep_response/05.plots_and_stats/gene_cat_counts.tsv", header = FALSE, sep = "\t")
names(gene_counts) <- c("genome", "Cloud", "Shell", "Core")
gene_counts_long <- melt(gene_counts, id.vars = c("genome"), variable.name = "Category")

gene_counts_long <- gene_counts_long %>%
  mutate(Category = fct_relevel(Category, 
                                "Cloud", "Shell", "Core"))

ordering <- gene_counts_long %>% 
  filter(Category == "Cloud") %>%
  arrange(desc(value))

gene_counts_long$genome <- factor(gene_counts_long$genome, levels=ordering$genome)
panelc <- ggplot(gene_counts_long, aes(x=genome, y=value, fill=Category)) + 
  geom_col() + 
  scale_fill_manual(values = c("#cccccc","#8c8c8c","#4d4d4d")) + 
  scale_y_continuous(expand=c(0,0))+
  xlab("Genome") + 
  ylab("Number of Genes") + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.title = element_blank(), 
        legend.position = c(0.8,0.25), 
        legend.margin = margin(c(0.5,0.5,0.5,0.5)), 
        legend.spacing.x = unit(1, "mm"), 
        legend.spacing.y = unit(1, "mm"), 
        legend.background = element_rect(color = "black", fill = "white", linetype="solid"))
panelc

write.table(gene_counts_long, "/Users/dylanmaghini/Downloads/sourcedata/ed8c.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_grid(panela, panelb, panelc, nrow = 1, align = "h")
ggsave("/Users/dylanmaghini/scg4/projects/awigen2/00.plots/treponema_upper.pdf", dpi = 300, w = 8, h = 2.8)

gene_counts_long_core <- gene_counts_long %>% filter(Category == "Core")

gene_total <- gene_counts_long %>% group_by(genome) %>% 
  summarise(total = sum(value))
gene_total <- merge(gene_total, gene_counts_long_core, by="genome")
gene_total <- gene_total %>% mutate(PercentCore = value/total * 100)
mean(gene_total$PercentCore)
sd(gene_total$PercentCore)
