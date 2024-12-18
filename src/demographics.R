library(ggplot2) # load plotting library
library(dplyr)   # load data table manipulation library
library(forcats) # library for ordering factors
library(RColorBrewer)
library(ggpubr)
library(cowplot)

pheno_data <- read.table("/Users/dylanmaghini/Documents/Bhatt/AWIGen/metadata_combination/AWIGen_mb_metadata_full.tsv", sep = "\t", header=TRUE)
pheno_data <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/metadata/AWIGen_mb_metadata_full.tsv", sep = "\t", header=TRUE)


hiv_participants <- pheno_data %>% filter(HIVComparisons == 1) %>% filter(gene_site != 2)

hiv_participants <- hiv_participants %>% mutate(hiv_cat = ifelse(p2_hiv_status_calculated_2gener == 0 & (p2_hiv_arv_meds_now_gener != 1 | is.na(p2_hiv_arv_meds_now_gener)), "HIV-", 
                                                                 ifelse(p2_hiv_status_calculated_2gener == 0, "exclude",
                                                                        ifelse(p2_hiv_status_calculated_2gener == 1 & p2_hiv_arv_meds_now_gener == 1, "PLWH",
                                                                               ifelse(p2_hiv_status_calculated_2gener == 1 & p2_hiv_arv_meds_now_gener == 0, "ART-", "other")))))


hiv_med <- c("#7F7776","#007C92", "#8C1515")
site_df <- data.frame(gene_site = c(1, 3, 6), 
                      sitename = c("Agincourt", "Nairobi", "Soweto"))
hiv_participants <- merge(hiv_participants, site_df, by = c("gene_site"), all.x = TRUE)

table(hiv_participants$hiv_cat)
hiv_participants <- hiv_participants %>% filter(hiv_cat != "exclude")
hiv_participants$gene_site <- as.factor(hiv_participants$gene_site)

hiv_participants <- hiv_participants %>% mutate(sitename = fct_relevel(sitename, "Agincourt", "Soweto", "Nairobi"))


hiv_participants <- hiv_participants %>% filter(metaG_id != "KY049" & metaG_id != "KY204")

hiv_participants %>% group_by(SiteCode, hiv_cat) %>% 
  summarise(n = n())

write.csv(hiv_participants, "/Users/dylanmaghini/Downloads/HIVdata.csv", quote = FALSE, row.names = FALSE)


hiv_age <- ggplot(hiv_participants %>% filter(hiv_cat != "ART-"), aes(x=sitename, y=demo_age_at_collection, fill = hiv_cat)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), color = "#D3D3D3") + 
  geom_boxplot(outlier.shape = NA, alpha = 0.6, col = "black") + 
  ylab("Age (years)") +
  scale_fill_manual(values = hiv_med) + 
  theme_bw() + 
  stat_compare_means(method = "wilcox", aes(group = hiv_cat), label = "p.format") + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.title = element_blank(), 
        text = element_text(color = "black", size = 12), 
        panel.border = element_rect(color = "black", size = 1), 
        legend.position = c(0.8, 0.8))

hiv_bmi <- ggplot(hiv_participants %>% filter(hiv_cat != "ART-"), aes(x=sitename, y=anth_bmi_c, fill = hiv_cat)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), color = "#D3D3D3") + 
  geom_boxplot(outlier.shape = NA, alpha = 0.6, col = "black") + 
  ylab("BMI") +
  scale_fill_manual(values = hiv_med) + 
  theme_bw() + 
  stat_compare_means(method = "wilcox", aes(group = hiv_cat), label = "p.format") + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.title = element_blank(), 
        text = element_text(color = "black", size = 12), 
        panel.border = element_rect(color = "black", size = 1), 
        legend.position = c(0.8, 0.8))

hiv_whr <- ggplot(hiv_participants %>% filter(hiv_cat != "ART-", anth_waist_hip_ratio > 0), aes(x=sitename, y=anth_waist_hip_ratio, fill = hiv_cat)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), color = "#D3D3D3") + 
  geom_boxplot(outlier.shape = NA, alpha = 0.6, col = "black") + 
  ylab("Waist Hip Ratio") +
  scale_fill_manual(values = hiv_med) + 
  theme_bw() + 
  stat_compare_means(method = "wilcox", aes(group = hiv_cat), label = "p.format") + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.title = element_blank(), 
        text = element_text(color = "black", size = 12), 
        panel.border = element_rect(color = "black", size = 1), 
        legend.position = c(0.8, 0.8))

hiv_waist <- ggplot(hiv_participants %>% filter(hiv_cat != "ART-", anth_waist_circumf > 0), aes(x=sitename, y=anth_waist_circumf, fill = hiv_cat)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), color = "#D3D3D3") + 
  geom_boxplot(outlier.shape = NA, alpha = 0.6, col = "black") + 
  ylab("Waist Circumference") +
  scale_fill_manual(values = hiv_med) + 
  theme_bw() + 
  stat_compare_means(method = "wilcox", aes(group = hiv_cat), label = "p.format") + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.title = element_blank(), 
        text = element_text(color = "black", size = 12), 
        panel.border = element_rect(color = "black", size = 1), 
        legend.position = c(0.8, 0.8))

hiv_waist

hiv_cholesterol <- ggplot(hiv_participants %>% filter(hiv_cat != "ART-"), aes(x=sitename, y=lipids_cholesterol, fill = hiv_cat)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), color = "#D3D3D3") + 
  geom_boxplot(outlier.shape = NA, alpha = 0.6, col = "black") + 
  ylab("Cholesterol (mmol/L)") +
  scale_fill_manual(values = hiv_med) + 
  theme_bw() + 
  stat_compare_means(method = "wilcox", aes(group = hiv_cat), label = "p.format") + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.title = element_blank(), 
        text = element_text(color = "black", size = 12), 
        panel.border = element_rect(color = "black", size = 1), 
        legend.position = c(0.8, 0.8))


hiv_glucose <- ggplot(hiv_participants %>% filter(hiv_cat != "ART-"), aes(x=sitename, y=glucose, fill = hiv_cat)) + 
 geom_jitter(position = position_jitterdodge(jitter.width = 0.1), color = "#D3D3D3") + 
  geom_boxplot(outlier.shape = NA, alpha = 0.6, col = "black") + 
  ylab("Glucose (mmol/L)") +
  scale_fill_manual(values = hiv_med) + 
  theme_bw() + 
  stat_compare_means(method = "wilcox", aes(group = hiv_cat), label = "p.format") + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title.x = element_blank(), 
        legend.title = element_blank(), 
        text = element_text(color = "black", size = 12), 
        panel.border = element_rect(color = "black", size = 1), 
        legend.position = c(0.8, 0.8))

plot_grid(hiv_age, hiv_whr, hiv_cholesterol, hiv_glucose, nrow = 1, ncol = 4)
ggsave("/Users/dylanmaghini/scg4/projects/awigen2/00.plots/supplement_hiv_pheno.pdf", dpi = 300, w = 13, h = 3.5)

write.table(hiv_participants %>% filter(hiv_cat != "ART-") %>% select(sitename, demo_age_at_collection, anth_waist_hip_ratio, lipids_cholesterol, glucose, hiv_cat), 
            "/Users/dylanmaghini/Downloads/sourcedata/ed10g.tsv", sep = "\t", 
            quote = FALSE, row.names = FALSE)


