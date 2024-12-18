library(ggplot2)
library(dplyr)
library(ggpubr)
library(stringr)
library(purrr)
library(ggbeeswarm)
library(forcats)
library(cowplot)
library(tidyr)
library(paletteer)

# read in classification data
load('/Users/dylanmaghini/Downloads/tempclass/all_classification_tables.RData')
species <- lst.motus$species
species_rel <- prop.table(as.matrix(species), 2)
species_rel <- species_rel[rowSums(species_rel==0)!=ncol(species_rel),]
species_rel <- species_rel[!str_detect(rownames(species_rel), 'unassigned'),]
species_rel <- species_rel[, !str_detect(colnames(species_rel), 'Bu|Zy')]
species_rel_filt <- species_rel[str_detect(rownames(species_rel), 'Treponema'),]
species_rel_filt <- log10(species_rel_filt + 1e-04)

# read metadata 
load("/Users/dylanmaghini/Downloads/metadata_clean.RData")

species_rel_filt <- data.frame(t(species_rel_filt))
species_rel_filt$metaG_id <- rownames(species_rel_filt)

merge_df <- merge(metadata.clean, species_rel_filt, by.y = "metaG_id", 
                  by.x = "general_sample_id") %>% 
  mutate(Tsucc = Treponema.succinifaciens..ref_mOTU_v3_11997. > -4)

merge_df %>% 
  group_by(general_site, Tsucc) %>% 
  tally() %>% 
  group_by(general_site) %>% 
  mutate(n.all=sum(n)) %>% 
  mutate(freq=n/n.all) %>% 
  filter(Tsucc)


merge_df_filt <- merge_df %>% 
  filter(general_site %in% c('Nanoro', 'Navrongo'))

fit <- lmerTest::lmer(Treponema.succinifaciens..ref_mOTU_v3_11997. ~ anthropometric_bmi + (1|general_site), merge_df)
summary(fit)



pheno_names <- c("microbiome_antibiotics", "anthropometric_bmi", "anthropometric_waist_hip_ratio", "anthropometric_waist_circumference",
                 "anthropometric_hip_circumference", "household_cattle", "household_other_livestock", "household_poultry", "card_hypertension_status",
                 "lipids_ldl", "lipids_hdl", "lipids_cholesterol", "lipids_triglycerides", "lab_insulin",
                 "lab_glucose", "demographic_age")

x <- 'Treponema.succinifaciens..ref_mOTU_v3_11997.'

results_table <- tibble(
  species = character(),
  pheno_name = character(),
  estimate = numeric(),
  p_value = numeric()
)

for (i in 1:length(pheno_names)) {
  message(i)
  x <- 'Treponema.succinifaciens..ref_mOTU_v3_11997.'
  y <- pheno_names[i]
  # fit <- lmerTest::lmer(formula=paste0(x, ' ~ ', y, 
  #                                      ' + (1|general_site) + (1|microbiome_antibiotics)'), 
  #                       data=merge_df)
  # res <-  coefficients(summary(fit))
  # results_table <- results_table %>% 
  #   add_row(
  #     species = paste0(x, '_abx'),
  #     pheno_name = y,
  #     estimate = res[2,1],
  #     p_value = res[2,5]
  #   )
  
  fit <- lmerTest::lmer(formula=paste0('Tsucc ~ ', y, 
                                       ' + (1|general_site)'), 
                        data=merge_df_filt)
  res <-  coefficients(summary(fit))
  results_table <- results_table %>% 
    add_row(
      species = "presence",
      pheno_name = y,
      estimate = res[2,1],
      p_value = res[2,5]
    )
  if (y != 'microbiome_antibiotics'){
    fit <- lmerTest::lmer(formula=paste0(
      'Tsucc ~ ', y, ' + (1|general_site) + (1|microbiome_antibiotics)'), 
                          data=merge_df_filt)
    res <-  coefficients(summary(fit))
    results_table <- results_table %>% 
      add_row(
        species = "presence_abx",
        pheno_name = y,
        estimate = res[2,1],
        p_value = res[2,5]
      )
  }
}

results_table %>% 
  group_by(species) %>% 
  mutate(q_value=p.adjust(p_value, method='fdr')) %>% 
  ungroup() %>% View


results_table %>% 
  select(pheno_name, species, p_value) %>% 
  mutate(p_value=-log10(p_value)) %>% 
  pivot_wider(names_from=species, values_from=p_value) %>% 
  ggplot(aes(x=presence, y=presence_abx)) + 
      geom_point() + 
    geom_abline(slope = 1, intercept = 0)

circum_stat <- data.frame(group1 = c("FALSE"), 
                          group2 = c("TRUE"), 
                          p = c("0.012"), 
                          y.position = 1)

hc <- merge_df_filt %>% 
  ggplot(aes(x=Tsucc, y=anthropometric_hip_circumference)) + 
    geom_beeswarm(color = "grey") +
    geom_boxplot(alpha=0.8, outlier.shape = NA, width = 0.4) + 
  ylab("Hip Circumference (cm)") + 
  ylim(70, 126) +
  xlab(expression(italic("T. succinifaciens"))) + 
  scale_x_discrete(labels = c("Absent", "Present")) + 
  stat_pvalue_manual(circum_stat, y.position = 125) +
    theme_bw() +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major.x = element_blank(), 
          panel.border = element_rect(color = "black", linewidth = 1),
          text = element_text(color = "black"), 
          axis.text = element_text(color = "black")) +
  #facet_wrap(~general_site) + 
  NULL

hc

write.table(merge_df_filt, "/Users/dylanmaghini/Downloads/sourcedata/fig2e.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

table(merge_df$Tsucc, merge_df$microbiome_antibiotics)

abx_stat <- data.frame(group1 = c("FALSE"), 
                          group2 = c("TRUE"), 
                          p = c("0.0075"), 
                          y.position = 1)

pal <- paletteer_d("Redmonder::sPBIBu")
#pal <- c(pal[2:9], "#000000")
abx_plot <- merge_df_filt %>% 
  group_by(Tsucc, general_site, microbiome_antibiotics) %>% 
  tally() %>% 
  group_by(Tsucc) %>% 
  mutate(n.all=sum(n)) %>% 
  mutate(freq=n/n.all) %>% 
  ggplot(aes(x=Tsucc, y=freq, fill=microbiome_antibiotics)) + 
    geom_bar(stat='identity') + 
  xlab(expression(italic("T. succinifaciens"))) + 
  ylab("Proportion of Individuals") + 
  scale_x_discrete(labels = c("Absent", "Present")) + 
  scale_fill_manual(values = rev(pal), na.value = "grey") + 
  #scale_fill_brewer(palette = "Blues", na.value = "grey", direction = -1) +
  guides(fill=guide_legend(title="Antibiotic Usage")) + 
  stat_pvalue_manual(abx_stat, y.position = 1.05, inherit.aes = FALSE) +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.border = element_rect(color = "black", linewidth = 1),
        text = element_text(color = "black"), 
        axis.text = element_text(color = "black"))
abx_plot

abx_leg <- get_legend(abx_plot)
plot_grid(hc, abx_plot + theme(legend.position = "none"), nrow = 1, ncol = 2, 
          rel_widths = c(1, 1))

ggsave("/Users/dylanmaghini/scg4/projects/awigen2/binning/05.treponema/04.trep_response/06.metadata_association/plots.pdf", dpi = 300, width = 3, h = 2.4)





#### looking at abx by site ####
metadata.clean %>% 
  group_by(general_site, microbiome_antibiotics) %>% 
  tally() %>% 
  group_by(general_site) %>% 
  mutate(n.all=sum(n)) %>% 
  mutate(freq=n/n.all) %>% 
  ggplot(aes(x=general_site, y=freq, fill=microbiome_antibiotics)) + 
  geom_bar(stat='identity') + 
  ylab("Proportion of Individuals") + 
  scale_fill_brewer(palette = "Blues", na.value = "grey", direction = -1) +
  guides(fill=guide_legend(title="Antibiotic Usage")) + 
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.border = element_rect(color = "black", linewidth = 1),
        text = element_text(color = "black"), 
        axis.text = element_text(color = "black"))

dev.off()
