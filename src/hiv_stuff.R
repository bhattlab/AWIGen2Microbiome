# ##############################################################################
#
## HIV analyses
#
# ##############################################################################

library("tidyverse")
library("here")
library("vegan")
library("labdsv")
library("lmerTest")
library("yaml")

source(here('src', 'utils.R'))

plot.params <- yaml.load_file(here('files', 'params.yml'))
site.dict <- unlist(plot.params$sites)
site.colours <- unlist(plot.params$sites_colours)

load(here('data', 'metadata', 'metadata_clean.RData'))
# load(here('data', 'classification', 'all_classification_tables.RData'))
load('./data/classification/new_db_tables.RData')

feat.red <- motus.gtdb.lvls$motus

meta.hiv <- metadata.clean %>% 
  filter(meta_hiv_comparison) %>% 
  filter(general_site !='Dimamo') %>% 
  filter(general_sample_id %in% colnames(feat.red)) %>% 
  filter(hiv_status %in% c("Yes", 'No')) %>%
  mutate(hiv_detail=case_when(
    hiv_status=='No' & hiv_medication=='Yes'~'exclude',
    hiv_status=='No' ~ 'HIV Negative(Control)',
    hiv_status=='Yes' & hiv_medication=='Yes' ~ 'HIV Positive(ART+)',
    hiv_status=='Yes' & hiv_medication=='No' ~ 'HIV Positive(ART-)',
    TRUE ~NA_character_)) %>% 
  filter(!is.na(hiv_detail)) %>% 
  filter(hiv_detail!='exclude')
feat.red <- feat.red[,meta.hiv$general_sample_id]

hiv_med <- c("#7F7776","#007C92", "#8C1515")

# export metadata for Dylan & Luicer
# meta.hiv %>% 
#  select(general_sample_id, general_study_id, hiv_status, 
#         hiv_medication, hiv_detail) %>% 
#  write_tsv(file='./files/hiv_metadata_red.tsv')

# ##############################################################################
# number of samples
g <- meta.hiv %>% 
  group_by(general_site, hiv_detail) %>% tally() %>% 
  ggplot(aes(x=general_site, y=n, fill=hiv_detail)) + 
  geom_bar(stat='identity', position = position_dodge()) + 
  scale_fill_manual(values=hiv_med) + 
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major.y = element_blank()) + 
  ylab('Number of samples') + xlab('') +
  coord_flip()
g$data %>% 
  write_csv('./files/source_data/Fig4a.csv')
ggsave(g, filename=here('figures','hiv','number_of_samples.pdf'), 
       width = 5, height = 3, useDingbats=FALSE)

# ##############################################################################
# Alpha
seed <- 1933
alpha <- .f_alpha(feat.red, l=7500)
df.alpha <- alpha$alpha_df %>% 
  rename(general_sample_id=Sample_ID) %>% 
  full_join(meta.hiv, by='general_sample_id')

g <- df.alpha %>% 
  mutate(general_site='all') %>% 
  mutate(facet='all') %>% 
  bind_rows(df.alpha %>% 
              mutate(facet='sites')) %>% 
  filter(hiv_detail!='HIV Positive(ART-)') %>% 
  mutate(general_site=factor(general_site, 
                             levels=c('all', 'Agincourt', 
                                      'Soweto', 'Nairobi'))) %>% 
  ggplot(aes(x=general_site, y=inv_simpson, fill=hiv_detail)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), 
              colour='#D3D3D3') + 
  geom_boxplot(alpha=0.6, outlier.shape = NA, col='black') + 
  facet_grid(~facet, scales = 'free', space='free') + 
  scale_fill_manual(values=hiv_med, name='') +
  xlab('') + ylab("Alpha diversity") +
  theme_bw() + theme(panel.grid.minor = element_blank(), 
                     panel.grid.major.x=element_blank(),
                     strip.background = element_blank(), 
                     strip.text = element_blank()) +
  NULL
g$data %>% 
  select(general_site, inv_simpson, general_sample_id, hiv_detail) %>% 
  write_csv('./files/source_data/Fig4b.csv')
ggsave(g, filename=here('figures','hiv','hiv_alpha_newdb.pdf'), 
       width = 10, height = 5, useDingbats=FALSE)

# testing
# individual
df.alpha %>% 
  group_by(general_site) %>% 
  bind_rows(df.alpha %>% 
              mutate(facet='sites')) %>% 
  filter(hiv_detail!='HIV Positive(ART-)') %>% 
  group_map(.f=function(.x, .y){
    res <- coefficients(summary(lm(data=.x, inv_simpson~hiv_status)))
    tibble(p.val=res[2,4], .y)
    }) %>% bind_rows()

#' # A tibble: 3 × 2
#' p.val general_site
#' <dbl> <fct>       
#' 1 0.00119  Agincourt  
#' 2 0.000153 Soweto      
#' 3 0.334    Nairobi

# combined
summary(lmer(inv_simpson~hiv_status + (1|general_site), 
             data=df.alpha %>% bind_rows(df.alpha %>% 
                                           mutate(facet='sites')) %>% 
               filter(hiv_detail!='HIV Positive(ART-)')))
#' Fixed effects:
#' Estimate Std. Error       df t value Pr(>|t|)    
#' (Intercept)   35.532      2.026    2.222  17.525  0.00199 ** 
#' hiv_status.L  -4.773      1.027 1692.856  -4.675 3.17e-06 ***

# ##############################################################################
# Beta

tmp <- .f_beta(feat.red, dist='bray', log.n0=1e-04)

df.plot <- tmp$coords %>% 
  left_join(meta.hiv %>% rename(Sample=general_sample_id)) %>% 
  filter(hiv_detail!='HIV Positive(ART-)')

g0 <- df.plot %>%
  ggplot(aes(x=PCo1, y=PCo2, col=general_site, alpha=hiv_detail)) +
    geom_point(pch=16) +
    theme_bw() + theme(panel.grid = element_blank()) + 
    scale_colour_manual(values=site.colours, guide='none') + 
    scale_alpha_manual(values=c(0.33, 0.85), guide='none')
g1 <- df.plot %>% 
  ggplot(aes(x=hiv_detail, y=PCo1, fill=hiv_detail)) + 
    geom_boxplot() + 
    theme_bw() + 
    theme(panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank()) + 
    coord_flip() + 
    scale_fill_manual(values=hiv_med, guide='none') 
g2 <- df.plot %>% 
  ggplot(aes(x=hiv_detail, y=PCo2, fill=hiv_detail)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_manual(values=hiv_med, guide='none') 
g <- cowplot::plot_grid(g1, NULL, g0, g2, ncol = 2, rel_widths = c(0.8, 0.2),
                        rel_heights = c(0.2, 0.8))
g0$data %>% 
  select(Sample, PCo1, PCo2, site.city, hiv_detail) %>% 
  write_csv('./files/source_data/Fig4c.csv')
ggsave(g, filename=here('figures', 'hiv', 'hiv_beta_newdb.pdf'),
       width = 6, height = 5, useDingbats=FALSE)

# test
df.plot %>% 
  select(Sample, PCo1, PCo2, site) %>% 
  pivot_longer(c(PCo1, PCo2)) %>% 
  group_by(name) %>% 
  group_map(.f=function(.x, .y){
    t <- kruskal.test(data=.x, value~site)
    tibble(p.val=t$p.value, .y, type='site')
  }) %>% 
  bind_rows()

#' # A tibble: 2 × 3
#' p.val name  type 
#' <dbl> <chr> <chr>
#' 1 1.43e- 9 PCo1  site 
#' 2 2.30e-28 PCo2  site 
df.plot %>% 
  select(Sample, PCo1, PCo2, hiv_detail) %>% 
  pivot_longer(c(PCo1, PCo2)) %>% 
  group_by(name) %>% 
  group_map(.f=function(.x, .y){
    t <- kruskal.test(data=.x, value~hiv_detail)
    tibble(p.val=t$p.value, .y, type='hiv_detail')
  }) %>% 
  bind_rows()

#' # A tibble: 2 × 3
#' p.val name  type      
#' <dbl> <chr> <chr>     
#' 1 1.03e- 1 PCo1  hiv_detail
#' 2 6.22e-12 PCo2  hiv_detail

# ##############################################################################
# testing
feat.rel <- prop.table(feat.red, 2)
feat.rel <- feat.rel[rowMeans(feat.rel > 1e-04) > 0.05,]
meta.red <- meta.hiv %>% 
  filter(hiv_detail!='HIV Positive(ART-)') %>% 
  select(general_sample_id, general_site, microbiome_antibiotics, 
         microbiome_diarrhea_last, microbiome_deworming_treatment, 
         microbiome_probiotics, household_size, demographic_age, 
         ses_employment, anthropometric_bmi, card_hypertension_status, 
         card_high_cholesterol, hiv_status) %>% 
  as.data.frame()
rownames(meta.red) <- meta.red$general_sample_id
meta.red$general_sample_id <- NULL


pb <- progress::progress_bar$new(total=nrow(feat.rel))
df.volcano <- map(rownames(feat.rel), .f = function(x){
  tmp <- enframe(feat.rel[x,], name='general_sample_id') %>% 
    full_join(meta.hiv, by='general_sample_id') %>% 
    mutate(value=log10(value + 5e-05)) %>% 
    filter(hiv_detail!='HIV Positive(ART-)')
  fit <- suppressMessages(lmer(value~hiv_detail + 
                                 (1|general_site) + 
                                 (1|microbiome_antibiotics) + 
                                 (1|microbiome_diarrhea_last), 
                               data=tmp))
  res <- coefficients(summary(fit))
  t <- tibble(species=x, p.value=res[2,5], effect.size=res[2,1])
  sites <- list()
  for (s in unique(tmp$general_site)){
    tmp2 <- tmp %>% filter(general_site==s)
    gfc <- tmp2 %>% group_by(hiv_detail) %>% 
      group_map(.f=function(.x, .y){
        enframe(quantile(.x$value, probs=seq(from=0.05, to=0.95, by=0.05))) %>% 
          mutate(.y)}) %>% 
      bind_rows() %>% 
      pivot_wider(names_from = hiv_detail, values_from = value) %>% 
      mutate(diff=`HIV Positive(ART+)`-`HIV Negative(Control)`) %>% 
      pull(diff) %>% mean
    
    
    t[[s]] <- gfc
  }
  tmp$hiv_detail <- as.factor(tmp$hiv_detail)
  
  res <- map(setdiff(colnames(tmp), 
                     c('general_sample_id', 'value', 'meta_site_comparison', 
                       'meta_genome_catalogue', 'meta_hiv_comparison', 
                       'meta_long_term_follow_up', 'meta_ltfu_reason', 
                       'general_study_id', 'general_tube_id', 
                       'general_enrollment_date', 
                       'anthropometric_waist_circumference',
                       'anthropometric_hip_circumference', 
                       'card_cholesterol_treatment',
                       'card_diabetes_treatment', 'card_esr_crp', 
                       'card_hypertension_treatment',
                       'card_rheumatoid_factor', 'ultrasound_vat', 
                       'ultrasound_scat',
                       'general_country', 'general_region', 
                       'demographic_gender',
                       'hiv_status', 'hiv_medication', 
                       'lab_fasting_confirmed',
                       'card_diabetes_treatment_spec', 
                       'household_pottable_water', 
                       'card_esr_crp')),
             .f=function(meta){
               tmp[[meta]] <- as.numeric(tmp[[meta]])
               x2 <- anova(lm(formula=paste0('value~', meta), data=tmp))
               tibble(meta=meta, var=x2$`Sum Sq`[1]/sum(x2$`Sum Sq`))
             }) %>% bind_rows()
  res <- res %>% 
    mutate(species=x, p.value=t$p.value, effect.size=t$effect.size, 
           Agincourt=t$Agincourt, Nairobi=t$Nairobi, Soweto=t$Soweto)
  pb$tick()
  return(res)}) %>% bind_rows()

df.volcano <- df.volcano %>% 
  pivot_wider(names_from = meta, values_from = var) %>% 
  mutate(p.adj=p.adjust(p.value, method='fdr')) %>% 
  arrange(p.adj) %>% 
  pivot_longer(-c(species, p.value, p.adj, effect.size, 
                  Agincourt, Nairobi, Soweto)) 

# # export table as supplement
# df.volcano %>% select(species, p.value, p.adj, effect.size,
#                       Agincourt, Nairobi, Soweto) %>%
#   distinct() %>%
#   filter(p.adj < 0.01) %>% write_tsv(here('files', 'hiv_table_newdb.tsv'))


g1 <- df.volcano %>% 
  ungroup() %>% 
  select(species, p.adj, effect.size) %>% 
  distinct() %>% 
  filter(p.adj < 0.01) %>% 
  mutate(order=sign(effect.size) *log10(p.adj)) %>% 
  arrange(order) %>% 
  mutate(species=factor(species, levels=species)) %>% 
  mutate(type=p.adj < 1e-05) %>% 
  ggplot(aes(x=species, y=-log10(p.adj), fill=type)) + 
  geom_col() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + 
  scale_fill_manual(values=c('#DAD7CB', '#7F7776'), guide='none')

g2 <- df.volcano %>% 
  ungroup() %>% 
  select(species, p.adj, effect.size, Agincourt, Soweto, Nairobi) %>% 
  distinct() %>% 
  filter(p.adj < 0.01) %>% 
  mutate(order=sign(effect.size) *log10(p.adj)) %>% 
  arrange(order) %>% 
  mutate(species=factor(species, levels=species)) %>% 
  select(species, Agincourt, Soweto, Nairobi) %>% 
  pivot_longer(-species) %>% 
  mutate(name=factor(name, levels=c('Agincourt','Soweto', 'Nairobi'))) %>% 
  ggplot(aes(x=species, y=name, fill=-value)) + 
  geom_tile() + 
  ggthemes::scale_fill_gradient2_tableau(palette='Red-Blue-White Diverging',
                                         name='Effect size') + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank()) + 
  ylab('')

g3 <- df.volcano %>% 
  ungroup() %>% 
  select(species, p.adj, effect.size, Agincourt, Soweto, Nairobi) %>% 
  distinct() %>% 
  filter(p.adj < 0.01) %>% 
  mutate(order=sign(effect.size) *log10(p.adj)) %>% 
  arrange(order) %>% 
  mutate(novel=str_detect(species, 'NEWDB')) %>% 
  mutate(species=factor(species, levels=species)) %>% 
  ggplot(aes(x=species, y=1, fill=novel)) + 
    geom_tile() + 
    scale_fill_manual(values=c('grey', 'black')) + 
    theme_bw() + theme(panel.grid = element_blank(),
                       axis.ticks = element_blank(), 
                       axis.text = element_blank())

df.volcano %>% 
  filter(p.adj < 1e-05) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  mutate(order=sign(effect.size) *log10(p.adj)) %>% 
  arrange(order) %>% 
  mutate(species=factor(species, levels=species)) %>% 
  select(-c(Agincourt, Nairobi, Soweto, p.value))%>% 
  pivot_longer(-c(species, p.adj, effect.size, hiv_detail)) %>% 
  mutate(value=value > hiv_detail) %>% 
  ggplot(aes(x=species, y=name, fill=value)) + 
  geom_tile()

# all the household variables are a bit odd... its because it has data for
# Soweto, but is non-applicable for AG/SW --> so it catches a lot of variation

# how about an actual volcano plot
g <- df.volcano %>% 
  select(species, p.value, effect.size, p.adj) %>% 
  distinct() %>% 
  mutate(type=case_when(p.adj < 1e-05 ~ '007C92', 
                        p.adj < 0.01 ~ 'medium',
                        TRUE~'other')) %>% 
  ggplot(aes(x=effect.size, y=-log10(p.adj), colour=type)) + 
    geom_point() +
    theme_bw() + theme(panel.grid.minor = element_blank()) + 
    xlab('Linear model effect size') + ylab('-log10(q-value)') +
  scale_colour_manual(values=c('#007C92','#AABEC6','#7F7776'))
ggsave(g, filename='./figures/hiv/bacteria_volcano.pdf',
       width = 4, height = 3, useDingbats=FALSE)


# add a phylum strip to the heatmap?
gtdb.tax.all <- read_tsv('./files/mOTUs_ext_GTDB_tax.tsv', 
  col_names = c('motus', 'domain', 'phylum', 'class', 
                'order', 'family', 'genus', 'species'),
  col_types = cols()) 

df.volcano %>% 
  ungroup() %>% 
  select(species, p.adj, effect.size, Agincourt, Soweto, Nairobi) %>% 
  distinct() %>% 
  filter(p.adj < 0.01) %>% 
  mutate(motus=str_extract(species, '\\[(ref|ext|meta|NEW).*\\]')) %>% 
  mutate(motus=str_remove(motus, '\\[')) %>% 
  mutate(motus=str_remove(motus, '\\]')) %>% 
  left_join(gtdb.tax.all %>% 
              transmute(motus, species_gtdb=species), by='motus') %>% 
  select(-motus) %>% 
  write_csv('./files/hiv_species_associations.csv')

g <- cowplot::plot_grid(g1, g2, g3, ncol=1, rel_heights = c(0.6, .2, 0.2))
ggsave(g, filename=here('figures', 'hiv', 'hiv_effect_size_heatmap_newdb.pdf'),
       width = 12, height = 6)
# export source data
g1$data %>% 
  write_csv('./files/source_data/Fig4d.csv')
g2$data %>% 
  write_csv('./files/source_data/Fig4d2.csv')
g3$data %>% 
  write_csv('./files/source_data/Fig4d3.csv')

df.volcano %>% 
  select(species, p.value, p.adj, effect.size) %>% 
  distinct() %>% 
  filter(p.adj < 0.01) %>% 
  dim
# 131
df.volcano %>% 
  select(species, p.value, p.adj, effect.size) %>% 
  distinct() %>% 
  filter(p.adj < 1e-05) %>% 
  arrange(effect.size) %>% 
  mutate(motus=str_extract(species, '\\[(NEWDB|ref|meta|ext).*\\]')) %>% 
  mutate(motus=str_remove(motus, '\\[')) %>% 
  mutate(motus=str_remove(motus, '\\]')) %>% 
  left_join(gtdb.tax.all %>% select(motus, species), by='motus') %>% 
  select(species.y, p.value, p.adj, effect.size) 

#' # A tibble: 13 × 4
#' species                            p.value    p.adj effect.size
#' <chr>                                <dbl>    <dbl>       <dbl>
#' 1 s__Faecalibacterium prausnitzii_E 2.57e- 8 3.24e- 6      -0.453
#' 2 s__Limivicinus sp900547315        2.22e-10 6.07e- 8      -0.411
#' 3 s__CAG-170 sp003516765            1.65e-12 1.35e- 9      -0.400
#' 4 s__Ventricola sp900542395         2.74e-10 6.40e- 8      -0.400
#' 5 s__Faecalibacterium prausnitzii_I 7.81e-10 1.60e- 7      -0.389
#' 6 s__Faecalibacterium prausnitzii_D 8.99e- 9 1.27e- 6      -0.370
#' 7 s__Limivicinus sp003150355        9.33e- 9 1.27e- 6      -0.355
#' 8 s__Faecousia sp900549705          3.19e-11 1.04e- 8      -0.343
#' 9 s__Faecalibacterium sp900539885   6.96e- 8 8.14e- 6      -0.304

#' 10 s__Pyramidobacter piscolens       2.73e- 9 4.97e- 7       0.153
#' 11 s__Enterocloster sp000431375      7.57e- 9 1.24e- 6       0.327
#' 12 s__Fusobacterium_A mortiferum     1.59e-11 6.50e- 9       0.330
#' 13 s__Dysosmobacter welbionis        4.18e-12 2.28e- 9       0.350
#' 14 s__Megamonas funiformis           6.04e-19 9.89e-16       0.565

# ##############################################################################
# SIAMCAT?
# ML overall & ML per site
set.seed(1791)

library("SIAMCAT")
# feat.rel <- prop.table(motus.gtdb.lvls$species, 2)[,colnames(feat.rel)]
if (!file.exists(here('files','hiv_models_newdb.RData'))){
  
  sc.obj <- siamcat(feat=feat.rel, meta=meta.red, label='hiv_status', case='Yes')
  sc.obj <- filter.features(sc.obj, filter.method='pass', cutoff=0.05)
  # check.confounders(sc.obj, fn.plot = '~/Desktop/confounder.pdf', verbose=3)
  sc.obj <- normalize.features(sc.obj, norm.method='log.std', 
                               norm.param=list(log.n0=5e-05, sd.min.q=0))
  sc.obj <- create.data.split(sc.obj, num.folds = 5, num.resample = 5)
  sc.obj <- train.model(sc.obj, method='lasso', grid.size = 5)
  sc.obj <- make.predictions(sc.obj)
  sc.obj <- evaluate.predictions(sc.obj)
  # model.evaluation.plot(sc.obj)
  
  sc.list <- list()
  sc.list[['all']] <- sc.obj
  
  for (site in c('Agincourt', 'Nairobi', 'Soweto')){
    tmp <- siamcat(feat=feat.rel, meta=meta.red[meta.red$general_site==site,], 
                   label='hiv_status', case='Yes', verbose=0)
    tmp <- filter.features(tmp, filter.method='prevalence', cutoff=0.05,
                           verbose=0)
    tmp <- normalize.features(tmp, norm.method='log.std', 
                                 norm.param=list(log.n0=5e-05, sd.min.q=0),
                              verbose=0)
    tmp <- create.data.split(tmp, num.folds = 5, num.resample = 5, verbose=0)
    tmp <- train.model(tmp, method='lasso', verbose=0)
    tmp <- make.predictions(tmp, verbose=0)
    tmp <- evaluate.predictions(tmp, verbose=0)
    sc.list[[site]] <- tmp
    message(site, ': ', sprintf(fmt='%.2f', eval_data(tmp)$auroc))
    
    tmp <- siamcat(feat=feat.rel, meta=meta.red[meta.red$general_site!=site,], 
                   label='hiv_status', case='Yes', verbose=0)
    tmp <- filter.features(tmp, filter.method='prevalence', cutoff=0.05,
                           verbose=0)
    tmp <- normalize.features(tmp, norm.method='log.std', 
                              norm.param=list(log.n0=5e-05, sd.min.q=0),
                              verbose=0)
    tmp <- create.data.split(tmp, num.folds = 5, num.resample = 5, verbose=0)
    tmp <- train.model(tmp, method='lasso', verbose=0)
    tmp <- make.predictions(tmp, verbose=0)
    tmp <- evaluate.predictions(tmp, verbose=0)
    sc.list[[paste0(site, 'left_out')]] <- tmp
  }
  save(sc.list, file=here('files', 'hiv_models_newdb.RData'))
} else {
  load(here('files', 'hiv_models_newdb.RData'))
}
df.auc <- tibble(train=character(0), test=character(0), auc=double(0),
                 type=character(0))
df.auc <- df.auc %>% 
  add_row(train='all', test='all', 
          auc=as.numeric(sc.list$all@eval_data$auroc), type='CV')

for (site in c('Agincourt', 'Nairobi', 'Soweto')){
  tmp <- sc.list[[site]]
  for (site2 in c('Agincourt', 'Nairobi', 'Soweto')){
    if (site==site2){
      df.auc <- df.auc %>% 
        add_row(train=site, test=site, 
                auc=as.numeric(tmp@eval_data$auroc), type='CV')
    } else {
      tmp.test <- siamcat(feat=feat.rel, 
                          meta=meta.red[meta.red$general_site==site2,], 
                          label='hiv_status', case='Yes', verbose=0)
      tmp.test <- make.predictions(tmp, tmp.test, verbose=0)
      tmp.test <- evaluate.predictions(tmp.test, verbose=0)
      df.auc <- df.auc %>% 
        add_row(train=site, test=site2, 
                auc=as.numeric(tmp.test@eval_data$auroc), type='test')
    }
  }
  tmp <- sc.list[[paste0(site, 'left_out')]]
  tmp.test <- siamcat(feat=feat.rel, 
                      meta=meta.red[meta.red$general_site==site,], 
                      label='hiv_status', case='Yes', verbose=0)
  tmp.test <- make.predictions(tmp, tmp.test, verbose=0)
  tmp.test <- evaluate.predictions(tmp.test, verbose=0)
  df.auc <- df.auc %>% 
    add_row(train='left_out', test=site, 
            auc=as.numeric(tmp.test@eval_data$auroc), type='test')
}

model.evaluation.plot('all'=sc.list$all, 'AG'=sc.list$Agincourt, 
                      'KY'=sc.list$Nairobi, 'SW'=sc.list$Soweto,
                      colours = c('black', '#889A4D', '#4A4C59', '#82354E'),
                      fn.plot = here('figures', 'hiv', 'ml_auroc_newdb.pdf'))
tibble(model='all',
       sensitivity=sc.list$all@eval_data$roc$sensitivities,
       specificity=sc.list$all@eval_data$roc$specificities) %>% 
  bind_rows(tibble(
    model='Agincourt',
    sensitivity=sc.list$Agincourt@eval_data$roc$sensitivities,
    specificity=sc.list$Agincourt@eval_data$roc$specificities)) %>% 
  bind_rows(tibble(
    model='Nairobi',
    sensitivity=sc.list$Nairobi@eval_data$roc$sensitivities,
    specificity=sc.list$Nairobi@eval_data$roc$specificities)) %>% 
  bind_rows(tibble(
    model='Soweto',
    sensitivity=sc.list$Soweto@eval_data$roc$sensitivities,
    specificity=sc.list$Soweto@eval_data$roc$specificities)) %>% 
  write_csv('./files/source_data/Fig4e.csv')

g <- df.auc %>% 
  mutate(train=factor(train, levels=rev(c('all', 'Agincourt', 'Soweto', 
                                          'Nairobi', 'left_out')))) %>% 
  mutate(test=factor(test, levels=c('all', 'Agincourt', 'Soweto', 
                                    'Nairobi'))) %>% 
  mutate(label=sprintf(fmt='%.2f', auc)) %>% 
  ggplot(aes(y=train, x=test, fill=auc)) + 
    geom_tile() + 
    scale_fill_viridis_c(limits=c(0.5, 1), name='AUROC') + 
    geom_text(aes(label=label)) + 
    theme_bw() + theme(panel.grid = element_blank()) + 
    ylab('Training set') + xlab('Test set')
ggsave(g, filename=here('figures', 'hiv', 'ml_model_transfer_newdb.pdf'),
       width = 6, height = 5, useDingbats=FALSE)
g$data %>% 
  write_csv('./files/source_data/Fig4f.csv')

# apply to other sites?
test.ext <- prop.table(motus.gtdb.lvls$motus, 2)
sc.test <- siamcat(feat = test.ext, verbose=0)
df.cross <- list()
for (model.name in names(sc.list)){
  model <- sc.list[[model.name]]
  sc.test <- make.predictions(model, sc.test)
  cutoff.idx <- which(model@eval_data$roc$specificities >= 0.95)[1]
  cutoff <- model@eval_data$roc$thresholds[cutoff.idx]
  df.cross[[model.name]] <- enframe(rowMeans(pred_matrix(sc.test)), 
                                    name = 'general_sample_id') %>% 
    mutate(positive=value >= cutoff) %>% 
    right_join(metadata.clean %>% filter(meta_site_comparison)) %>% 
    filter(general_site %in% c('Nanoro', 'Navrongo', 'Dimamo')) %>% 
    mutate(hiv_status=case_when(is.na(hiv_status)~'No', TRUE~hiv_status)) %>% 
    filter(!is.na(value)) %>%
    group_by(general_site, hiv_status) %>% 
    summarise(m=mean(positive)) %>% 
    mutate(model=model.name)
}
g <- df.cross$all %>% 
  ggplot(aes(x=hiv_status, y=m, fill=hiv_status)) + 
    geom_col() +
    facet_grid(~general_site, space = 'free', scales = 'free') + 
    geom_hline(yintercept = 0.05) + 
    scale_fill_manual(values=hiv_med) + 
    theme_bw() + xlab('') + ylab('FPR/TPR') + 
    theme(panel.grid = element_blank(), axis.ticks.x = element_blank(),
          axis.text.x = element_blank(), strip.background = element_blank()) 
g$data %>% 
  write_csv('./files/source_data/Fig4g.csv')
ggsave(g, filename=here('figures', 'hiv', 'hiv_fpr_newdb.pdf'),
       width = 4, height = 6)

df.cross$all
#' # A tibble: 4 × 4
#' # Groups:   general_site [3]
#' general_site hiv_status      m model
#' <fct>        <chr>       <dbl> <chr>
#' 1 Nanoro       No         0.0314 all  
#' 2 Navrongo     No         0.0367 all  
#' 3 Dimamo       No         0.0256 all  
#' 4 Dimamo       Yes        0.667  all 

# model introspection?
# relate this to the univariate testing?
df.feat.weight <- map(names(sc.list), .f = function(x){
  as_tibble(feature_weights(sc.list[[x]]), rownames='species') %>% 
    select(species, mean.rel.weight) %>% 
    mutate(name=x)
}) %>% bind_rows()

df.volcano %>% select(species, p.adj, effect.size) %>% 
  distinct() %>% full_join(df.feat.weight, by='species') %>% 
  filter(p.adj < 0.01) %>%
  mutate(order=sign(effect.size) *log10(p.adj)) %>% 
  arrange(order) %>% 
  mutate(species=factor(species, levels=unique(species))) %>% 
  filter(name %in% c('all', 'Agincourt', 'Nairobi', 'Soweto')) %>% 
  ggplot(aes(x=species, y=name, fill=mean.rel.weight)) + 
    geom_tile() + 
    ggthemes::scale_fill_gradient2_tableau(palette='Green-Blue-White Diverging',
                                           name='Effect size') +
    theme_bw() + 
    theme(axis.text = element_blank(), axis.ticks = element_blank()) 
  

df.weight.red <- df.feat.weight %>% 
  group_by(name) %>% 
  pivot_wider(names_from = name, values_from = mean.rel.weight) %>% 
  arrange(desc(abs(all))) %>% 
  full_join(df.volcano %>% select(species, p.value, effect.size, p.adj) %>% 
              distinct()) %>% 
  mutate(Agincourt=rank(-abs(Agincourt)),
         Agincourtleft_out=rank(-abs(Agincourtleft_out)),
         Nairobi=rank(-abs(Nairobi)),
         Nairobileft_out=rank(-abs(Nairobileft_out)),
         Soweto=rank(-abs(Soweto)),
         Sowetoleft_out=rank(-abs(Sowetoleft_out)),
         all.rank=rank(-abs(all))) %>% 
  head(n=50)
g <- df.weight.red %>% 
  arrange(all) %>% 
  mutate(species=factor(species, levels=species)) %>% 
  pivot_longer(c(Nairobi, Nairobileft_out, Agincourt, Agincourtleft_out, 
                 Soweto, Sowetoleft_out, all.rank)) %>% 
  filter(value < 50) %>% 
  mutate(name=factor(name, levels=c('all.rank', 'Agincourtleft_out', 
                                    'Nairobileft_out', 'Sowetoleft_out', 
                                    'Agincourt', 'Nairobi', 'Soweto'))) %>% 
  ggplot(aes(x=name, y=species, fill=value)) + 
    geom_tile()
  

g2.1 <- df.weight.red %>% 
  arrange(all) %>% 
  mutate(species=factor(species, levels=species)) %>% 
  ggplot(aes(x=species, y=-log10(p.adj))) +
    geom_col() + coord_flip() + 
    theme_bw() + theme(axis.text.y=element_blank())

g2 <- df.weight.red %>% 
  arrange(all) %>% 
  mutate(species=factor(species, levels=species)) %>% 
  ggplot(aes(x=species, y=1, fill=effect.size)) +
    geom_tile() + coord_flip() + 
    theme_bw() + theme(axis.text.y=element_blank()) +
    ggthemes::scale_fill_gradient2_tableau(palette='Red-Blue-White Diverging', 
                                           guide='none')


cowplot::plot_grid(g, g2, g2.1, nrow=1, rel_widths = c(0.6, 0.2, 0.2))

# ##############################################################################
# phage differences?
# phanta.rel <- read.table(here('data', 'classification', 'phanta', 
#                               'relative_taxonomic_abundance.txt'),
#                          sep='\t', header = TRUE)
# rownames(phanta.rel) <- phanta.rel$Taxon_Lineage_with_Names
# phanta.red <- phanta.rel[,meta.hiv$general_sample_id]
# 
# tmp <- which(str_detect(rownames(phanta.red), 'species_'))
# phanta.red <- phanta.red[tmp,]
# tmp <- which(str_detect(rownames(phanta.red), 'superkingdom_Viruses'))
# phanta.red <- phanta.red[tmp,]
# phanta.red <- prop.table(as.matrix(phanta.red), 2)
# 
# hist(rowMeans(phanta.red > 1e-05), 100)
# 
# phanta.filt <- phanta.red[rowMeans(phanta.red>1e-05) > 0.1,]

load(here('data', 'classification', 'all_classification_tables.RData'))
tbl.phanta <- lst.phanta.vir$species
df.alpha.phage <- enframe(colSums(tbl.phanta > 1e-04), 
                          value='richness', 
                    name='general_sample_id') %>% 
  right_join(meta.hiv, by='general_sample_id')
phanta.filt <- tbl.phanta[rowMeans(
  tbl.phanta[,meta.hiv$general_sample_id]>1e-05) > 0.1,
  meta.hiv$general_sample_id]
# alpha
seed <- 1933
# alpha <- .f_alpha(phanta.filt, rarefy = FALSE)
# df.alpha.phage <- alpha$alpha_df %>% 
  # rename(general_sample_id=Sample_ID) %>% 
  # full_join(meta.hiv, by='general_sample_id')

g <- df.alpha.phage %>% 
  mutate(general_site='all') %>% 
  mutate(facet='all') %>% 
  bind_rows(df.alpha.phage %>% 
              mutate(facet='sites')) %>% 
  filter(hiv_detail!='HIV Positive(ART-)') %>% 
  mutate(general_site=factor(general_site, 
                             levels=c('all', 'Agincourt', 
                                      'Soweto', 'Nairobi'))) %>% 
  ggplot(aes(x=general_site, y=richness, fill=hiv_detail)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), 
              colour='#D3D3D3') + 
  geom_boxplot(alpha=0.6, outlier.shape = NA, col='black') + 
  facet_grid(~facet, scales = 'free', space='free') + 
  scale_fill_manual(values=hiv_med, name='') +
  xlab('') + ylab("Alpha diversity") +
  theme_bw() + theme(panel.grid.minor = element_blank(), 
                     panel.grid.major.x=element_blank(),
                     strip.background = element_blank(), 
                     strip.text = element_blank()) +
  NULL

ggsave(g, filename=here('figures','hiv','hiv_alpha_phanta.pdf'), 
       width = 10, height = 5, useDingbats=FALSE)

df.alpha.phage %>% 
  group_by(general_site) %>% 
  bind_rows(df.alpha.phage %>% 
              mutate(facet='sites')) %>% 
  filter(hiv_detail!='HIV Positive(ART-)') %>% 
  group_map(.f=function(.x, .y){
    res <- coefficients(summary(lm(data=.x, richness~hiv_status)))
    tibble(p.val=res[2,4], .y)
  }) %>% bind_rows()

#'
#' # A tibble: 3 × 2 
#' p.val general_site
#' <dbl> <fct>       
#' 1 0.0757    Agincourt    
#' 2 0.0000203 Soweto    
#' 3 0.844     Nairobi    

# combined
summary(lmer(richness~hiv_status + (1|general_site), 
             data=df.alpha.phage %>% bind_rows(df.alpha.phage %>% 
                                           mutate(facet='sites')) %>% 
               filter(hiv_detail!='HIV Positive(ART-)')))
#'              Estimate Std. Error       df t value Pr(>|t|)    
#'  (Intercept)   389.180     19.154    2.070  20.318 0.002049 ** 
#'  hiv_status.L  -17.675      4.684 1693.273  -3.774 0.000166 ***

# ##############################################################################
# Beta

tmp <- .f_beta(phanta.filt, dist='bray', log.n0=1e-05)

df.plot <- tmp$coords %>% 
  left_join(meta.hiv %>% rename(Sample=general_sample_id)) %>% 
  filter(hiv_detail!='HIV Positive(ART-)')

g0 <- df.plot %>%
  ggplot(aes(x=PCo1, y=PCo2, col=general_site, alpha=hiv_detail)) +
  geom_point(pch=16) +
  theme_bw() + theme(panel.grid = element_blank()) + 
  scale_colour_manual(values=site.colours, guide='none') + 
  scale_alpha_manual(values=c(0.33, 0.85), guide='none')
g1 <- df.plot %>% 
  ggplot(aes(x=hiv_detail, y=PCo1, fill=hiv_detail)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  coord_flip() + 
  scale_fill_manual(values=hiv_med, guide='none') 
g2 <- df.plot %>% 
  ggplot(aes(x=hiv_detail, y=PCo2, fill=hiv_detail)) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank()) + 
  scale_fill_manual(values=hiv_med, guide='none') 
g <- cowplot::plot_grid(g1, NULL, g0, g2, ncol = 2, rel_widths = c(0.8, 0.2),
                        rel_heights = c(0.2, 0.8))
ggsave(g, filename=here('figures', 'hiv', 'hiv_beta_phanta.pdf'),
       width = 6, height = 5, useDingbats=FALSE)

# species associations?
df.volcano.phage <- map(rownames(phanta.filt), .f = function(x){
  df.x <- enframe(phanta.filt[x,], name='general_sample_id') %>% 
    mutate(value=log10(value+1e-05)) %>% 
    left_join(meta.hiv, by='general_sample_id') %>% 
    filter(hiv_detail!='HIV Positive(ART-)')
  fit <- lmerTest::lmer(value~hiv_detail + (1|general_site), data=df.x)
  coefs <- coefficients(summary(fit))
  tibble(species=x, coef=coefs[2,1], pval=coefs[2,5])
}) %>% bind_rows() %>% 
  mutate(qval=p.adjust(pval, method='fdr')) 

df.volcano.phage %>% 
  filter(qval < 0.01) %>% dim

df.volcano.phage %>%
  filter(qval < 0.01) %>%
  arrange(qval) %>%
  write_csv('./files/hiv_phage_associations.csv')

g <- df.volcano.phage %>% 
  mutate(type=case_when(qval < 1e-05 ~ '007C92', 
                        qval < 0.01 ~ 'medium',
                        TRUE~'other')) %>% 
  ggplot(aes(x=coef, y=-log10(qval), colour=type)) + 
  geom_point() +
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  xlab('Linear model effect size') + ylab('-log10(q-value)') +
  scale_colour_manual(values=c('#007C92','#AABEC6','#7F7776'))
ggsave(g, filename='./figures/hiv/phanta_volcano.pdf',
       width = 4, height = 3, useDingbats=FALSE)


# 
x <- df.volcano.phage %>% arrange(pval) %>% slice_head(n=1) %>% pull(species)
df.x <- enframe(phanta.filt[x,], name='general_sample_id') %>% 
  mutate(value=log10(value+1e-05)) %>% 
  left_join(meta.hiv, by='general_sample_id') %>% 
  filter(hiv_detail!='HIV Positive(ART-)')
df.x %>% 
  ggplot(aes(x=general_site, y=value, fill=hiv_detail)) + 
    geom_boxplot() + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) 


g.ef <- 
  df.volcano.phage %>% 
  mutate(p.adj=p.adjust(pval, method='fdr')) %>% 
  mutate(data='phage') %>% 
  bind_rows(df.volcano %>% transmute(species, p.adj, coef=effect.size) %>% 
              distinct() %>% 
              mutate(data='motus')) %>% 
  mutate(type=case_when(p.adj < 1e-05 ~ '007C92', 
                        p.adj < 0.01 ~ 'medium',
                        TRUE~'other')) %>% 
  ggplot(aes(x=coef, y=-log10(p.adj), colour=type)) + 
  geom_point() +
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  xlab('Linear model effect size') + ylab('-log10(q-value)') +
  scale_colour_manual(values=c('#007C92','#AABEC6','#7F7776')) +
  facet_grid(~data) +
    geom_hline(yintercept = c(-log10(0.01), -log10(0.00001)))
ggsave(g.ef, filename='./figures/hiv/hiv_motus_phanta_volcano.pdf',
       width = 6, height = 4, useDingbats=FALSE)  
g.ef$data %>% 
  select(species, coef, p.adj, type, data) %>% 
  write_csv('./files/source_data/ED10d.csv')

# ##############################################################################
# phage models as well
library("SIAMCAT")

if (!file.exists(here('files','hiv_models_phanta.RData'))){
  meta.red <- meta.hiv %>% 
    filter(hiv_detail!='HIV Positive(ART-)') %>% 
    select(general_sample_id, general_site, microbiome_antibiotics, 
           microbiome_diarrhea_last, microbiome_deworming_treatment, 
           microbiome_probiotics, household_size, demographic_age, 
           ses_employment, anthropometric_bmi, card_hypertension_status, 
           card_high_cholesterol, hiv_status) %>% 
    as.data.frame()
  rownames(meta.red) <- meta.red$general_sample_id
  meta.red$general_sample_id <- NULL
  
  sc.obj <- siamcat(feat=phanta.filt, meta=meta.red, label='hiv_status', case='Yes')
  sc.obj <- filter.features(sc.obj, filter.method='pass', cutoff=0.05)
  # check.confounders(sc.obj, fn.plot = '~/Desktop/confounder.pdf', verbose=3)
  sc.obj <- normalize.features(sc.obj, norm.method='log.std', 
                               norm.param=list(log.n0=1e-05, sd.min.q=0))
  sc.obj <- create.data.split(sc.obj, num.folds = 5, num.resample = 5)
  sc.obj <- train.model(sc.obj, method='lasso', grid.size = 5)
  sc.obj <- make.predictions(sc.obj)
  sc.obj <- evaluate.predictions(sc.obj)
  # model.evaluation.plot(sc.obj)
  
  sc.list <- list()
  sc.list[['all']] <- sc.obj
  
  for (site in c('Agincourt', 'Nairobi', 'Soweto')){
    tmp <- siamcat(feat=phanta.filt, 
                   meta=meta.red[meta.red$general_site==site,], 
                   label='hiv_status', case='Yes', verbose=0)
    tmp <- filter.features(tmp, filter.method='prevalence', cutoff=0.05,
                           verbose=0)
    tmp <- normalize.features(tmp, norm.method='log.std', 
                              norm.param=list(log.n0=5e-05, sd.min.q=0),
                              verbose=0)
    tmp <- create.data.split(tmp, num.folds = 5, num.resample = 5, verbose=0)
    tmp <- train.model(tmp, method='lasso', verbose=0)
    tmp <- make.predictions(tmp, verbose=0)
    tmp <- evaluate.predictions(tmp, verbose=0)
    sc.list[[site]] <- tmp
    message(site, ': ', sprintf(fmt='%.2f', eval_data(tmp)$auroc))
    
    tmp <- siamcat(feat=phanta.filt, 
                   meta=meta.red[meta.red$general_site!=site,], 
                   label='hiv_status', case='Yes', verbose=0)
    tmp <- filter.features(tmp, filter.method='prevalence', cutoff=0.05,
                           verbose=0)
    tmp <- normalize.features(tmp, norm.method='log.std', 
                              norm.param=list(log.n0=5e-05, sd.min.q=0),
                              verbose=0)
    tmp <- create.data.split(tmp, num.folds = 5, num.resample = 5, verbose=0)
    tmp <- train.model(tmp, method='lasso', verbose=0)
    tmp <- make.predictions(tmp, verbose=0)
    tmp <- evaluate.predictions(tmp, verbose=0)
    sc.list[[paste0(site, 'left_out')]] <- tmp
  }
  save(sc.list, file=here('files', 'hiv_models_phanta.RData'))
} else {
  load(here('files', 'hiv_models_phanta.RData'))
}
df.auc <- tibble(train=character(0), test=character(0), auc=double(0),
                 type=character(0))
df.auc <- df.auc %>% 
  add_row(train='all', test='all', 
          auc=as.numeric(sc.list$all@eval_data$auroc), type='CV')

for (site in c('Agincourt', 'Nairobi', 'Soweto')){
  tmp <- sc.list[[site]]
  for (site2 in c('Agincourt', 'Nairobi', 'Soweto')){
    if (site==site2){
      df.auc <- df.auc %>% 
        add_row(train=site, test=site, 
                auc=as.numeric(tmp@eval_data$auroc), type='CV')
    } else {
      tmp.test <- siamcat(feat=phanta.filt, 
                          meta=meta.red[meta.red$general_site==site2,], 
                          label='hiv_status', case='Yes', verbose=0)
      tmp.test <- make.predictions(tmp, tmp.test, verbose=0)
      tmp.test <- evaluate.predictions(tmp.test, verbose=0)
      df.auc <- df.auc %>% 
        add_row(train=site, test=site2, 
                auc=as.numeric(tmp.test@eval_data$auroc), type='test')
    }
  }
  tmp <- sc.list[[paste0(site, 'left_out')]]
  tmp.test <- siamcat(feat=phanta.filt, 
                      meta=meta.red[meta.red$general_site==site,], 
                      label='hiv_status', case='Yes', verbose=0)
  tmp.test <- make.predictions(tmp, tmp.test, verbose=0)
  tmp.test <- evaluate.predictions(tmp.test, verbose=0)
  df.auc <- df.auc %>% 
    add_row(train='left_out', test=site, 
            auc=as.numeric(tmp.test@eval_data$auroc), type='test')
}

model.evaluation.plot('all'=sc.list$all, 'AG'=sc.list$Agincourt, 
                      'KY'=sc.list$Nairobi, 'SW'=sc.list$Soweto,
                      colours = c('black', '#889A4D', '#4A4C59', '#82354E'),
                      fn.plot = here('figures', 'hiv', 'ml_auroc_phanta.pdf'))
tibble(model='all',
       sensitivity=sc.list$all@eval_data$roc$sensitivities,
       specificity=sc.list$all@eval_data$roc$specificities) %>% 
  bind_rows(tibble(
    model='Agincourt',
    sensitivity=sc.list$Agincourt@eval_data$roc$sensitivities,
    specificity=sc.list$Agincourt@eval_data$roc$specificities)) %>% 
  bind_rows(tibble(
    model='Nairobi',
    sensitivity=sc.list$Nairobi@eval_data$roc$sensitivities,
    specificity=sc.list$Nairobi@eval_data$roc$specificities)) %>% 
  bind_rows(tibble(
    model='Soweto',
    sensitivity=sc.list$Soweto@eval_data$roc$sensitivities,
    specificity=sc.list$Soweto@eval_data$roc$specificities)) %>% 
  write_csv('./files/source_data/ED10e.csv')

g <- df.auc %>% 
  mutate(train=factor(train, levels=rev(c('all', 'Agincourt', 'Soweto', 
                                          'Nairobi', 'left_out')))) %>% 
  mutate(test=factor(test, levels=c('all', 'Agincourt', 'Soweto', 
                                    'Nairobi'))) %>% 
  mutate(label=sprintf(fmt='%.2f', auc)) %>% 
  ggplot(aes(y=train, x=test, fill=auc)) + 
  geom_tile() + 
  scale_fill_viridis_c(limits=c(0.5, 1), name='AUROC') + 
  geom_text(aes(label=label)) + 
  theme_bw() + theme(panel.grid = element_blank()) + 
  ylab('Training set') + xlab('Test set')
ggsave(g, filename=here('figures', 'hiv', 'ml_model_transfer_phanta.pdf'),
       width = 6, height = 5, useDingbats=FALSE)
# export source data
g$data %>% 
  write_csv('./files/source_data/ED10f.csv')

# ##############################################################################
# Other datasets?

.f_get_motus_level <- function(lvl='phylum', motus.tbl){
  stopifnot(lvl %in% colnames(motus.tax.gtdb)[2:8])
  rownames(motus.tbl) <- str_extract(
    rownames(motus.tbl), '((ref|meta|ext)_mOTU_v3_[0-9]{5}|(unassigned)|(NEWDB_[0-9]+))')
  new.level <- unique(motus.tax.gtdb[[lvl]])
  new.mat <- matrix(0, nrow = length(new.level)+1, ncol=ncol(motus.tbl),
                    dimnames = list(c(new.level, 'unassigned'), 
                                    colnames(motus.tbl)))
  for (x in new.level){
    incl.motus <- motus.tax.gtdb %>% 
      filter(eval(sym(lvl)) ==x) %>% 
      filter(motus %in% rownames(motus.tbl)) %>% 
      pull(motus)
    if (length(incl.motus) > 0){
      new.mat[x,] <- colSums(motus.tbl[incl.motus,,drop=FALSE])
    }
  }
  new.mat['unassigned',] <- motus.tbl['unassigned',]
  new.mat <- new.mat[rowMeans(new.mat!=0) > 0,]
  return(new.mat)
}
model <- sc.list$all
cutoff.idx <- which(model@eval_data$roc$specificities >= 0.95)[1]
cutoff <- model@eval_data$roc$thresholds[cutoff.idx]
df.other <- list()
for (i in c('Fulcher', 'NogueraJulian', 'Zhang_Frontiers', 'Bai', 'Rocafort')){
  message(i)
  f <- list.files('/Volumes/lab_asbhatt/data/public_data', pattern=i)
  motus.other <- read.table(
    paste0('/Volumes/lab_asbhatt/data/public_data/', 
           f, '/classification/motus_all_v3.0.3_AWIGEN.tsv'),
    sep='\t',stringsAsFactors = FALSE, check.names = FALSE,
    row.names = 1, header = TRUE, quote = '', 
    comment.char = '', skip=2)
  motus.other <- as.matrix(motus.other[rowMeans(motus.other!=0) > 0,])
  feat.other <- .f_get_motus_level(lvl='species', motus.other)
  
  x <- motus.gtdb.lvls$species
  u <- union(rownames(x), rownames(feat.other))
  mat.new <- matrix(0, nrow=length(u), ncol=ncol(x)+ncol(feat.other),
         dimnames = list(u, c(colnames(x), colnames(feat.other))))
  mat.new[rownames(x), colnames(x)] <- x
  mat.new[rownames(feat.other), colnames(feat.other)] <- feat.other
  df.pco <- labdsv::pco(vegan::vegdist(vegan::rrarefy(t(mat.new), 3000)))
  df.pco.plot <- as.data.frame(df.pco$points) %>% 
    as_tibble(rownames='sample_id') %>% 
    mutate(study=case_when(sample_id %in% colnames(x)~'AWIGEN', TRUE~i))
  print(df.pco.plot %>% 
    ggplot(aes(x=V1, y=V2, col=study)) + 
      geom_point())
  
  
  rest <- setdiff(rownames(get.orig_feat.matrix(model)), 
                  rownames(feat.other))
  
  print(feature_weights(model) %>% as_tibble(rownames='species') %>% 
    filter(species %in% rest) %>% 
    pull(mean.rel.weight) %>% 
    abs() %>% sum)
  
  feat.other <- rbind(prop.table(feat.other, 2), 
                      matrix(data=0, nrow=length(rest),
                             ncol=ncol(feat.other),
                             dimnames = list(rest, colnames(feat.other))))
  sc.test <- siamcat(feat = feat.other, verbose=0)
  sc.test <- make.predictions(model, sc.test)
  df.other[[i]] <- enframe(rowMeans(pred_matrix(sc.test)), 
                           name = 'sample_id') %>% 
    mutate(positive=value >= cutoff) 
}

df.other <- df.other %>% bind_rows(.id = 'study')

meta.others <- read_tsv(
  './files/hiv_other/meta_bai.tsv', col_types = cols()) %>% 
    select(sample_id, HIV_status) %>% 
    mutate(Country='Sweden') %>% 
    bind_rows(read_tsv('./files/hiv_other/meta_fulcher.tsv', 
                       col_types = cols()) %>% 
                select(sample_id, HIV_status, Country)) %>% 
    bind_rows(read_tsv('./files/hiv_other/meta_rocafort.tsv', 
                       col_types = cols()) %>% 
                mutate(Country='Mozambique', HIV_status='unclear') %>% 
                select(sample_id, HIV_status, Country)) %>% 
    bind_rows(read_tsv('./files/hiv_other/meta_noguera_julian.tsv', 
                       col_types = cols()) %>% 
                mutate(Country=paste0('Spain-', MSM, '-', Sex)) %>% 
                select(sample_id, HIV_status, Country)) %>% 
    bind_rows(read_tsv('./files/hiv_other/meta_zhang.tsv', 
                       col_types = cols()) %>% 
                mutate(Country='Netherlands') %>% 
                rename(HIV_status=HIV_Status) %>% 
                select(sample_id, HIV_status, Country)) %>% 
  mutate(HIV_status=case_when(HIV_status %in% c('Positive', 'positive', 'PLWH')~'PLWH',
                              HIV_status %in% c('Healthy', 'negative', 'Negative')~'HIV-negative',
                              TRUE~'unclear'))

g <- df.other %>% 
  left_join(meta.others, by='sample_id') %>% 
  mutate(study=case_when(study=='Zhang_Frontiers'~'Zhang', 
                         TRUE~study)) %>% 
  # bind_rows(enframe(rowMeans(model@pred_matrix), name='sample_id') %>% 
              # left_join(meta.hiv %>% transmute(sample_id=general_sample_id, HIV_status=hiv_status))) %>% 
  ggplot(aes(x=study, y=value, fill=HIV_status)) + 
    geom_boxplot(outlier.shape = NA, alpha=0.8) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) +
    geom_hline(yintercept = cutoff) + 
    theme_bw() + theme(panel.grid.minor = element_blank()) + 
    scale_fill_manual(values=c("#7F7776", "#007C92", "#D3D3D3")) + 
    xlab('') + ylab('Model prediction score')
ggsave(g, filename='~/Desktop/hiv_other_datasets.png',
       width = 7, height = 4)


df.other %>% 
  left_join(meta.others, by='sample_id') %>% 
  group_by(study, HIV_status, positive) %>% 
  tally() %>% 
  group_by(study, HIV_status) %>% 
  mutate(n.all=sum(n)) %>% 
  mutate(freq=n/n.all) %>% 
  filter(positive)
#' # A tibble: 6 × 6
#' # Groups:   study, HIV_status [6]
#' study           HIV_status   positive     n n.all   freq
#' <chr>           <chr>        <lgl>    <int> <int>  <dbl>
#' 1 Fulcher         HIV-negative TRUE        15    98 0.153 
#' 2 Fulcher         PLWH         TRUE         8    41 0.195 
#' 3 NogueraJulian   HIV-negative TRUE         2    27 0.0741
#' 4 NogueraJulian   PLWH         TRUE        23   127 0.181 
#' 5 Rocafort        unclear      TRUE         2    54 0.0370
#' 6 Zhang_Frontiers PLWH         TRUE        24   143 0.168 

# ##############################################################################
# ART- samples

# alpha
g.alpha.all <- df.alpha %>% 
  mutate(general_site='all') %>% 
  mutate(facet='all') %>% 
  bind_rows(df.alpha %>% 
              mutate(facet='sites')) %>% 
  ggplot(aes(x=general_site, y=inv_simpson, fill=hiv_detail)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), 
              colour='#D3D3D3') + 
  geom_boxplot(alpha=0.6, outlier.shape = NA, col='black') + 
  facet_grid(~facet, scales = 'free', space='free') + 
  scale_fill_manual(values=hiv_med[c(1,3,2)], name='') +
  xlab('') + ylab("Prokaryotic diversity") +
  theme_bw() + theme(panel.grid.minor = element_blank(), 
                     panel.grid.major.x=element_blank(),
                     strip.background = element_blank(), 
                     strip.text = element_blank()) +
  NULL
g.alpha.all$data %>% 
  select(general_sample_id, general_site, hiv_detail, inv_simpson) %>% 
  write_csv('./files/source_data/ED10a.csv')
ggsave(g.alpha.all, filename='./figures/hiv/hiv_alpha_all_newdb.pdf',
       width = 7, height = 4, useDingbats=FALSE)
# test
df.alpha %>% 
  group_by(general_site) %>% 
  bind_rows(df.alpha %>% 
              mutate(facet='sites')) %>% 
  group_map(.f=function(.x, .y){
    res.neg <- coefficients(summary(
      lm(data=.x %>% filter(hiv_detail!='HIV Positive(ART+)'), 
         inv_simpson~hiv_detail)))
    res.pos <- coefficients(summary(
      lm(data=.x %>% filter(hiv_detail!='HIV Negative(Control)'), 
         inv_simpson~hiv_detail)))
    tibble(p.val=c(res.neg[2,4], res.pos[2,4]),
           type=c('neg.art-', 'art-.art+'),
           .y)
  }) %>% bind_rows()
#' 
#' # A tibble: 6 × 3
#' p.val type      general_site
#' <dbl> <chr>     <fct>       
#' 1 0.0210  neg.art-  Agincourt   
#' 2 0.803   art-.art+ Agincourt  
#' 3 0.787   neg.art-  Soweto  
#' 4 0.440   art-.art+ Soweto            
#' 5 0.0235  neg.art-  Nairobi    
#' 6 0.00694 art-.art+ Nairobi  

# combined
summary(lmer(inv_simpson~hiv_detail + (1|general_site), 
             data=df.alpha %>% 
               filter(hiv_detail!='HIV Positive(ART+)')))
#'                              Estimate Std. Error       df t value Pr(>|t|)   
#'  (Intercept)                   39.073      1.743   1.811  22.417  0.00318 **
#'  hiv_detailHIV Positive(ART-)  -3.241      4.271 744.976  -0.759  0.44828   

summary(lmer(inv_simpson~hiv_detail + (1|general_site), 
             data=df.alpha %>% 
               filter(hiv_detail!='HIV Negative(Control)')))
#'                              Estimate Std. Error      df t value Pr(>|t|)   
#' (Intercept)                   36.567      4.934   6.646   7.411 0.000191 ***
#' hiv_detailHIV Positive(ART+)  -3.714      4.270 154.936  -0.870 0.385811    


g.alpha.all.phage <- df.alpha.phage %>% 
  mutate(general_site='all') %>% 
  mutate(facet='all') %>% 
  bind_rows(df.alpha.phage %>% 
              mutate(facet='sites')) %>% 
  ggplot(aes(x=general_site, y=richness, fill=hiv_detail)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), 
              colour='#D3D3D3') + 
  geom_boxplot(alpha=0.6, outlier.shape = NA, col='black') + 
  facet_grid(~facet, scales = 'free', space='free') + 
  scale_fill_manual(values=hiv_med[c(1,3,2)], name='') +
  xlab('') + ylab("Phage richness") +
  theme_bw() + theme(panel.grid.minor = element_blank(), 
                     panel.grid.major.x=element_blank(),
                     strip.background = element_blank(), 
                     strip.text = element_blank()) +
  NULL
g.alpha.all$data %>% 
  select(general_sample_id, richness, hiv_detail, general_site) %>% 
  write_csv('./files/source_data/ED10a2.csv')
ggsave(g.alpha.all.phage, filename='./figures/hiv/hiv_alpha_all_phanta.pdf',
       width = 7, height = 4, useDingbats=FALSE)
# test
df.alpha.phage %>% 
  group_by(general_site) %>% 
  bind_rows(df.alpha %>% 
              mutate(facet='sites')) %>% 
  group_map(.f=function(.x, .y){
    res.neg <- coefficients(summary(
      lm(data=.x %>% filter(hiv_detail!='HIV Positive(ART+)'), 
         richness~hiv_detail)))
    res.pos <- coefficients(summary(
      lm(data=.x %>% filter(hiv_detail!='HIV Negative(Control)'), 
         richness~hiv_detail)))
    tibble(p.val=c(res.neg[2,4], res.pos[2,4]),
           type=c('neg.art-', 'art-.art+'),
           .y)
  }) %>% bind_rows()
#' 
#' # A tibble: 6 × 3
#' p.val type      general_site
#' <dbl> <chr>     <fct>       
#' 1 0.104 neg.art-  Agincourt   
#' 2 0.334 art-.art+ Agincourt   
#' 3 0.259 neg.art-  Soweto      
#' 4 0.876 art-.art+ Soweto      
#' 5 0.164 neg.art-  Nairobi     
#' 6 0.144 art-.art+ Nairobi     
# combined
summary(lmer(richness~hiv_detail + (1|general_site), 
             data=df.alpha.phage %>% 
               filter(hiv_detail!='HIV Positive(ART+)')))
#'                              Estimate Std. Error       df t value Pr(>|t|)   
#'  (Intercept)                   402.4648    18.4982   2.0342  21.757  0.00194 **
#'  hiv_detailHIV Positive(ART-)  -0.7271    18.6158 743.6480  -0.039  0.96885   

summary(lmer(richness~hiv_detail + (1|general_site), 
             data=df.alpha.phage %>% 
               filter(hiv_detail!='HIV Negative(Control)')))
#'                              Estimate Std. Error      df t value Pr(>|t|)   
#' (Intercept)                   402.41      29.50   5.71   13.64 1.42e-05 ***
#' hiv_detailHIV Positive(ART+)  -19.06      21.91 154.72   -0.87    0.386    



# effect size correlation
feat.rel <- prop.table(feat.red, 2)
feat.rel <- feat.rel[rowMeans(feat.rel > 1e-04) > 0.05,]
meta.red <- meta.hiv %>% 
  select(general_sample_id, general_site, microbiome_antibiotics, 
         microbiome_diarrhea_last, microbiome_deworming_treatment, 
         microbiome_probiotics, household_size, demographic_age, 
         ses_employment, anthropometric_bmi, card_hypertension_status, 
         card_high_cholesterol, hiv_status) %>% 
  as.data.frame()
rownames(meta.red) <- meta.red$general_sample_id
meta.red$general_sample_id <- NULL

pb <- progress::progress_bar$new(total=nrow(feat.rel))
df.volcano.art <- map(rownames(feat.rel), .f = function(x){
  tmp <- enframe(feat.rel[x,], name='general_sample_id') %>% 
    full_join(meta.hiv, by='general_sample_id') %>% 
    mutate(value=log10(value + 5e-05)) 
  lm.ctr.hiv <- coefficients(summary(
    lm(value~hiv_detail, data=tmp %>% 
         filter(hiv_detail!='HIV Positive(ART-)'))))
  lm.ctr.art <- coefficients(summary(
    lm(value~hiv_detail, data=tmp %>% 
         filter(hiv_detail!='HIV Positive(ART+)'))))
  lm.hiv.art <- coefficients(summary(
    lm(value~hiv_detail, data=tmp %>% 
         filter(hiv_detail!='HIV Negative(Control)'))))
  pb$tick()
  tibble(ctr.hiv=lm.ctr.hiv[2,1],
         ctr.art=lm.ctr.art[2,1],
         art.hiv=lm.hiv.art[2,1],
         species=x)
}) %>% bind_rows()


g.gfc <- df.volcano %>% 
  select(species, p.adj) %>% 
  distinct() %>% 
  full_join(df.volcano.art, by='species') %>% 
  mutate(type=case_when(p.adj < 1e-05~'hiv',
                        p.adj < 1e-01~'hiv2',
                        TRUE~'rest')) %>% 
  arrange(desc(p.adj)) %>% 
  mutate(label=case_when((abs(ctr.art) > 0.3  & abs(ctr.hiv) < 0.3)~species, 
                         p.adj < 1e-05 ~ species,
                         TRUE~'')) %>% 
  ggplot(aes(x=ctr.hiv, y=ctr.art, col=type)) + 
    geom_abline(slope = 1, intercept = 0) +
    geom_point() + 
    ggrepel::geom_text_repel(aes(label=label)) + 
    xlab('gFC control vs HIV (ART+)') + 
    ylab("gFC control vs HIV (ART-)") + 
    theme_bw() + theme(panel.grid.minor = element_blank()) + 
    scale_color_manual(values=c('#007C92', '#AABEC6', '#7F7776'))
g.gfc$data %>% 
  write_csv('./files/source_data/ED10b.csv')
ggsave(g.gfc, filename='./figures/hiv/hiv_art_effect_size_cor.pdf',
       width = 7, height = 5, useDingbats=FALSE)

# model application
load('./files/hiv_models_newdb.RData')
model <- sc.list$all
cutoff.idx <- which(model@eval_data$roc$specificities >= 0.95)[1]
cutoff <- model@eval_data$roc$thresholds[cutoff.idx]
sc.test <- siamcat(feat=feat.rel)
sc.test <- make.predictions(model, sc.test)
df.pred <- enframe(rowMeans(pred_matrix(sc.test)), 
                   name='general_sample_id', value='pred') %>% 
  full_join(meta.hiv) %>% 
  filter(hiv_detail=='HIV Positive(ART-)') %>% 
  bind_rows(enframe(rowMeans(pred_matrix(model)), 
                    name='general_sample_id', value='pred') %>% 
              full_join(meta.hiv) %>% 
              filter(hiv_detail!='HIV Positive(ART-)')) 
df.pred %>% group_by(hiv_detail) %>% 
  group_map(.f=function(.x, .y){tibble(.y, pr=mean(.x$pred > cutoff))}) %>% 
  bind_rows()
g.pred <- df.pred %>% 
  ggplot(aes(x=hiv_detail, y=pred, fill=hiv_detail)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.1, aes(colour=hiv_detail)) +
    geom_hline(yintercept = cutoff) + 
    xlab('') + ylab('Model predictions') + 
    theme_bw() + theme(panel.grid.minor = element_blank()) + 
    theme(panel.grid.major.x = element_blank()) +
    scale_colour_manual(values=hiv_med, guide='none') +
    scale_fill_manual(values=hiv_med, guide='none')
ggsave(g.pred, filename=here('figures', 'hiv', 'hiv_art_pred.pdf'),
       width = 5, height = 5, useDingbats=FALSE)
g.pred$data %>% 
  select(general_sample_id, pred, hiv_detail) %>% 
  write_csv('./files/source_data/ED10c.csv')

# together with the other studies
df.pred.everything <- df.other %>% 
  left_join(meta.others, by='sample_id') %>% 
  mutate(study=case_when(study=='Zhang_Frontiers'~'Zhang', 
                         TRUE~study)) %>% 
  bind_rows(df.pred %>%
              transmute(study='AWI-Gen 2', value=pred, 
                        sample_id=general_sample_id, 
                        positive=pred>cutoff,
                        HIV_status=hiv_detail, 
                        Country='Africa'))
g.pred.all <- df.pred.everything %>% 
  # bind_rows(enframe(rowMeans(model@pred_matrix), name='sample_id') %>% 
  # left_join(meta.hiv %>% transmute(sample_id=general_sample_id, HIV_status=hiv_status))) %>% 
  ggplot(aes(x=study, y=value, fill=HIV_status)) + 
  geom_boxplot(outlier.shape = NA, alpha=0.8) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1)) +
  geom_hline(yintercept = cutoff) + 
  theme_bw() + theme(panel.grid.minor = element_blank()) + 
  scale_fill_manual(values=c("#7F7776", '#8C1515', "#007C92", '#007C92', 
                             "#7F7776", "#D3D3D3")) + 
  xlab('') + ylab('Model prediction score')
ggsave(g.pred.all, filename='./figures/hiv/hiv_predictions_all_combined.pdf',
       width = 7, height = 5)

df.pred.everything %>% 
  group_by(study, HIV_status) %>% 
  summarise(m=mean(positive))

#' 
#' # A tibble: 11 × 3
#' # Groups:   study [6]
#' study         HIV_status                 m
#' <chr>         <chr>                  <dbl>
#' 1 AWI-Gen 2     HIV Negative(Control) 0.0499
#' 2 AWI-Gen 2     HIV Positive(ART+)    0.519 
#' 3 AWI-Gen 2     HIV Positive(ART-)    0.393 
#' 4 Bai           HIV-negative          0     
#' 5 Bai           PLWH                  0     
#' 6 Fulcher       HIV-negative          0.153 
#' 7 Fulcher       PLWH                  0.195 
#' 8 NogueraJulian HIV-negative          0.0741
#' 9 NogueraJulian PLWH                  0.181 
#' 10 Rocafort      unclear               0.0370
#' 11 Zhang         PLWH                  0.168 
