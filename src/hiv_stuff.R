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
load(here('data', 'classification', 'all_classification_tables.RData'))

feat.red <- lst.motus.gtdb$species

meta.hiv <- metadata.clean %>% 
  filter(meta_hiv_comparison) %>% 
  filter(general_site !='Dimamo') %>% 
  filter(general_sample_id %in% colnames(feat.red)) %>% 
  filter(hiv_status %in% c("Yes", 'No')) %>% 
  mutate(hiv_detail=case_when(
    hiv_status=='No' ~ 'HIV Negative(Control)',
    hiv_status=='Yes' & hiv_medication=='Yes' ~ 'HIV Positive(ART+)',
    hiv_status=='Yes' & hiv_medication=='No' ~ 'HIV Positive(ART-)',
    TRUE ~NA_character_)) %>% 
  filter(!is.na(hiv_detail))
feat.red <- feat.red[,meta.hiv$general_sample_id]

hiv_med <- c("#7F7776","#007C92", "#8C1515")

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

ggsave(g, filename=here('figures','hiv','hiv_alpha.pdf'), 
       width = 10, height = 5, useDingbats=FALSE)

# testing
# individual
df.alpha %>% 
  group_by(general_site) %>% 
  bind_rows(df.alpha %>% 
              mutate(facet='sites')) %>% 
  filter(hiv_detail!='HIV Positive(ART-)') %>% 
  group_map(.f=function(.x, .y){
    res <- coefficients(summary(lm(data=.x, richness~hiv_status)))
    tibble(p.val=res[2,4], .y)
    }) %>% bind_rows()

# combined
summary(lmer(richness~hiv_status + (1|general_site), 
             data=df.alpha %>% bind_rows(df.alpha %>% 
                                           mutate(facet='sites')) %>% 
               filter(hiv_detail!='HIV Positive(ART-)')))

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
ggsave(g, filename=here('figures', 'hiv', 'hiv_beta.pdf'),
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
df.plot %>% 
  select(Sample, PCo1, PCo2, hiv_detail) %>% 
  pivot_longer(c(PCo1, PCo2)) %>% 
  group_by(name) %>% 
  group_map(.f=function(.x, .y){
    t <- kruskal.test(data=.x, value~hiv_detail)
    tibble(p.val=t$p.value, .y, type='hiv_detail')
  }) %>% 
  bind_rows()


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
  fit <- suppressMessages(lmer(value~hiv_detail + (1|general_site) + 
                                 (1|microbiome_antibiotics) + (1|microbiome_diarrhea_last), 
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
                       'general_enrollment_date', 'anthropometric_waist_circumference',
                       'anthropometric_hip_circumference', 'card_cholesterol_treatment',
                       'card_diabetes_treatment', 'card_esr_crp', 
                       'card_hypertension_treatment',
                       'card_rheumatoid_factor', 'ultrasound_vat', 'ultrasound_scat',
                       'general_country', 'general_region', 'demographic_gender',
                       'hiv_status', 'hiv_medication', 'lab_fasting_confirmed',
                       'card_diabetes_treatment_spec', 'household_pottable_water', 
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

# export table as supplement
df.volcano %>% select(species, p.value, p.adj, effect.size, 
                      Agincourt, Nairobi, Soweto) %>% 
  distinct() %>% 
  filter(p.adj < 0.01) %>% write_tsv(here('files', 'hiv_table.tsv'))

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


# add a phylum strip to the heatmap?
gtdb.tax.all <- read_tsv('./files/mOTUs_3.0.0_GTDB_tax.tsv', 
                         col_names = c('motus', 'domain', 'phylum', 'class', 
                                       'order', 'family', 'genus', 'species'),
                         col_types = cols()) 
g3 <- df.volcano %>% 
  ungroup() %>% 
  select(species, p.adj, effect.size, Agincourt, Soweto, Nairobi) %>% 
  distinct() %>% 
  filter(p.adj < 0.01) %>% 
  mutate(ordering=sign(effect.size) *log10(p.adj)) %>% 
  select(species, ordering, Agincourt, Soweto, Nairobi) %>% 
  left_join(gtdb.tax.all %>% select(-motus), by='species') %>% distinct() %>%
  arrange(ordering) %>% 
  mutate(species=factor(species, levels=species)) %>% 
  ggplot(aes(x=species, y=1, fill=phylum)) +
  geom_tile() +
  scale_fill_manual(values=c('#FB8072', '#FDB462', '#FFFFB3', '#8DD3C7', 
                             '#8DD3C7', '#8DD3C7', '#BC80BD', '#B3DE69', 
                             '#BEBADA')) + 
  theme_bw() + theme(panel.grid = element_blank(), 
                     axis.text = element_blank(),
                     axis.ticks = element_blank())

g <- cowplot::plot_grid(g1, g2, g3, ncol=1, rel_heights = c(0.4, .4, 0.2))
ggsave(g, filename=here('figures', 'hiv', 'hiv_effect_size_heatmap.pdf'),
       width = 12, height = 6)


# ##############################################################################
# SIAMCAT?
# ML overall & ML per site
set.seed(1791)

library("SIAMCAT")

if (!file.exists(here('files','hiv_models.RData'))){
  
  sc.obj <- siamcat(feat=feat.rel, meta=meta.red, label='hiv_status', case='Yes')
  sc.obj <- filter.features(sc.obj, filter.method='pass', cutoff=0.05)
  # check.confounders(sc.obj, fn.plot = '~/Desktop/confounder.pdf', verbose=3)
  sc.obj <- normalize.features(sc.obj, norm.method='log.std', 
                               norm.param=list(log.n0=5e-05, sd.min.q=0))
  sc.obj <- create.data.split(sc.obj, num.folds = 5, num.resample = 5)
  sc.obj <- train.model(sc.obj, method='lasso', grid.size = 5)
  sc.obj <- make.predictions(sc.obj)
  sc.obj <- evaluate.predictions(sc.obj)
  model.evaluation.plot(sc.obj)
  
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
  save(sc.list, file=here('files', 'hiv_models.RData'))
} else {
  load(here('files', 'hiv_models.RData'))
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
                      fn.plot = here('figures', 'hiv', 'ml_auroc.pdf'))

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
ggsave(g, filename=here('figures', 'hiv', 'ml_model_transfer.pdf'),
       width = 6, height = 5, useDingbats=FALSE)

# apply to other sites?
test.ext <- prop.table(lst.motus.gtdb$species, 2)
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
ggsave(g, filename=here('figures', 'hiv', 'hiv_fpr.pdf'),
       width = 4, height = 6)

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
  full_join(df.volcano %>% select(species, p.value, effect.size, p.adj) %>% distinct()) %>% 
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
  mutate(name=factor(name, levels=c('all.rank', 'Agincourtleft_out', 'Nairobileft_out', 'Sowetoleft_out', 'Agincourt', 'Nairobi', 'Soweto'))) %>% 
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
    scale_fill_gradient2_tableau(palette='Red-Blue-White Diverging', guide='none')


cowplot::plot_grid(g, g2, g2.1, nrow=1, rel_widths = c(0.6, 0.2, 0.2))


# ##############################################################################
# ART- samples

# alpha
hiv_med <- c("#7F7776","#8C1515", "#007C92")
g.alpha.all <- df.alpha %>% 
  mutate(general_site='all') %>% 
  mutate(facet='all') %>% 
  bind_rows(df.alpha %>% 
              mutate(facet='sites')) %>% 
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
# test
df.alpha %>% 
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

# combined
summary(lmer(richness~hiv_detail + (1|general_site), 
             data=df.alpha %>% bind_rows(df.alpha %>% 
                                           mutate(facet='sites')) %>% 
               filter(hiv_detail!='HIV Positive(ART+)')))
summary(lmer(richness~hiv_detail + (1|general_site), 
             data=df.alpha %>% bind_rows(df.alpha %>% 
                                           mutate(facet='sites')) %>% 
               filter(hiv_detail!='HIV Negative(Control)')))

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



# sp <- df.volcano.art %>% mutate(diff=abs(ctr.hiv - ctr.art)) %>% arrange(desc(diff)) %>% filter(diff > 0.4) %>% pull(species)
# sp <- df.volcano.art %>% full_join(df.volcano %>% 
#                                      select(species, p.adj) %>% 
#                                      distinct()) %>%  filter(p.adj < 1e-05) %>% pull(species)
# as_tibble(feat.rel[sp,], rownames='species') %>% 
#   pivot_longer(-species, names_to = 'general_sample_id') %>% 
#   mutate(value=log10(value+1e-04)) %>% 
#   left_join(meta.hiv) %>% 
#   ggplot(aes(x=hiv_detail, y=value)) + 
#     geom_boxplot(outlier.shape = NA) +
#     geom_jitter(width = 0.1) +
#     facet_wrap(~species)

# model application
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

g <- cowplot::plot_grid(g.alpha.all, g.gfc, g.pred)
ggsave(g, filename=here('figures', 'hiv', 'hiv_art_supp.pdf'),
       width = 12, height = 12, useDingbats=FALSE)
# ##############################################################################
# OLD STUFF
# ##############################################################################



# some examples?
tmp <- df.feat.weight %>% group_by(name) %>% 
  arrange(desc(abs(mean.rel.weight))) %>% 
  mutate(rank=row_number(-abs(mean.rel.weight))) %>% 
  left_join(df.volcano %>% select(species, p.adj, effect.size) %>% distinct,
            by='species') %>% 
  mutate(rank.l=case_when(rank<=50~rank, TRUE~NA)) %>% 
  select(species, mean.rel.weight, name, p.adj, effect.size) %>% 
  pivot_wider(names_from = name, values_from = mean.rel.weight) %>% 
  mutate(rank.all=row_number(-abs(all))) %>% View
  filter(rank.all < 20) %>% 
  mutate(rank.all=sign(all)*rank.all) %>% 
  pivot_longer(-c(species, p.adj, effect.size, rank.all)) %>% 
  arrange(rank.all) %>% 
  mutate(species=factor(species, levels=unique(species))) %>% 
  filter(name %in% c('Agincourt', 'Nairobi', 'Soweto', 'all')) %>% 
  mutate(name=factor(name, levels=c('all', 'Agincourt', 'Soweto', 'Nairobi')))
tmp %>% 
  ggplot(aes(x=species, y=value, fill=name)) + 
    geom_bar(stat='identity', position = position_dodge()) + 
    coord_flip()
tmp %>% select(species, p.adj) %>% distinct() %>% 
  ggplot(aes(x=species, y=-log10(p.adj))) + 
    geom_col() + coord_flip()
  
  
  ggplot(aes(x=effect.size, y=mean.rel.weight)) + 
    geom_point() + 
    facet_wrap(~name)
  
  
  filter(rank.l > 20) %>% 
  filter(!str_detect(name, 'left')) %>% 
  mutate(ordering=sign(effect.size) *log10(p.adj)) %>% 
  arrange(ordering) %>%
  mutate(species=factor(species, levels=unique(species))) %>% 
  ggplot(aes(x=species, y=name, fill=mean.rel.weight)) + 
    geom_tile() + 
    scale_fill_gradientn(colours=c('#9e3d22', '#e36621', '#fcad52', '#ffffff', '#91d183', '#539e52', '#24693d')) + 
    theme_bw() + 
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) + 
    geom_text(aes(label=rank.l))
    



  '#9e3d22'



  
  
x <- 's__Dysosmobacter welbionis'
x <- 's__Bifidobacterium ruminantium'
x <- 's__Fusobacterium_A mortiferum'
x <- 's__Limivicinus sp900547315'
x <- 's__UMGS1375 sp900066615'
tmp <- enframe(feat.rel[x,], name='general_sample_id') %>% 
  full_join(meta.hiv, by='general_sample_id') %>% 
  mutate(value=log10(value + 5e-05)) %>% 
  filter(hiv_detail!='HIV Positive(ART-)')
tmp %>% 
  ggplot(aes(x=general_site, y=value, fill=hiv_detail)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1))

tmp %>% 
  ggplot(aes(x=household_toilet, y=value)) + 
  geom_boxplot()

   
  select(species, p.adj, effect.size, Agincourt, Soweto, Nairobi) %>% 
  distinct() %>% 
  filter(p.adj < 0.01) %>% 
  mutate(order=sign(effect.size) *log10(p.adj)) %>% 
  arrange(order) %>% 
  mutate(species=factor(species, levels=species)) %>% 
  select(species, Agincourt, Soweto, Nairobi) %>% 
  pivot_longer(-species) %>% 
  mutate(name=factor(name, levels=c('Agincourt','Soweto', 'Nairobi'))) %>% 
  ggplot(aes(x=species, y=name, fill=value)) + 
  geom_tile() + 
  ggthemes::scale_fill_gradient2_tableau(palette='Red-Blue-White Diverging',
                                         name='Effect size') + 
  theme_bw() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid = element_blank()) + 
  ylab('')




df.volcano %>% 
  mutate(p.adj=p.adjust(p.value, method='fdr')) %>% 
  arrange(p.adj) %>%
  ggplot(aes(x=effect.size, y=-log10(p.adj))) + 
    geom_point()
df.volcano %>% 
  mutate(p.adj=p.adjust(p.value, method='fdr')) %>% 
  arrange(p.adj) %>% 
  filter(p.adj < 1e-05) 

# - variance explained for all meta-variables (take SIAMCAT code?)
# - association between HIV and meta-variables (confounding, also take SIAMCAT code?)

# meta-deconfounder????
# 



g1 <- df.volcano %>% 
  mutate(p.adj=p.adjust(p.value, method='fdr')) %>% 
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
  mutate(p.adj=p.adjust(p.value, method='fdr')) %>% 
  filter(p.adj < 0.01) %>% 
  mutate(order=sign(effect.size) *log10(p.adj)) %>% 
  arrange(order) %>% 
  mutate(species=factor(species, levels=species)) %>% 
  select(species, Agincourt, Soweto, Nairobi) %>% 
  pivot_longer(-species) %>% 
  mutate(name=factor(name, levels=c('Agincourt','Soweto', 'Nairobi'))) %>% 
  ggplot(aes(x=species, y=name, fill=value)) + 
    geom_tile() + 
    ggthemes::scale_fill_gradient2_tableau(palette='Red-Blue-White Diverging',
                                           name='Effect size') + 
    theme_bw() + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid = element_blank()) + 
    ylab('')
g <- cowplot::plot_grid(g1, g2, ncol = 1)
ggsave(g, filename='./figures/hiv/effect_size_heatmap.pdf',width = 8, height = 4)






tmp <- enframe(feat.rel['s__Coprococcus eutactus',], name='general_sample_id') %>% 
  full_join(meta.hiv, by='general_sample_id') %>% 
  mutate(value=log10(value + 1e-04)) %>% 
  filter(hiv_detail!='HIV Positive(ART-)')

tmp %>% 
  ggplot(aes(x=hiv_detail, y=value)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.1) +
    facet_grid(~general_site)



tmp <- map(rownames(feat.rel), .f=function(x){
  df.tmp <- enframe(feat.rel[x,], name='general_sample_id') %>% 
    full_join(meta.hiv, by='general_sample_id') %>% 
    mutate(value=log10(value + 1e-04)) %>% 
    filter(hiv_detail!='HIV Positive(ART-)')
  df.tmp <- df.tmp %>% group_by(general_site) %>% group_map(~coefficients(summary(lm(value~hiv_status, data=.x))))
  matrix(map(df.tmp, .f = function(x){x[2,1]}) %>% unlist()) %>% t() %>% as.data.frame()
}) %>% bind_rows()



















