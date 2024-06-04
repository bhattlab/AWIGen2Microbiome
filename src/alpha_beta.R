# ##############################################################################
#
## Alpha and beta diversity plots
#
# ##############################################################################

library('tidyverse')
library('here')
library("yaml")
library("vegan")
library("ggthemes")
library("infotheo")

source(here('src', 'utils.R'))

plot.params <- yaml.load_file(here('files', 'params.yml'))
site.dict <- unlist(plot.params$sites)
site.colours <- unlist(plot.params$sites_colours)

load(here('data', 'metadata', 'metadata_clean.RData'))
load(here('data', 'classification', 'all_classification_tables.RData'))

# ##############################################################################
# main stuff
set.seed(1951)

tbl.bacteria <- lst.motus.gtdb$species
meta.red <- metadata.clean %>% 
  filter(meta_site_comparison) %>% 
  filter(general_sample_id %in% colnames(tbl.bacteria))

tbl.bacteria <- tbl.bacteria[,meta.red$general_sample_id]

meta.red %>% 
  select(general_sample_id, general_study_id, general_site) %>% 
  write_tsv(file='./files/site_comp_metadata_red.tsv')

# alpha
alpha.motus <- .f_alpha(tbl.bacteria, l=7500, rarefy=TRUE)
ggsave(alpha.motus$plot, filename = here('figures', 'classification_analyses',
                                         'alpha_motus.pdf'),
       width = 5, height = 8, useDingbats=FALSE)

# beta
x <- .f_beta(t(vegan::rrarefy(t(tbl.bacteria), 7500)), dist='bray')
ggsave(x$plot, filename = here('figures', 'classification_analyses',
                               'beta_motus.pdf'),
       width = 6, height = 6, useDingbats=FALSE)

# correlation with phyla and alpha diversity
df.alpha <- alpha.motus$alpha_df
df <- as_tibble(prop.table(lst.motus.gtdb$phylum, 2), rownames='phylum') %>% 
  pivot_longer(-phylum, names_to = 'Sample', values_to = 'value') %>% 
  full_join(x$coords, by='Sample') %>% 
  filter(!is.na(PCo1)) %>% 
  mutate(value=log10(value+1e-04)) %>% 
  full_join(df.alpha %>% rename(Sample=Sample_ID), by='Sample')

sp.cor <- df %>% 
  group_by(phylum) %>% 
  group_map(.f=function(.x, .y){
    tibble(c1=cor(.x$PCo1, .x$value, method='spearman'),
           c2=cor(.x$PCo2, .x$value, method='spearman'),
           calpha=cor(.x$richness, .x$value, method='spearman'),
           .y)}) %>% 
  bind_rows()

g <- sp.cor %>% select(-calpha) %>% 
  filter(!str_detect(phylum, 'Not_annotated|unassigned')) %>% 
  mutate(name=case_when((abs(c1) > 0.5 | abs(c2) > 0.5)~phylum, TRUE~'')) %>% 
  add_row(c1=cor(df$richness, df$PCo1, method='spearman'),
          c2=cor(df$richness, df$PCo2, method='spearman'),
          name='alpha', phylum='alpha') %>% 
  mutate(colouring=name=='') %>% 
  arrange(desc(colouring)) %>% 
  ggplot(aes(x=0, y=0, xend=c1, yend=c2)) + 
    geom_segment(aes(group=phylum, col=colouring)) + 
    xlim(-1, 1) + ylim(-1,1) + 
    ggrepel::geom_text_repel(aes(label=name, x=c1, y=c2)) + 
    theme_bw() + 
    theme(panel.grid=element_blank()) + 
    scale_colour_manual(values=c('#007C92', '#D3D3D3')) + 
    xlab('Correlation with PCo1') + 
    ylab('Correlation with PCo2')
ggsave(g, filename=here('figures','classification_analyses',
                        'phylum_pco_correlation.pdf'),
       width = 5, height = 5, useDingbats=FALSE)

g <- df %>% 
  filter(phylum%in% 
           c(sp.cor %>% select(-calpha) %>% 
               filter(!str_detect(phylum, 'Not_annotated|unassigned')) %>% 
               mutate(name=case_when((abs(c1) > 0.5 | abs(c2) > 0.5)~phylum, 
                                     TRUE~'')) %>% 
               filter(name!='') %>% pull(name))) %>% 
  mutate(site.city=factor(site.city, levels=names(site.colours))) %>% 
  select(phylum, Sample, value, PCo1, PCo2, site.city) %>% 
  pivot_longer(c(PCo1, PCo2), names_to='pc', values_to='pc_value') %>% 
  ggplot(aes(x=pc_value, y=value, col=site.city)) + 
    geom_point() + 
    theme_bw() + theme(panel.grid.minor = element_blank()) + 
    scale_color_manual(values=site.colours, name='Site') +
    facet_grid(pc~phylum) + ylab('log10 rel. ab') + xlab('PCo value')
ggsave(g, filename=here('figures','classification_analyses',
                        'phylum_pco_correlation_plots.pdf'),
       width = 12, height = 6, useDingbats=FALSE)

# phage richness
cluster.loc <- paste0('/Volumes/lab_asbhatt/data/bhatt_lab_sequencing/',
                      '23-03-08_awigen2/phageannotation/blast/')

phage.clusters <- read_tsv(paste0(cluster.loc, 'my_clusters.tsv'), 
                           col_names = c('cluster_rep', 'cluster_members'), 
                           col_types = cols())
phage.table <- phage.clusters %>% 
  mutate(cluster_members=str_split(cluster_members, ',')) %>% 
  unnest(cluster_members) %>% 
  mutate(cluster_members=str_remove(cluster_members, '_k141_.*')) %>% 
  mutate(count=1) %>% distinct %>% 
  pivot_wider(names_from = cluster_members, values_from = count, 
              values_fill = 0)

phage.table <- as.data.frame(phage.table)
rownames(phage.table) <- phage.table$cluster_rep
phage.table$cluster_rep <- NULL

df.alpha.phage <- enframe(colSums(phage.table), value='phage_richness', 
                          name='general_sample_id')
g <- df.alpha.phage %>% 
  right_join(meta.red, by='general_sample_id') %>% 
  ggplot(aes(x=general_site, y=phage_richness)) + 
  geom_jitter(width = 0.2, colour='#D3D3D3') +
  geom_boxplot(aes(fill=general_site), outlier.shape = NA) +
  scale_fill_manual(values=alpha(site.colours, 0.85), 
                    name='Site', guide='none') + 
  scale_colour_manual(values=site.colours, name='Site', guide='none') + 
  xlab('') + theme_bw() + ylab("Phage richness") +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank())
ggsave(g, filename=here('figures','classification_analyses',
                        'alpha_phage.pdf'),
       width = 4, height = 4, useDingbats=FALSE)

# compare to bacterial richness
df.combine <- full_join(df.alpha.phage, 
          df.alpha %>% rename(general_sample_id=Sample_ID),
          by='general_sample_id') %>% 
  right_join(meta.red, by='general_sample_id') 
g.alpha.cor <- df.combine %>% 
  ggplot(aes(x=richness, y=phage_richness, col=general_site)) + 
    geom_point() +
    scale_colour_manual(values=site.colours, name='Site', guide='none') + 
    theme_bw() +
    geom_smooth(method = 'lm', se=FALSE)
ggsave(g.alpha.cor, filename=here('figures','classification_analyses',
                                  'phage_bac_alpha_scatter.pdf'),
       width = 5, height = 5, useDingbats=FALSE)

df.combine %>% 
  group_by(general_site) %>% 
  group_map(.f=function(.x,.y){
    tibble(c=cor(.x$phage_richness, .x$richness, method='spearman'),
           .y)
  }) %>% bind_rows()

# stats
kruskal.test(richness~general_site, data=df.combine)$p.value
kruskal.test(phage_richness~general_site, data=df.combine)$p.value

# site-to-site comparison
df.res <- tibble(site1=character(0), site2=character(0), 
                 p=double(0), ef=double(0), type=character(0))
for (i in seq_along(levels(df.combine$general_site))){
  i.site <- levels(df.combine$general_site)[i]
  for (j in i:length(levels(df.combine$general_site))){
    j.site <- levels(df.combine$general_site)[j]
    if (i==j) next()
    df.test.tmp <- df.combine %>% 
      filter(general_site %in% c(i.site,j.site)) 
    fit <- lm(richness~general_site, data=df.test.tmp)
    tmp <- coefficients(summary(fit))
    df.res <- df.res %>% 
      add_row(site1=i.site, site2=j.site, p=tmp[2,4], ef=tmp[2,1]/10,
              type='bac')
    fit <- lm(phage_richness~general_site, data=df.test.tmp)
    tmp <- coefficients(summary(fit))
    df.res <- df.res %>% 
      add_row(site1=j.site, site2=i.site, p=tmp[2,4], ef=-tmp[2,1], 
              type='phage')
  }
}
g <- df.res %>% 
  mutate(site1=factor(site1, levels=levels(df.combine$general_site))) %>% 
  mutate(site2=factor(site2, levels=levels(df.combine$general_site))) %>% 
  mutate(p.adj = p.adjust(p, method='fdr')) %>% 
  mutate(p.nice=sprintf(fmt='%.2E', p.adj)) %>% 
  mutate(p.label=case_when(p.adj < 1e-10~'***', p.adj < 1e-05~'**', 
                           p.adj<1e-1~'*', TRUE~'')) %>% 
  ggplot(aes(x=site1, y=site2, fill=ef)) + 
    geom_abline(slope = 1, intercept = 0, colour='#D3D3D3') +
    geom_tile() + 
    geom_text(aes(label=p.label)) + 
    scale_fill_gradient2(low='#3B528BFF', high='#27AD81FF', mid = 'white', 
                       limits=c(-31, 31)) + 
    theme_minimal() + theme(panel.grid = element_blank()) + xlab('') + 
    ylab('')

ggsave(g, filename=here('figures','classification_analyses',
                        'phage_bac_alpha_comp.pdf'),
       width = 5, height = 4, useDingbats=FALSE)

# ##############################################################################
# double-check that the phage stuff is real

df.stats <- read_tsv(here('data','misc','read_counts.tsv'))
quast.rep <- paste0('/Volumes/lab_asbhatt/data/bhatt_lab_sequencing/',
                    '23-03-08_awigen2/assembly/')
df.assembly.stats <- read_tsv(paste0(quast.rep, 'quast_report.tsv'),
                              col_types = cols()) %>% 
  pivot_longer(-stats, names_to = 'Sample', values_to = 'assembly') %>% 
  filter(stats=='Total length') %>% 
  select(-stats)

df.tmp <- df.combine %>% 
  left_join(df.stats %>% rename(general_sample_id=Sample), 
            by='general_sample_id') %>% 
  left_join(df.assembly.stats %>% rename(general_sample_id=Sample),
            by='general_sample_id') %>% 
  mutate(all_reads=hostremoved_reads + orphan_reads) %>% 
  select(general_sample_id, phage_richness, assembly, all_reads, richness) %>% 
  pivot_longer(c(assembly, all_reads), values_to = 'value_read', 
               names_to = 'name_read') %>% 
  pivot_longer(c(phage_richness, richness), 
               values_to = 'rich', names_to = 'type')
  
g.phage.cor <- df.tmp %>% 
  ggplot(aes(x=log10(value_read), y=rich)) + 
    geom_point(alpha=0.3) +
    facet_grid(type~name_read, scales='free') + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) + 
    geom_smooth(method='lm') + 
    ggpubr::stat_regline_equation(label.y = 100, 
                                  aes(label = after_stat(rr.label)))

df.tmp %>% 
  group_by(name_read, type) %>% 
  group_map(.f=function(.x, .y){
    t.cor <- cor.test(.x$value_read, .x$rich, method='spearman')
    tibble(p.val=t.cor$p.value, c=t.cor$estimate, .y)
    }) %>% 
  bind_rows()

# also check phanta
tbl.phanta <- lst.phanta.vir$species
g.phanta <- enframe(colSums(tbl.phanta > 1e-04), value='phanta_richness', 
             name='general_sample_id') %>% 
  right_join(meta.red, by='general_sample_id') %>% 
  ggplot(aes(x=general_site, y=phanta_richness)) + 
  geom_jitter(width = 0.2, colour='#D3D3D3') +
  geom_boxplot(aes(fill=general_site), outlier.shape = NA) +
  scale_fill_manual(values=alpha(site.colours, 0.85), 
                    name='Site', guide='none') + 
  scale_colour_manual(values=site.colours, name='Site', guide='none') + 
  xlab('') + theme_bw() + ylab("Phage richness") +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank())

g.bottom <- cowplot::plot_grid(g.alpha.cor, g.phanta)
g.sub <- cowplot::plot_grid(g.phage.cor, g.bottom, ncol=1, 
                            rel_heights = c(0.6, 0.4))
ggsave(g.sub, filename=here('figures', 'classification_analyses',
                            'alpha_supplement.pdf'),
       width = 8, height = 9, useDingbats=FALSE)

# ##############################################################################
# metadata stuff

# this is only for SiteComparison, so we filter the master metadata df
df.meta <- metadata.clean %>% 
  filter(meta_site_comparison)

# selection criteria: (number of NAs and so on)

# minimum number of non-NA values per column (at least 100 non-NA values)
# minimum number of unique values (at least two groups)
# minimum variation in each group (at least entropy of 0.1)
# maximum number of unique values (filter out general_study_id, etc.)
# remove meta-variabels (those are only for us)

info.metadata <- map(colnames(df.meta), .f = function(x){
  n.non.na <- sum(!is.na(df.meta[[x]]))
  n.unique <- length(na.omit(unique(df.meta[[x]])))
  data.type <- typeof(df.meta[[x]])
  if (data.type == 'double'){
    ent <- NA_real_
  } else {
    ent <- entropy(df.meta[[x]])
  }
  return(tibble(variable=x, n.non.na=n.non.na, n.unique=n.unique, 
                entropy=ent, type=data.type))}) %>% bind_rows()

info.metadata.red <- info.metadata %>% 
  filter(type != 'character') %>% 
  filter(!str_detect(variable, 'meta_')) %>% 
  filter(n.unique >= 2) %>% 
  filter(n.non.na >= 100) %>% 
  filter(entropy >= 0.2 | is.na(entropy)) %>% 
  filter(variable != 'general_enrollment_date')


tmp <- metadata.clean[,info.metadata.red$variable] %>% 
  mutate_all(as.numeric)
x <- cor(tmp,use='pairwise.complete.obs', method='pearson')
g.cor <- x %>% 
  as_tibble(rownames='var1') %>% 
  pivot_longer(-var1) %>% 
  mutate(var1=case_when(var1%in%c('insulin', 'glucose')~paste0('anth_', var1),
                        TRUE~var1)) %>% 
  mutate(name=case_when(name%in%c('insulin', 'glucose')~paste0('anth_', name),
                        TRUE~name)) %>% 
  mutate(label=case_when(var1==name~'', 
                         abs(value) > 0.8~'*',
                         TRUE~'')) %>% 
  separate(var1, into=c('type1', 'var1'), sep='_', extra='merge') %>% 
  separate(name, into=c('type2', 'var2'), sep='_', extra='merge') %>% 
  mutate(type1=factor(type1, levels = rev(sort(unique(type1))))) %>% 
  
  ggplot(aes(x=var1, y=var2, fill=value)) + 
  geom_tile() + 
  theme_classic() +
  scale_y_discrete(position = 'right') +
  facet_grid(type2~type1, scales = 'free', space='free', switch='x') + 
  scale_fill_gradient2_tableau(palette='Red-Blue-White Diverging',
                               trans='reverse', limits=c(1,-1), 
                               name="Pearson's r", na.value='white') +
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
        axis.ticks = element_blank()) +
  xlab('') + ylab('') +
  geom_text(aes(label=label))

ggsave(g.cor, filename = here('figures', 'general', 
                              'metadata_association.pdf'),
       width = 12, height = 10, useDingbats=FALSE)  

# variance explained by each variable
df.meta.red <- df.meta %>% 
  filter(general_sample_id %in% colnames(tbl.bacteria)) %>% 
  select(c('general_sample_id', info.metadata.red$variable))

# fix NA issue
for (x in colnames(df.meta.red)){
  if (is_double(df.meta.red[[x]])){
    tmp <- df.meta.red[[x]]
    if (any(is.na(tmp))){
      tmp[is.na(tmp)] <- mean(tmp, na.rm=TRUE)
    }
    df.meta.red[[x]] <- tmp
  } else if (is.factor(df.meta[[x]])){
    tmp <- df.meta.red[[x]]
    if (any(is.na(tmp))){
      tmp <- addNA(tmp)
    }
    if (any(tmp == 'NotApplicable')){
      tmp[tmp=='NotApplicable'] <- NA
    }
    if (any(tmp == 'Not Applicable')){
      tmp[tmp=='Not Applicable'] <- NA
    }
    df.meta.red[[x]] <- tmp
  } else {
    message(x)
  }
}

# get variance explained by each variable (save in a tbl)
bc.dist <- vegan::vegdist(t(tbl.bacteria[,df.meta.red$general_sample_id]))
df.plot <- tibble(var=character(0), var.explained=double(0), type=character(0))
df.meta.red <- as.data.frame(df.meta.red)
rownames(df.meta.red) <- df.meta.red$general_sample_id
df.meta.red$general_sample_id <- NULL
pb <- progress::progress_bar$new(total=ncol(df.meta.red))
for (x in colnames(df.meta.red)){
  form <- paste0('bc.dist ~ ', x)
  res <- vegan::dbrda(formula = as.formula(form), data=df.meta.red)
  tmp <- summary(res)
  df.plot <- df.plot %>% 
    add_row(var=x, 
            var.explained=tmp$constr.chi/tmp$unconst.chi, 
            type='single')
  pb$tick()
}

# select the most variance-explainy variable for cases where we have r>0.8
high.cor <- 
  cor(metadata.clean[,info.metadata.red$variable] %>% 
      mutate_all(as.numeric),
    use='pairwise.complete.obs', method='pearson') %>% 
  as_tibble(rownames='var1') %>% 
  pivot_longer(-var1) %>% 
  filter(var1!=name) %>% 
  filter(abs(value) > 0.8)

# kicking out:
# - lipids_ldl
# - hig_medication
# - microbiome_deworming_treatment
# - microbiome_probiotics_period
# - household_cattle & household_poultry
# - anthropometric_bmi & anthropometric_waist_circumference

# afterwards, do the additive model (until we get to less than 0.5%)
iteration <- 0
formula.base <- 'bc.dist ~ 0 + '
left_variables <- colnames(df.meta.red)
left_variables <- setdiff(
  left_variables, c('lipids_ldl', 'hiv_medication', 'microbiome_probiotics',
                    'microbiome_deworming_treatment', 'household_poultry',
                    'household_cattle', 'anthropometric_bmi',
                    'anthropometric_waist_circumference'))
while (length(left_variables) > 0){
  message('iteration: ', iteration)
  pb <- progress::progress_bar$new(total=length(left_variables))
  var.tmp <- rep(0, length(left_variables))
  names(var.tmp) <- left_variables
  for (x in left_variables){
    form <- paste0(formula.base, x)
    # message(form)
    res <- vegan::dbrda(formula = as.formula(form), data=df.meta.red)
    tmp <- summary(res)
    var.tmp[x] <- tmp$constr.chi/tmp$unconst.chi
    pb$tick()
  }
  
  # adjustments!
  idx.var.selected <- which(var.tmp == max(var.tmp))
  selected.var <- names(var.tmp)[idx.var.selected]
  s <- df.plot %>% filter(type=='iterative') %>% pull(var.explained) %>% sum
  var.explained.added <- var.tmp[idx.var.selected] - s
  df.plot <- df.plot %>% 
    add_row(var=selected.var, 
            var.explained=var.explained.added, type='iterative')
  message("Best variable: ", selected.var,
          ', explaining ', sprintf(fmt='%.2f', 100*var.explained.added), 
          '% of variance')
  left_variables <- left_variables[-which(left_variables == selected.var)]
  formula.base <- paste0(formula.base, selected.var, '+')
  iteration <- iteration + 1
}


g <- df.plot %>% 
  pivot_wider(values_from = var.explained, names_from = type) %>% 
  arrange(desc(iterative), desc(single)) %>% 
  mutate(var=factor(var, levels=rev(var))) %>% 
  pivot_longer(-var) %>% 
  ggplot(aes(x=var, y=value, fill=name)) + 
    geom_col(position = position_dodge()) + 
    coord_flip() + 
    ylab('Variance explained [%]') + xlab('') + 
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       panel.grid.major.y=element_blank()) + 
    scale_fill_manual(values=c('#F28E2B', '#4E79A7'), name='')
ggsave(g, filename=here('figures', 'general', 'metadata_beta_assoc.pdf'),
       width = 6, height = 8, useDingbats=FALSE)

# ##############################################################################
# sanity check: other tax levels

# TODO

# ##############################################################################
# sanity check: other classification methods

# TODO





# ##############################################################################
# old

alpha.motus <- .f_alpha(lst.motus.gtdb$species_full_name[,included.samples], 
                        l=7500, rarefy=TRUE)
ggsave(alpha.motus$plot, filename = here('figures', 'classification_analyses',
                                         'alpha_motus.pdf'),
       width = 5, height = 8, useDingbats=FALSE)

alpha.motus.gtdb <- .f_alpha(lst.motus.gtdb$species[,included.samples], 
                        l=7500, rarefy=TRUE)
ggsave(alpha.motus.gtdb$plot, filename = here('figures', 'classification_analyses',
                                         'alpha_motus_gtdb.pdf'),
       width = 5, height = 8, useDingbats=FALSE)
alpha.mpa4 <- .f_alpha(lst.mpa$species, rarefy=FALSE)
ggsave(alpha.mpa4$plot, filename = here('figures', 'classification_analyses',
                                         'alpha_mpa4.pdf'),
       width = 5, height = 8, useDingbats=FALSE)
# alpha.phanta <- .f_alpha(lst.phanta.tax$species, rarefy=FALSE)
# ggsave(alpha.phanta$plot, filename = here('figures', 'classification_analyses',
#                                          'alpha_phanta.pdf'),
       #width = 5, height = 8, useDingbats=FALSE)
#alpha.phanta.vir <- .f_alpha(lst.phanta.vir$species, rarefy=FALSE)
#ggsave(alpha.phanta.vir$plot, filename = here('figures', 'classification_analyses',
#                                          'alpha_phanta_vir.pdf'),
#       width = 5, height = 8, useDingbats=FALSE)

# ##############################################################################
# compare alpha diversity across tools

.f_pairwise_scatter <- function(df){
  # stopifnot(all(c('motus', 'mpa', 'phanta', 'phanta_vir') %in% colnames(df)))
  bind_rows(df %>% transmute(x=motus, y=mpa) %>% mutate(xt='motus', yt='mpa'),
            df %>% transmute(x=motus, y=gtdb) %>% mutate(xt='motus', yt='gtdb'),
            df %>% transmute(x=mpa, y=gtdb) %>% mutate(xt='mpa', yt='gtdb')) %>% 
    # df %>% transmute(x=motus, y=phanta) %>% mutate(xt='motus', yt='phanta'),
    # df %>% transmute(x=motus, y=phanta_vir) %>% 
      # mutate(xt='motus', yt='phanta_vir'),
    # df %>% transmute(x=mpa, y=phanta) %>% mutate(xt='mpa', yt='phanta'),
    # df %>% transmute(x=mpa, y=phanta_vir) %>% mutate(xt='mpa', yt='phanta_vir'),
    # df %>% transmute(x=phanta, y=phanta_vir) %>% 
      # mutate(xt='phanta', yt='phanta_vir')) %>% 
    ggplot(aes(x=x, y=y)) + 
      geom_abline(slope = 1, intercept = 0) + 
      geom_point() + 
      facet_grid(xt~yt) + 
      theme_bw() + 
      theme(panel.grid.minor = element_blank())
    
}

df.alpha.comp <- alpha.motus$alpha_df %>% 
  pivot_longer(-Sample_ID) %>% mutate(type='motus') %>% 
  bind_rows(alpha.mpa4$alpha_df %>% pivot_longer(-Sample_ID) %>% 
              mutate(type='mpa')) %>% 
  bind_rows(alpha.motus.gtdb$alpha_df %>% pivot_longer(-Sample_ID) %>%
              mutate(type='gtdb')) 
  # bind_rows(alpha.phanta$alpha_df %>% pivot_longer(-Sample_ID) %>% 
              # mutate(type='phanta')) %>% 
  # bind_rows(alpha.phanta.vir$alpha_df %>% pivot_longer(-Sample_ID) %>% 
              # mutate(type='phanta_vir')) 

# pdf(here('figures', 'classification_analyses', 'alpha_comp.pdf'),
    # width = 8, height = 8, useDingbats = FALSE)
for (i in c('inv_simpson', 'richness', 'shannon', 'simpson')){
  print(.f_pairwise_scatter(df.alpha.comp %>% 
    filter(name==i) %>% 
    pivot_wider(names_from = 'type', values_from = 'value')) +
      ggtitle(i))
}
# dev.off()  

# ##############################################################################
# beta
for (d in c('bray', 'jaccard', 'log-euclidean')){
  pdf(here('figures', 'classification_analyses', 
           paste0('beta_motus_', d,'.pdf')),
      width = 6, height = 5, useDingbats = FALSE)
  for (lvl in names(lst.motus)){
    if (d=='log-euclidean'){
      x <- .f_beta(prop.table(lst.motus[[lvl]], 2), dist=d, log.n0=1e-05)
    } else {
      x <- .f_beta(t(vegan::rrarefy(t(lst.motus[[lvl]][,included.samples]), 7500)), dist=d)
    }
    print(x$plot)
    message(lvl)
  }
  dev.off()
}

for (d in c('bray', 'jaccard', 'log-euclidean')){
  pdf(here('figures', 'classification_analyses', 
           paste0('beta_gtdb_', d,'.pdf')),
      width = 6, height = 5, useDingbats = FALSE)
  for (lvl in names(lst.motus)){
    if (d=='log-euclidean'){
      x <- .f_beta(prop.table(lst.motus.gtdb[[lvl]], 2), dist=d, log.n0=1e-05)
    } else {
      x <- .f_beta(t(vegan::rrarefy(t(lst.motus.gtdb[[lvl]][,included.samples]), 7500)), dist=d)
    }
    print(x$plot)
    message(lvl)
  }
  dev.off()
}

for (d in c('bray', 'jaccard', 'log-euclidean')){
  pdf(here('figures', 'classification_analyses', 
           paste0('beta_mpa4_', d,'.pdf')),
      width = 6, height = 5, useDingbats = FALSE)
  for (lvl in names(lst.mpa)){
    x <- .f_beta(lst.mpa[[lvl]], dist=d, log.n0 = 1e-06)
    print(x$plot)
    message(lvl)
  }
  dev.off()
}
# for (d in c('bray', 'jaccard', 'log-euclidean')){
#   pdf(here('figures', 'classification_analyses', 
#            paste0('beta_phanta_', d,'.pdf')),
#       width = 6, height = 5, useDingbats = FALSE)
#   for (lvl in names(lst.phanta.tax)){
#     x <- .f_beta(lst.phanta.tax[[lvl]], dist=d, log.n0 = 1e-06)
#     print(x$plot)
#     message(lvl)
#   }
#   dev.off()
# }
# for (d in c('bray', 'jaccard', 'log-euclidean')){
#   pdf(here('figures', 'classification_analyses', 
#            paste0('beta_phanta_vir_', d,'.pdf')),
#       width = 6, height = 5, useDingbats = FALSE)
#   for (lvl in names(lst.phanta.vir)){
#     x <- .f_beta(lst.phanta.vir[[lvl]], dist=d, log.n0 = 1e-06)
#     print(x$plot)
#     message(lvl)
#   }
#   dev.off()
# }

# ##############################################################################
# top phyla in gtdb?
phyla <- prop.table(lst.motus.gtdb$phylum, 2)

phyla.order <- sort(rowMeans(phyla))
as_tibble(phyla, rownames='phylum') %>% 
  pivot_longer(-phylum) %>% 
  mutate(site=str_remove(name, '[0-9]*$')) %>% 
  mutate(site.city=site.dict[site]) %>% 
  mutate(site.city=factor(site.city, levels = names(site.colours))) %>% 
  mutate(value=log10(value + 1e-03)) %>% 
  mutate(phylum=factor(phylum, levels=names(tail(phyla.order)))) %>% 
  filter(phylum %in% names(tail(phyla.order, n=8))) %>% 
  ggplot(aes(x=site.city, y=value, fill=phylum)) + 
    geom_boxplot()

# difference between sites?











# ##############################################################################
# sourmash
sourmash.dist <- read_csv('./data/sourmash/sourmash_k21.csv')
colnames(sourmash.dist) <- str_remove(colnames(sourmash.dist), 'all_reads_')
colnames(sourmash.dist) <- str_remove(colnames(sourmash.dist), '.fq.abundtrim')
sourmash.dist <- as.matrix(sourmash.dist)
rownames(sourmash.dist) <- colnames(sourmash.dist)

beta.div <- vegan::vegdist(vegan::rrarefy(t(lst.motus.gtdb$species), 7500), 
                           method = 'jaccard')
beta.div <- as.matrix(beta.div)


diag(sourmash.dist) <- NA_real_
sourmash.dist[lower.tri(sourmash.dist)] <- NA_real_

diag(beta.div) <- NA_real_
beta.div[lower.tri(beta.div)] <- NA_real_

df.comp <- sourmash.dist %>% 
  as_tibble(rownames='sample_1') %>% 
  pivot_longer(-sample_1, names_to = 'sample_2', values_to = 'sourmash') %>% 
  filter(!is.na(sourmash)) %>% 
  right_join(as_tibble(beta.div, rownames='sample_1') %>% 
              pivot_longer(-sample_1, names_to = 'sample_2', 
                           values_to = 'beta') %>% 
              filter(!is.na(beta)),
            by=c('sample_1', 'sample_2')) %>% 
  mutate(sourmash=1-sourmash)
g <- df.comp %>% 
  ggplot(aes(x=sourmash, y=beta)) + 
    geom_point(alpha=0.1) + 
    geom_abline(slope = 1, intercept = 0) + 
    theme_bw()
ggsave(g, filename='~/Desktop/test.png', width = 6, height = 6)


df.comp %>% 
  mutate(site1=str_remove(sample_1, '[0-9]{3}$')) %>% 
  mutate(site2=str_remove(sample_2, '[0-9]{3}$')) %>% 
  mutate(type=case_when(site1==site2~site1,
                        TRUE~paste0(site1,'-', site2))) %>% 
  mutate(type2=str_detect(type, '-')) %>% 
  select(type, type2, sourmash, beta) %>% 
  pivot_longer(cols=c(sourmash, beta)) %>% 
  ggplot(aes(x=type, y=value, fill=name)) + 
    geom_boxplot() + 
    facet_grid(~type2, scales = 'free_x', space = 'free')

g <- df.comp %>% 
  mutate(site1=str_remove(sample_1, '[0-9]{3}$')) %>% 
  mutate(site2=str_remove(sample_2, '[0-9]{3}$')) %>% 
  mutate(type=case_when(site1==site2~site1,
                        TRUE~paste0(site1,'-', site2))) %>% 
  filter(!str_detect(type, '-')) %>% 
  pivot_longer(cols=c(sourmash, beta)) %>% 
  mutate(site=site.dict[type]) %>% 
  mutate(site=factor(site, levels=names(site.colours))) %>% 
  ggplot(aes(x=site, y=value, fill=site)) + 
    geom_violin( ) +
    geom_boxplot(outlier.shape = NA,fill='white', width=0.2) +
    facet_wrap(~name, scales='free_y') + 
    theme_bw() + 
    scale_fill_manual(values=site.colours) +
    scale_colour_manual(values=site.colours) +
    xlab('') + ylab("Distance")
ggsave(g, filename='./figures/classification_analyses/distances_sourmash.pdf',
       width = 7, height = 4, useDingbats=FALSE)

tmp <- df.comp %>% 
  mutate(site1=str_remove(sample_1, '[0-9]{3}$')) %>% 
  mutate(site2=str_remove(sample_2, '[0-9]{3}$')) %>% 
  mutate(type=case_when(site1==site2~site1,
                        TRUE~paste0(site1,'-', site2))) %>% 
  filter(!str_detect(type, '-')) 


# ##############################################################################
# which taxa are different across sites?

sites <- c('BF', 'GH', 'DM', 'AG', 'KY', 'SW')

phi <- function(t){  
  # expects: t is a 2 x 2 matrix or a vector of length(4)
  stopifnot(prod(dim(t)) == 4 || length(t) == 4)
  if(is.vector(t)) t <- matrix(t, 2)
  r.sum <- rowSums(t)
  c.sum <- colSums(t)
  total <- sum(r.sum)
  r.sum <- r.sum/total
  c.sum <- c.sum/total
  v <- prod(r.sum, c.sum)
  phi <- (t[1,1]/total - c.sum[1]*r.sum[1]) /sqrt(v)
  names(phi) <- NULL
  return(phi)
}

all.res <- list()
for (lvl in names(lst.motus.gtdb)){
  message(lvl)
  tmp <- prop.table(lst.motus.gtdb[[lvl]], 2) %>% 
    as_tibble(rownames='taxa') %>% 
    pivot_longer(-taxa) %>% 
    mutate(site=str_remove(name, '[0-9]{3}$'))
  
  
  for (s1 in seq_along(sites[-1])){
    for (s2 in seq(from=s1+1, to=length(sites))){
      message(sites[s1], '-', sites[s2])
      x <- tmp %>% 
        filter(site %in% c(sites[s1], sites[s2])) %>% 
        mutate(value=value > 1e-04)
      df.res <- map(unique(x$taxa), .f = function(t){
        tab <- table(x %>% filter(taxa==t) %>% pull(value),
                     x %>% filter(taxa==t) %>% pull(site))
        if (any(dim(tab) ==1)){
          return(tibble(taxon=t, effect.size=NA_real_, p.val=NA_real_))
        } else {
          res <- chisq.test(tab)
          ef <- phi(tab)
          return(tibble(taxon=t, effect.size=ef, p.val=res$p.value))
        }})  %>% bind_rows() %>% 
        mutate(comp=paste0(sites[s1], '-', sites[s2])) %>% 
        mutate(level=lvl)
      all.res[[length(all.res) + 1]] <- df.res
    }
  }
}

df.plot <- all.res %>% 
  bind_rows() %>% 
  filter(!is.na(p.val))
g <- df.plot %>% 
  group_by(level, comp) %>% 
  summarise(m=sum(p.val < 1e-03), .groups = 'drop') %>% 
  full_join(df.plot %>% group_by(level) %>% 
              summarise(n.all=length(unique(taxon))), by='level') %>% 
  mutate(frac=m/n.all) %>% 
  separate(comp, into = c('site1', 'site2'), sep='-') %>% 
  mutate(site1=factor(site1, levels=sites)) %>% 
  mutate(site2=factor(site2, levels=sites)) %>% 
  mutate(level=factor(level, levels=names(lst.motus.gtdb))) %>% 
  filter(level!='kingdom') %>%
  
  ggplot(aes(x=site1, y=site2, fill=frac)) + 
    geom_tile() + 
    facet_wrap(~level) + 
    theme_bw() + theme(panel.grid = element_blank()) +
    scale_fill_gradientn(colours=viridis::viridis(n=20), 
                         name='Fraction of\ndiff. ab. taxa') +
    xlab('') + ylab('')

















lvl <- 'species'
d <- 'bray'
x <- .f_beta(t(vegan::rrarefy(t(lst.motus.gtdb[[lvl]][,included.samples]), 7500)), dist=d)


tmp <- x$coords %>% 
  left_join(prop.table(lst.motus.gtdb$phylum, 2) %>% 
              as_tibble(rownames='phylum') %>% 
              pivot_longer(-phylum, names_to = 'Sample'), by='Sample') %>% 
  mutate(value=log10(value + 1e-04))

tmp.r <- enframe(colSums(prop.table(as.matrix(lst.motus.gtdb$species_full_name), 2)>1e-04),
        name='Sample', value='richness') %>% 
  right_join(x$coords, by='Sample') 
  

df.cor <- tmp %>% 
  group_by(site, phylum) %>% 
  mutate(med=median(value)) %>% 
  group_by(phylum) %>% 
  summarise(correlation2=cor(PCo2, value, method='pearson'), 
            correlation1=cor(PCo1, value, method='pearson'), 
            m=max(value))

df.cor %>% 
  mutate(label=case_when(m > -1~phylum, TRUE~'')) %>% 
  ggplot(aes(x=abs(correlation1), y=abs(correlation2), size=m)) + 
    geom_point() + 
    ggrepel::geom_text_repel(aes(label=label))

tmp %>% 
  ggplot(aes(x=PCo2, y=value)) + 
    geom_point() + 
    facet_wrap(~phylum)


tmp %>% 
  filter(phylum=='f__Rikenellaceae') %>% 
  ggplot(aes(x=PCo2, y=value, colour=site)) + 
    geom_point()


# get the eta for each phylum per site
tmp %>% group_by(phylum) %>% 
  group_map(.f=function(.x, .y){
    res <- anova(lm(value~site, data=.x))
    eta <- res$`Sum Sq`[1]/sum(res$`Sum Sq`)
    tibble(eta=eta, .y)
    }) %>% bind_rows() %>% full_join(df.cor, by='phylum') %>% View
  ggplot(aes(x=correlation1, y=correlation2, colour=eta)) + 
    geom_point() + 
    scale_colour_viridis_c()


  tmp <- bind_cols(
    metadata.clean %>% 
      select(general_sample_id), 
    metadata.clean %>% 
      select_if(is.factor)) %>% 
    left_join(x$coords %>% rename(general_sample_id=Sample))
df.cor <- tibble(pc=character(0), meta=character(0), p=double(0), sp=double(0))
for (x in colnames(metadata.clean %>% select_if(is.factor))){
  pearson <- cor(tmp$PCo1, as.numeric(tmp[[x]]), use = 'pairwise.complete.obs')
  spearman <- cor(tmp$PCo1, as.numeric(tmp[[x]]), use = 'pairwise.complete.obs', method='spearman')
  df.cor <- df.cor %>% 
    add_row(pc='1', meta=x, p=pearson, sp=spearman)
  pearson <- cor(tmp$PCo2, as.numeric(tmp[[x]]), use = 'pairwise.complete.obs')
  spearman <- cor(tmp$PCo2, as.numeric(tmp[[x]]), use = 'pairwise.complete.obs', method='spearman')
  df.cor <- df.cor %>% 
    add_row(pc='2', meta=x, p=pearson, sp=spearman)
  
  
}
  
  