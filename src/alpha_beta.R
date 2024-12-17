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
load('./data/classification/new_db_tables.RData')

# ##############################################################################
# main stuff
set.seed(1951)

tbl.bacteria <- motus.gtdb.lvls$species
meta.red <- metadata.clean %>% 
  filter(meta_site_comparison) %>% 
  filter(general_sample_id %in% colnames(tbl.bacteria))

tbl.bacteria <- tbl.bacteria[,meta.red$general_sample_id]

meta.red %>% 
  select(general_sample_id, general_study_id, general_site) %>% 
  write_tsv(file='./files/site_comp_metadata_red.tsv')

# alpha
alpha.motus <- .f_alpha(tbl.bacteria, l=7500, rarefy=TRUE)
alpha.motus$plot.inv$data %>% 
  select(Sample_ID, site.city, inv_simpson) %>% 
  write_csv('./files/source_data/Fig1c.csv')
ggsave(alpha.motus$plot, filename = here('figures', 'classification_analyses',
                                         'alpha_motus_newdb.pdf'),
       width = 5, height = 8, useDingbats=FALSE)

# beta
x <- .f_beta(t(vegan::rrarefy(t(tbl.bacteria), 7500)), dist='bray')
x$coords %>% 
  write_csv('./files/source_data/Fig1b.csv')
ggsave(x$plot, filename = here('figures', 'classification_analyses',
                               'beta_motus_newdb.pdf'),
       width = 6, height = 6, useDingbats=FALSE)

# correlation with phyla and alpha diversity
df.alpha <- alpha.motus$alpha_df
df <- as_tibble(prop.table(motus.gtdb.lvls$phylum, 2), rownames='phylum') %>% 
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
                        'phylum_pco_correlation_newdb.pdf'),
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
g$data %>% 
  write_csv('./files/source_data/ED4a.csv')
ggsave(g, filename=here('figures','classification_analyses',
                        'phylum_pco_correlation_plots_newdb.pdf'),
       width = 12, height = 6, useDingbats=FALSE)

prop.table(as.matrix(motus.gtdb.lvls$phylum),2)[
  sp.cor %>% 
    select(-calpha) %>% 
    filter(!str_detect(phylum, 'Not_annotated|unassigned')) %>% 
    mutate(name=case_when((abs(c1) > 0.5 | abs(c2) > 0.5)~phylum, 
                          TRUE~'')) %>% 
    filter(name!='') %>% pull(name),] %>% 
  as_tibble(rownames='phylum') %>% 
  pivot_longer(-phylum, names_to = 'general_sample_id') %>% 
  right_join(meta.red, by='general_sample_id') %>% 
  ggplot(aes(x=general_site, y=log10(value+1e-04), fill=general_site)) + 
    geom_boxplot() + 
    facet_grid(~phylum) + 
    scale_fill_manual(values=site.colours, guide='none') + 
    theme_bw() + theme(panel.grid = element_blank(), 
                       axis.ticks.x = element_blank(),
                       axis.text.x = element_blank()) 

df.bar <- prop.table(as.matrix(motus.gtdb.lvls$phylum),2) %>% 
  as_tibble(rownames='phylum') %>% 
  pivot_longer(-phylum, names_to = 'general_sample_id') %>%  
  right_join(meta.red, by='general_sample_id')
df.bar.order <- df.bar %>% 
  filter(phylum=='p__Firmicutes_A') %>% 
  arrange(desc(value))
df.bar.select <- df.bar %>% 
  group_by(phylum) %>% 
  summarise(med_ab=median(value)) %>% 
  arrange(desc(med_ab)) %>% 
  filter(phylum!='unassigned') %>% 
  slice_head(n=10)
g <-df.bar %>% 
  mutate(general_sample_id=factor(general_sample_id, 
                                  levels=df.bar.order$general_sample_id)) %>% 
  mutate(phylum=case_when(phylum%in%df.bar.select$phylum~phylum,
                          TRUE~'other')) %>% 
  mutate(phylum=factor(phylum, levels=c(df.bar.select$phylum, 'other'))) %>% 
  ggplot(aes(x=general_sample_id, y=value, fill=phylum)) + 
    geom_bar(stat='identity') + 
    facet_grid(~general_site, space='free', scales='free') + 
    theme_bw() + theme(panel.grid = element_blank(), 
                       axis.ticks.x = element_blank(),
                       axis.text.x = element_blank(),
                       strip.background = element_blank()) +
    scale_fill_manual(values=RColorBrewer::brewer.pal(11, 'Paired'),
                      name='Phylum') + 
  xlab('') + ylab('Relative abundance')


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
g$data %>% 
  select(general_sample_id, phage_richness, general_site) %>% 
  write_csv('./files/source_data/Fig3c.csv')
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
                                  'phage_bac_alpha_scatter_newdb.pdf'),
       width = 5, height = 5, useDingbats=FALSE)

df.combine %>% 
  group_by(general_site) %>% 
  group_map(.f=function(.x,.y){
    tibble(c=cor(.x$phage_richness, .x$richness, method='spearman'),
           .y)
  }) %>% bind_rows()

#' # A tibble: 6 × 2
#' c general_site
#' <dbl> <fct>       
#' 1 0.287 Nanoro      
#' 2 0.320 Navrongo    
#' 3 0.400 Dimamo      
#' 4 0.272 Agincourt   
#' 5 0.498 Soweto      
#' 6 0.470 Nairobi   

# stats
kruskal.test(richness~general_site, data=df.combine)$p.value # 1.177042e-159
kruskal.test(phage_richness~general_site, data=df.combine)$p.value # 2.506611e-22

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
                        'phage_bac_alpha_comp_newdb.pdf'),
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

#'  A tibble: 4 × 4
#' p.val       c name_read type          
#' <dbl>   <dbl> <chr>     <chr>         
#' 1 1.05e- 26  0.249  all_reads phage_richness
#' 2 6.51e-  2 -0.0435 all_reads richness      
#' 3 2.23e-136  0.540  assembly  phage_richness
#' 4 0          0.768  assembly  richness 

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
g.phage.cor$data %>% 
  write_csv('./files/source_data/ED3b.csv')
g.alpha.cor$data %>% 
  select(general_sample_id, phage_richness, richness, general_site) %>% 
  write_csv('./files/source_data/ED9d.csv')
g.phanta$data %>% 
  select(general_sample_id, phanta_richness, general_site) %>% 
  write_csv('./files/source_data/ED9e.csv')
ggsave(g.sub, filename=here('figures', 'classification_analyses',
                            'alpha_supplement_newdb.pdf'),
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
g.cor$data %>% 
  write_csv('./files/source_data/ED5a.csv')

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
g$data %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  write_csv('./files/source_data/ED5b.csv')
ggsave(g, filename=here('figures', 'general', 'metadata_beta_assoc_newdb.pdf'),
       width = 6, height = 8, useDingbats=FALSE)

write_tsv(df.plot, file='./files/var_explained_newdb.tsv')
