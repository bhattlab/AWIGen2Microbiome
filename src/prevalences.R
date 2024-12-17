# ##############################################################################
#
## Plot the prevalence across sites for all species of interest
#
# ##############################################################################

library("tidyverse")
library("here")
library("ggthemes")
library("pheatmap")
library("yaml")

plot.params <- yaml.load_file(here('files', 'params.yml'))
site.dict <- unlist(plot.params$sites)
site.colours <- unlist(plot.params$sites_colours)
site.colours <- c(site.colours[1:2], site.colours[6], site.colours[3:5])

# ##############################################################################
# get data
#load(here('data', 'classification', 'all_classification_tables.RData'))
load(here('data', 'metadata', 'metadata_clean.RData'))
load('./data/classification/new_db_tables.RData')

# ##############################################################################
# Differences between sites

tbl.bacteria <- motus.gtdb.lvls$motus

df.meta <- metadata.clean %>% 
  filter(meta_site_comparison) %>% 
  filter(general_sample_id %in% colnames(tbl.bacteria))
tbl.bacteria <- tbl.bacteria[,df.meta$general_sample_id]
tbl.bacteria <- prop.table(tbl.bacteria, 2)
tbl.bacteria <- tbl.bacteria[rowMeans(tbl.bacteria!=0) > 0, ]
# write.csv(tbl.bacteria, '../misc/classification_table_newdb.csv')


# calculate prevalence per site
prev.per.site <- tbl.bacteria %>% as_tibble(rownames='species') %>% 
  pivot_longer(-species) %>% 
  mutate(site=str_remove(name, '[0-9]*$')) %>% 
  group_by(species, site) %>% 
  summarise(prevalence=mean(value>1e-04)) %>% 
  pivot_wider(names_from = site, values_from = prevalence) %>% 
  as.data.frame()
rownames(prev.per.site) <- prev.per.site$species
prev.per.site$species <- NULL

# some filtering
prev.per.site <- prev.per.site[rowSums(prev.per.site > 0.05) > 2,]



# write.csv(prev.per.site, '../misc/prevalence_per_site.csv')

# site correlation by prevalence
g <- cor(prev.per.site, method='spearman') %>% 
  as_tibble() %>% 
  mutate(site1=colnames(prev.per.site)) %>%
  pivot_longer(-site1, names_to='site2') %>% 
  filter(site1!=site2) %>% 
  mutate(site1=factor(site1, c('BF', 'GH', 'KY', 'DM', 'AG', 'SW'))) %>% 
  mutate(site2=factor(site2, c('BF', 'GH', 'KY', 'DM', 'AG', 'SW'))) %>% 
  mutate(label=sprintf(fmt='%.2f', value)) %>% 
  ggplot(aes(x=site1, y=site2, fill=value)) + 
  geom_tile() + 
  scale_fill_gradientn(colours=viridis::viridis(9), limits=c(0,1)) + 
  geom_text(aes(label=label))
g$data %>% 
  write_csv('./files/source_data/ED6a2.csv')
ggsave(g, filename=here('figures','prevalence','site_correlation_newdb.pdf'),
       width = 5, height = 5, useDingbats=FALSE)

# variance/mean prevalence
prev.per.site %>% 
  as_tibble(rownames='species') %>% 
  pivot_longer(-species) %>% 
  group_by(species) %>% 
  summarise(m=mean(value), v=var(value)) %>% 
  mutate(type=m>0.15 & v>0.02) %>% 
  ggplot(aes(x=m, y=v, col=type)) + 
    geom_point()

species.included <- prev.per.site %>% 
  as_tibble(rownames='species') %>% 
  pivot_longer(-species) %>% 
  group_by(species) %>% 
  summarise(m=mean(value), v=var(value)) %>% 
  mutate(type=m>0.15) %>% 
  filter(type) %>% 
  pull(species)

dist <- 'maximum'

clustering <- hclust(dist(prev.per.site[species.included,], 
                          method = dist), 
                     method='ward.D2')
labels <- cutree(clustering, k=5)

prev.per.site %>% 
  as_tibble(rownames='species') %>% 
  full_join(enframe(labels, name='species', value='cluster')) %>% 
  filter(!is.na(cluster)) %>% 
  mutate(species=str_replace(species, 'Incongruent', 'Not_annotated')) %>% 
  write_tsv(here('files', 'prevalences.tsv'))

dev.off()
pdf(here('figures','prevalence','prev_heatmap_gtdb_newdb.pdf'), width = 6, 
    height = 7, useDingbats = FALSE)
pheatmap(prev.per.site[names(labels),c('BF', 'GH', 'KY', 'DM', 'AG', 'SW')],
         breaks=seq(from=0,to=1, length.out=99),
         show_rownames = FALSE, border_color = NA,
         color = colorRampPalette(
           colors = c('white', '#DAD7CB', '#7F7776', '#5F574F'))(101),
         cluster_cols = FALSE, clustering_method = 'ward.D2',
         clustering_distance_rows = dist, cutree_rows = 5)
dev.off()

write_csv(as_tibble(prev.per.site[names(labels),], rownames='species'),
          file = './files/source_data/ED6a.csv')

prev.per.site[names(labels),c('BF', 'GH')] %>% 
  as_tibble(rownames='species') %>% 
  full_join(enframe(labels, name='species', value='label')) %>% 
  ggplot(aes(x=BF, y=GH, col=as.factor(label))) + 
    geom_point() +
    geom_smooth(method='lm')

# for each cluster, prevalence and abundance in each site
g <- prev.per.site %>% 
  as_tibble(rownames='species') %>% 
  left_join(enframe(labels, name='species', value='cluster')) %>% 
  filter(!is.na(cluster)) %>% 
  pivot_longer(-c(species, cluster)) %>% 
  mutate(name=site.dict[name]) %>%
  mutate(name=factor(name, levels=names(site.colours))) %>% 
  ggplot(aes(x=name, y=value, fill=name)) + 
    geom_boxplot(outlier.color = '#D3D3D3') + 
    facet_grid(~cluster) + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major.x=element_blank()) + 
    scale_fill_manual(values=site.colours)
ggsave(g, filename=here('figures','prevalence','cluster_prevalence.pdf'),
       width = 6, height = 4, useDingbats=FALSE)  

prev.per.site %>% 
  as_tibble(rownames='species') %>% 
  left_join(enframe(labels, name='species', value='cluster')) %>% 
  filter(!is.na(cluster)) %>% 
  pivot_longer(-c(species, cluster)) %>% 
  mutate(name=site.dict[name]) %>%
  mutate(name=factor(name, levels=names(site.colours))) %>% 
  select(species, cluster) %>% 
  distinct() %>% 
  left_join(tbl.bacteria %>% 
              as_tibble(rownames='species') %>% 
              pivot_longer(-species) %>% 
              filter(species %in% names(labels))) %>% 
  group_by(cluster, name) %>% 
  summarise(sm=sum(value)) %>% 
  mutate(site=str_remove(name, '[0-9]*$')) %>% 
  mutate(site=site.dict[site]) %>%
  mutate(site=factor(site, levels=names(site.colours))) %>% 
  ggplot(aes(x=site, y=sm, fill=site)) +
    geom_boxplot() +
    facet_grid(~cluster) +
    scale_fill_manual(values=site.colours) +
    theme_bw() + 
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major.x=element_blank())

# over-representation of genera?
gtdb.tax.all <- read_tsv('./files/mOTUs_ext_GTDB_tax.tsv',
                     col_names = c('motus', 'domain', 'phylum', 'class', 
                                   'order', 'family', 'genus', 'species'),
                     col_types = cols()) 

# join prevalence with mean, GTDB classification, and save
tbl.bacteria %>% as_tibble(rownames='species') %>% 
  pivot_longer(-species) %>% 
  mutate(site=str_remove(name, '[0-9]*$')) %>% 
  filter(species %in% rownames(prev.per.site)) %>% 
  group_by(species, site) %>% 
  summarise(prevalence=mean(value>1e-04), mean=mean(log10(value + 1e-05))) %>% 
  mutate(site=site.dict[site]) %>% 
  pivot_wider(names_from = site, values_from = c(prevalence, mean)) %>% 
  mutate(motus=str_extract(species, '\\[(ref|meta|ext|NEW).*\\]')) %>% 
  mutate(motus=str_remove(motus, '\\[')) %>% 
  mutate(motus=str_remove(motus, '\\]')) %>% 
  left_join(gtdb.tax.all %>% transmute(motus, gtdb.species=species), 
            by='motus') %>% 
  select(-motus) %>% 
  write_csv('./files/mean_prev_per_site.csv')


species.included.ids <- species.included %>%
  str_remove('\\]$') %>% 
  str_remove('.*\\[')


gtdb.tax <- gtdb.tax.all %>% 
  filter(motus %in% species.included.ids) %>% 
  filter(str_detect(genus, '^g__')) %>% 
  select(genus, species, motus) %>% 
  distinct()

tax.cluster <- gtdb.tax %>% 
  full_join(enframe(labels, name='motus', value='cluster') %>% 
              mutate(motus=str_extract(motus, '\\[(ref|meta|ext|NEW).*\\]')) %>% 
              mutate(motus=str_remove(motus, '\\[')) %>% 
              mutate(motus=str_remove(motus, '\\]'))) %>% 
  filter(!is.na(genus))


# abundance difference between sites 
feat.rel <- prop.table(tbl.bacteria, 2)
site.gfc <- as_tibble(feat.rel[species.included,], rownames='species') %>% 
  pivot_longer(-species) %>% 
  mutate(value=log10(value+1e-04)) %>% 
  mutate(site=str_remove(name, '[0-9]*$')) %>% 
  group_by(species, site) %>% 
  group_map(.f=function(.x, .y){
    tmp <- quantile(.x$value, probs = seq(from=0.05, to=0.95, by=0.05))
    tibble(quant=tmp, quant_lvl=seq(from=0.05, to=0.95, by=0.05), .y)}) %>% 
  bind_rows() %>% 
  pivot_wider(names_from = site, values_from = quant) %>% 
  group_by(species) %>% 
  summarise(bf_sw=mean(BF-SW), bf_ag=mean(BF-AG), bf_gh=mean(BF-GH), 
            bf_ky=mean(BF-KY), bf_dm=mean(BF-DM), gh_sw=mean(GH-SW), 
            gh_ag=mean(GH-AG), gh_ky=mean(GH-KY), gh_dm=mean(GH-DM),
            ky_dm=mean(KY-DM), ky_ag=mean(KY-AG), ky_sw=mean(KY-SW),
            dm_ag=mean(DM-AG), dm_sw=mean(DM-SW), ag_sw=mean(AG-SW))

site.gfc %>% 
  mutate(motus=str_extract(species, '\\[(ref|meta|ext|NEW).*\\]')) %>% 
  mutate(motus=str_remove(motus, '\\[')) %>% 
  mutate(motus=str_remove(motus, '\\]')) %>% 
  left_join(gtdb.tax.all %>% 
              transmute(motus, gtdb.species=species), by='motus') %>%
  select(-motus) %>% 
  write_csv('./files/fc_per_site_newdb.csv')

# write.csv(site.gfc, '../misc/fc_per_site.csv')

site.gfc %>% 
  pivot_longer(-species) %>% 
  mutate(value=abs(value)) %>% 
  ggplot(aes(x=value, colour=name)) + 
    geom_density()

cutoff <- site.gfc %>% 
  pivot_longer(-species) %>% 
  mutate(value=abs(value)) %>% 
  pull(value) %>% quantile(probs=0.9)

# network!
site.gfc %>% 
  pivot_longer(-species) %>% 
  filter(abs(value)>cutoff) %>% 
  group_by(name) %>% tally() %>% arrange(desc(n))

# heatmap?
g <- site.gfc %>% 
  pivot_longer(-species) %>% 
  filter(abs(value)>cutoff) %>% 
  group_by(name) %>% tally() %>%
  separate(name, into=c('site', 'site2'), sep='_') %>% 
  mutate(site=site.dict[toupper(site)]) %>% 
  mutate(site2=site.dict[toupper(site2)]) %>% 
  mutate(site=factor(site, levels=names(site.colours))) %>% 
  mutate(site2=factor(site2, levels=names(site.colours))) %>% 
  mutate(white=n>100) %>% 
  ggplot(aes(x=site, y=site2, fill=n)) + 
    geom_tile() + 
    theme_bw() + theme(panel.grid=element_blank()) + 
    xlab('') + ylab('') + 
    scale_fill_gradientn(colours=viridis::viridis(12), guide='none') + 
    geom_text(aes(label=n, colour=white)) + 
    scale_colour_manual(values=c('white', 'black'), guide='none')
g$data %>% 
  write_csv('./files/source_data/Fig1d.csv')
ggsave(g, filename='./figures/prevalence/number_of_different_species_newdb.pdf',
       width = 4, height = 3, useDingbats=FALSE)

# select genera to highlight!
genera <- site.gfc %>% 
  pivot_longer(-species) %>% 
  mutate(motus=str_remove(species, '.*\\[')) %>% 
  mutate(motus=str_remove(motus, '\\]$')) %>% 
  left_join(gtdb.tax %>% rename(species.gtdb=species), by='motus') %>% 
  filter(abs(value) > cutoff) %>%
  filter(!is.na(genus)) %>% 
  group_by(genus) %>% 
  mutate(n.all=length(unique(species))) %>% 
  summarise(n.all=length(unique(species)), v=var(abs(value)), 
            m=mean(abs(value))) %>% 
  filter(n.all > 2) %>% 
  arrange(desc(v)) %>%
  filter(v>0.025)

genera <- as_tibble(prop.table(motus.gtdb.lvls$genus, 2)[
  genera$genus,colnames(feat.rel)], 
  rownames='genus') %>% 
  pivot_longer(-genus) %>% 
  mutate(value=log10(value+1e-04)) %>% 
  mutate(site=str_remove(name, '[0-9]*$')) %>% 
  mutate(site=site.dict[site]) %>%
  mutate(site=factor(site, levels=names(site.colours))) %>% 
  group_by(genus, site) %>% 
  summarise(med=median(value)) %>% 
  group_by(genus) %>% 
  summarise(z=var(med)) %>% 
  arrange(desc(z)) %>% 
  head(n=7) %>% pull(genus)



df.test <- as_tibble(feat.rel[species.included,], 
          rownames='species') %>% 
  pivot_longer(-species) %>% 
  mutate(motus=str_remove(species, '.*\\[')) %>% 
  mutate(motus=str_remove(motus, '\\]$')) %>% 
  left_join(gtdb.tax %>% rename(species.gtdb=species), by='motus') %>% 
  mutate(site=str_remove(name, '[0-9]{3}$')) %>% 
  filter(genus %in% genera) %>% 
  group_by(genus, species, site) %>% 
  summarise(m=mean(log10(value+1e-04)), species.gtdb=unique(species.gtdb), 
            .groups='drop') %>% 
  mutate(site=factor(site, levels=c('BF', 'GH', 'DM', 'AG', 'SW', 'KY')))
df.test.order <- df.test %>% 
  group_by(species) %>% 
  summarise(m=mean(m), g=unique(genus), s.g=unique(species.gtdb)) %>% 
  arrange(desc(m))
g <- df.test %>% 
  mutate(species=factor(species, df.test.order$species)) %>% 
  ggplot(aes(y=site, x=species, fill=m)) + 
    geom_tile() + 
    facet_grid(~genus, space='free', scales = 'free') + 
    scale_fill_gradientn(colors=RColorBrewer::brewer.pal(9, 'Blues'),
                         name='log10(mean rel. ab.)') + 
    theme_bw() + theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       panel.grid = element_blank()) + 
    xlab('') + ylab('')
g$data %>% 
  write_csv('./files/source_data/Fig1f.csv')
ggsave(g, filename='./figures/prevalence/species_ab.pdf',
       width = 12, height = 4, useDingbats=FALSE)

# specific species

# genus boxplots
g.genus.ab <- as_tibble(prop.table(motus.gtdb.lvls$genus, 2)[
  genera,colnames(feat.rel)], 
                        rownames='genus') %>% 
  pivot_longer(-genus) %>% 
  # filter(genus %in% unique(tax.cluster$genus)) %>% 
  mutate(value=log10(value+1e-04)) %>% 
  mutate(site=str_remove(name, '[0-9]*$')) %>% 
  mutate(site=site.dict[site]) %>%
  mutate(site=factor(site, levels=names(site.colours))) %>% 
  # mutate(genus=factor(genus, levels=genera)) %>% 
  ggplot(aes(x=site, y=value, fill=site)) + 
  geom_boxplot() + 
  facet_grid(~genus) + 
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major.x = element_blank()) + 
  scale_fill_manual(values=site.colours)
g.genus.ab$data %>% 
  write_csv('./files/source_data/Fig1e.csv')
ggsave(g.genus.ab, filename=here('figures','prevalence',
                                 'genus_abundances_newdb.pdf'),
       width = 12, height = 5, useDingbats=FALSE)

# ##############################################################################
# rest is old
# prevalence heatmap with those genera annotated?

genus.cols <- c(   
  'g__Treponema_D'='#F28E2B',
  'g__Dialister'='#EDC948',
  'g__Phocaeicola'='#76B7B2',
  'g__Bacteroides'='#4E79A7',
  'g__Oribacterium'='#FF9DA7',
  'g__Prevotella'='#59A14F',
  'g__Bifidobacterium'='#E15759',
  'g__Cryptobacteroides'='#B07AA1',
  'other'='#D3D3D3'
  )

df.comp <- as_tibble(feat.rel[species.included,], rownames='species') %>% 
  pivot_longer(-species) %>% 
  mutate(value=log10(value+1e-04)) %>% 
  mutate(site=str_remove(name, '[0-9]*$')) %>% 
  group_by(species, site) %>% 
  group_by(species, site) %>% summarize(med=mean(value)) %>% 
  pivot_wider(names_from=site, values_from = med) %>% 
  mutate(motus=str_extract(species, '\\[(ref|meta|ext|NEWDB).+\\]$')) %>% 
  mutate(motus=str_remove(motus, '\\[')) %>% 
  mutate(motus=str_remove(motus, '\\]')) %>% 
  full_join(gtdb.tax %>% select(motus, genus), by='motus') %>% 
  mutate(group=case_when(genus %in% genera~genus, TRUE~'other')) %>% 
  mutate(group=factor(group, levels=rev(names(genus.cols)))) %>%
  arrange(group)
g1 <- df.comp %>%   
  ggplot(aes(x=SW, y=BF, col=as.factor(group))) + 
    geom_abline(slope=1, intercept = 0, lty=2) + 
    geom_point() + 
    theme_bw() + theme(panel.grid.minor = element_blank()) +
    scale_colour_manual(values=genus.cols) +
    ylim(-4, -1.5) + xlim(-4, -1.5) +
    NULL
g2 <-  df.comp %>%   
  ggplot(aes(x=KY, y=BF, col=as.factor(group))) + 
  geom_abline(slope=1, intercept = 0, lty=2) + 
  geom_point() + 
  theme_bw() + theme(panel.grid.minor = element_blank()) +
  scale_colour_manual(values=genus.cols) +
  ylim(-4, -1.5) + xlim(-4, -1.5) +
  NULL
g3 <-  df.comp %>%   
  ggplot(aes(x=AG, y=BF, col=as.factor(group))) + 
  geom_abline(slope=1, intercept = 0, lty=2) + 
  geom_point() + 
  theme_bw() + theme(panel.grid.minor = element_blank()) +
  scale_colour_manual(values=genus.cols) +
  ylim(-4, -1.5) + xlim(-4, -1.5) +
  NULL
g <- cowplot::plot_grid(g1, g2, g3, nrow=1)
g1$data %>% 
  write_csv('./files/source_data/ED6b.csv')
ggsave(g, filename=here('figures','prevalence','abundance_comp_newdb.pdf'), 
       width = 15, height = 5, useDingbats=FALSE)

# prevtolla/bacteroides co-occurrence
df.pbp <- as_tibble(prop.table(motus.gtdb.lvls$genus, 2)[
  c('g__Prevotella', 'g__Bacteroides', 'g__Phocaeicola'),colnames(tbl.bacteria)], 
  rownames='genus') %>% 
  pivot_longer(-genus) %>% 
  mutate(value=log10(value+1e-04)) %>% 
  pivot_wider(names_from = genus, values_from = value) %>% 
  pivot_longer(-c(name, g__Prevotella), names_to = 'type') %>% 
  mutate(site=str_remove(name, '[0-9]*$')) %>% 
  mutate(site=site.dict[site]) %>%
  mutate(site=factor(site, levels=names(site.colours))) 
# prevalences
g.pbp <- df.pbp %>% 
  group_by(site, type) %>% 
  group_map(.f=function(.x, .y){
    tibble(.y, cooccur=sum(.x$g__Prevotella> -4 & .x$value > -4)/nrow(.x))}) %>% 
  bind_rows() %>% 
  mutate(label=sprintf(fmt='%.2f', cooccur)) %>% 
  ggplot(aes(y=site, x=type, fill=cooccur)) + 
    geom_tile() + geom_text(aes(label=label)) +
    theme_bw() + 
    scale_fill_viridis_c(limits=c(0,1), guide='none')
    
# plot
g.pbp.2 <- df.pbp %>% 
  ggplot(aes(x=g__Prevotella, y=value, colour=site)) + 
    geom_point() +
    facet_grid(~type) + 
    theme_bw() + theme(panel.grid.minor = element_blank()) +
    scale_colour_manual(values=site.colours, guide='none') +
    xlab('log10 rel. ab. Prevotella')
g <- cowplot::plot_grid(g.pbp.2, g.pbp, rel_widths = c(0.7, 0.3))  
g.pbp.2$data %>% 
  write_csv('./files/source_data/ED6c.csv')
g.pbp$data %>% 
  write_csv('./files/source_data/ED6d.csv')
ggsave(g, filename=here('figures', 'prevalence', 
                        'bacteroides_prevotella_newdb.pdf'),
       width = 8, height = 4, useDingbats=FALSE)

# bifidos?
# look at Bifido?
g <- as_tibble(prop.table(motus.gtdb.lvls$motus, 2)[species.included,], 
          rownames='species') %>% 
  pivot_longer(-species) %>% 
  filter(str_detect(species, 'Bifido')) %>% 
  mutate(value=log10(value+1e-04)) %>%
  mutate(site=str_remove(name, '[0-9]*$')) %>% 
  mutate(site=site.dict[site]) %>%
  mutate(site=factor(site, levels=names(site.colours))) %>% 
  ggplot(aes(x=site, y=value, fill=site)) + 
    geom_boxplot() +
    facet_wrap(~species) + 
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       panel.grid.major.x = element_blank(),
                       axis.text.x = element_blank(),
                       axis.ticks.x = element_blank()) + 
    scale_fill_manual(values=site.colours) + 
    ylab('log10 rel. ab.') + xlab('')
ggsave(g, filename=here('figures', 'prevalence', 'bifidos_newdb.pdf'),
       width = 8, height = 5, useDingbats=FALSE)


# other bacteria?
g <- as_tibble(prop.table(motus.gtdb.lvls$motus, 2)[species.included,], 
               rownames='species') %>% 
  pivot_longer(-species) %>% 
  filter(str_detect(species, 'Bac')) %>% 
  mutate(value=log10(value+1e-04)) %>%
  mutate(site=str_remove(name, '[0-9]*$')) %>% 
  mutate(site=site.dict[site]) %>%
  mutate(site=factor(site, levels=names(site.colours))) %>% 
  ggplot(aes(x=site, y=value, fill=site)) + 
  geom_boxplot() +
  facet_wrap(~species) + 
  theme_bw() + theme(panel.grid.minor = element_blank(),
                     panel.grid.major.x = element_blank(),
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank()) + 
  scale_fill_manual(values=site.colours) + 
  ylab('log10 rel. ab.') + xlab('')
