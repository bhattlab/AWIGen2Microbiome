# ##############################################################################
#
## Look at correlation of old vs new motus database
#
# ##############################################################################

library("tidyverse")

plot.params <- yaml::yaml.load_file(here::here('files', 'params.yml'))
site.dict <- unlist(plot.params$sites)
site.colours <- unlist(plot.params$sites_colours)

cluster.loc <- paste0('/Volumes/lab_asbhatt/wirbel/projects/AWIgen2/',
                      'motus_extension/profile_ext/')

# motus original
motus.ref <- read.table('./data/classification/motus.tsv',
                        sep='\t',stringsAsFactors = FALSE, check.names = FALSE,
                        row.names = 1, header = TRUE, quote = '', 
                        comment.char = '', skip=2)
motus.ref <- as.matrix(motus.ref)
motus.ref <- motus.ref[,colSums(motus.ref) > 3000]
motus.ref <- motus.ref[rowMeans(motus.ref!=0)>0,]
motus.ref.rel <- prop.table(motus.ref, 2)

# motus new database
motus.new <- read.table('./data/classification/motus_newdb.tsv',
                        sep='\t',stringsAsFactors = FALSE, check.names = FALSE,
                        row.names = 1, header = TRUE, quote = '', 
                        comment.char = '', skip=2)
motus.new <- as.matrix(motus.new)
motus.new <- motus.new[,colnames(motus.ref)]
motus.new <- motus.new[rowMeans(motus.new!=0)>0,]
motus.new.rel <- prop.table(motus.new, 2)

# relative abundance of unassigned calls
df.unassigned <- enframe(motus.ref['unassigned',], name = 'Sample',
                         value='unassigned_ref') %>% 
  full_join(enframe(motus.new['unassigned',], name = 'Sample',
                    value='unassigned_ext'), by='Sample') %>% 
  full_join(enframe(motus.ref.rel['unassigned',], name = 'Sample',
                    value='unassigned_ref_rel'), by='Sample') %>% 
  full_join(enframe(motus.new.rel['unassigned',], name = 'Sample',
                    value='unassigned_ext_rel'), by='Sample') 


g <- df.unassigned %>% 
  pivot_longer(-Sample) %>% 
  mutate(value=case_when(str_detect(name, 'rel')~log10(value + 1e-04), 
                         TRUE~log10(value))) %>% 
  mutate(name=str_remove(name, 'unassigned_')) %>% 
  mutate(rel_ab=case_when(str_detect(name, '_rel')~'Relative abundance',
                          TRUE~'Number of counts')) %>% 
  mutate(name=str_remove(name, '_rel')) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  mutate(Site=str_remove(Sample, '[0-9]*$')) %>% 
  mutate(Site=site.dict[Site]) %>% 
  filter(!is.na(Site)) %>% 
  ggplot(aes(x=ref, y=ext, col=Site)) +
    geom_abline(slope = 1, intercept = 0, colour='#D3D3D3', lty=2) +
    facet_wrap(~rel_ab, scales='free') + 
    geom_point(pch=16) + 
    scale_colour_manual(values=site.colours) + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) + 
    xlab('Original database') + 
    ylab('Extended database')
g.unassigned <- g
g.unassigned$data %>% 
  arrange(rel_ab) %>% 
  write_csv('./files/source_data/ED3c.csv')
ggsave(g, filename='./figures/revision/unassigned_counts.pdf',
       width = 8, height = 4, useDingbats=FALSE)


# cumulative fraction of abundance
rel.ab.new <- colSums(motus.new.rel[str_detect(
  rownames(motus.new.rel), 'NEWDB'),]) %>% enframe
g <- rel.ab.new %>% 
  mutate(site=str_remove(name, '[0-9]+$')) %>% 
  filter(!site %in% c('Buffer', 'Zymo')) %>% 
  mutate(site.nice=site.dict[site]) %>% 
  mutate(site.nice=factor(site.nice, levels=names(site.colours))) %>% 
  ggplot(aes(x=site.nice, y=value, fill=site.nice)) + 
    geom_jitter(width = 0.12, col='#D3D3D3', pch=16) +
    geom_boxplot(outlier.shape = NA) + 
    xlab('') + ylab("Cumulative rel. ab. from added MAGs") + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major.x = element_blank()) + 
    scale_fill_manual(values=site.colours, guide='none')
g.cumsum <- g
g.cumsum$data %>% 
  write_csv('./files/source_data/ED3e.csv')
ggsave(g, filename='./figures/revision/cum_abund_new_motus.pdf',
       width = 5, height = 4, useDingbats=FALSE)

# correlation
df.comp <- motus.ref.rel %>% as_tibble(rownames='motus') %>% 
  pivot_longer(-motus, names_to = 'Sample', values_to = 'value_old') %>% 
  full_join(motus.new.rel %>% as_tibble(rownames='motus') %>% 
              pivot_longer(-motus, names_to = 'Sample', 
                           values_to = 'value_new'),
            by=c('motus', 'Sample')) %>% 
  mutate(value_old=replace_na(value_old, 0)) %>% 
  mutate(value_new=replace_na(value_new, 0)) %>% 
  mutate(f=value_old==0 & value_new==0) %>% 
  filter(!f) %>% 
  mutate(Site=str_remove(Sample, '[0-9]*$')) %>% 
  mutate(Site=site.dict[Site]) %>% 
  filter(!is.na(Site))

g <- df.comp %>% 
  group_by(Sample, Site) %>% 
  summarise(m=cor(value_old, value_new, method='spearman'), 
            .groups='drop') %>% 
  mutate(Site=factor(Site, levels=names(site.colours))) %>% 
  ggplot(aes(x=Site, y=m, fill=Site)) + 
  geom_jitter(width = 0.12, col='#D3D3D3', pch=16) +
  geom_boxplot(outlier.shape = NA) + 
  xlab('') + ylab("Spearman's rho between original and extended database") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank()) + 
  scale_fill_manual(values=site.colours, guide='none')
g.cor <- g
g.cor$data %>% 
  write_csv('./files/source_data/ED3d.csv')
ggsave(g, filename='./figures/revision/correlation_new_motus.pdf',
       width = 5, height = 4, useDingbats=FALSE)

# ##############################################################################
# convert to GTDB and save
motus.tax.gtdb <- read_tsv(paste0(cluster.loc, '../assign_MAGs_v303/',
                                  'comparison/extension/NEWDB/db_mOTU/',
                                  'GTDB_tax.tsv'),
                           col_names=c('mOTU', 'domain', 'phylum', 'class', 
                                       'order', 'family', 'genus', 'species'),
                           col_types = cols())
.f_get_motus_level <- function(lvl='phylum', motus.tbl){
  message(lvl)
  stopifnot(lvl %in% colnames(motus.tax.gtdb)[2:8])
  rownames(motus.tbl) <- str_extract(
    rownames(motus.tbl), '((ref|meta|ext)_mOTU_v[0-9]*_[0-9]{5}|(unassigned))')
  new.level <- unique(motus.tax.gtdb[[lvl]])
  new.mat <- matrix(0, nrow = length(new.level)+1, ncol=ncol(motus.tbl),
                    dimnames = list(c(new.level, 'unassigned'), 
                                    colnames(motus.tbl)))
  for (x in new.level){
    incl.motus <- motus.tax.gtdb %>% 
      filter(eval(sym(lvl)) ==x) %>% 
      filter(mOTU %in% rownames(motus.tbl)) %>% 
      pull(mOTU)
    if (length(incl.motus) > 0){
      new.mat[x,] <- colSums(motus.tbl[incl.motus,,drop=FALSE])
    }
  }
  new.mat['unassigned',] <- motus.tbl['unassigned',]
  new.mat <- new.mat[rowMeans(new.mat!=0) > 0,]
  return(new.mat)
}

lvls <- c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')
motus.gtdb.lvls <- map(lvls[-1], .f = .f_get_motus_level, motus.tbl=motus.new)
names(motus.gtdb.lvls) <- lvls[-1]
motus.gtdb.lvls[['motus']] <- motus.new

save(motus.gtdb.lvls, file='./data/classification/new_db_tables.RData')

# ##############################################################################
# save it nicely with cowplot?

g <- cowplot::plot_grid(g.unassigned,
                   cowplot::plot_grid(g.cor, g.cumsum, nrow=1, align = 'vh', 
                                      axis='tblr'), nrow=2)
ggsave(g, filename='./figures/revision/new_db_figure.pdf',
       width = 8, height = 6, useDingbats=FALSE)


# ##############################################################################
# Compare this new database to other datasets and make some phylum plots, maybe?

df.other <- list()
df.other.gtdb <- list()
df.other.phylum <- list()
for (x in c('PRJNA319574_Schirmer_Cell_2016', 'PRJDB4176_Yachida_NatMed_2019', 
            'PRJNA400072_Franzosa_NatMed_2018')){
  dataset <- str_split(x, pattern='_')[[1]][2]
  motus.ext <- read.table(paste0(
    '/Volumes/lab_asbhatt/data/public_data/', x,
    '/classification/motus_all_v3.0.3_AWIGEN.tsv'),
    sep='\t',stringsAsFactors = FALSE, check.names = FALSE,
    row.names = 1, header = TRUE, quote = '', 
    comment.char = '', skip=2)
  motus.ext <- motus.ext[rowMeans(motus.ext!=0)>0,]
  motus.ext.phylum <- .f_get_motus_level(lvl='phylum', as.matrix(motus.ext))
  motus.ext.gtdb <- .f_get_motus_level(lvl='species', as.matrix(motus.ext))
  df.other[[dataset]] <- motus.ext
  df.other.gtdb[[dataset]] <- motus.ext.gtdb
  df.other.phylum[[dataset]] <- motus.ext.phylum
}

# combine everything into a single dataset
all.species <- union(rownames(motus.gtdb.lvls$species),
                     map(df.other.gtdb, rownames) %>% unlist) 
                     
mat.all <- matrix(0, nrow=length(all.species),
                  ncol=ncol(motus.gtdb.lvls$species) + sum(
                    map(df.other.gtdb, ncol) %>% unlist))
rownames(mat.all) <- all.species
colnames(mat.all) <- c(colnames(motus.gtdb.lvls$species), 
                       map(df.other.gtdb, colnames) %>% unlist())
mat.all[rownames(motus.gtdb.lvls$species),colnames(motus.gtdb.lvls$species)] <- 
  motus.gtdb.lvls$species
for (x in seq_along(df.other)){
  tmp <- df.other.gtdb[[x]]
  mat.all[rownames(tmp),colnames(tmp)] <- tmp 
}

# beta-diversity
hist(log10(colSums(mat.all)), 100)
pco.res <- labdsv::pco(vegan::vegdist(vegan::rrarefy(t(mat.all), 3000)))
df.pco <- as_tibble(as.data.frame(pco.res$points), rownames='Sample') %>% 
  mutate(
    study=case_when(Sample %in% colnames(motus.gtdb.lvls$species)~'AWI-Gen 2',
                    Sample %in% colnames(df.other$Schirmer)~'Schirmer',
                    Sample %in% colnames(df.other$Franzosa)~'Franzosa',
                    Sample %in% colnames(df.other$Yachida)~'Yachida'))
labels <- head(pco.res$eig/sum(pco.res$eig[pco.res$eig>0]), n=2)
g <- df.pco %>% 
  ggplot(aes(x=V1, y=V2, col=study)) + 
    geom_point() + 
    theme_bw() + theme(panel.grid = element_blank()) + 
    xlab(paste0('PCo 1 [', sprintf(fmt='%.2f', labels[1]*100), '%]')) + 
    ylab(paste0('PCo 2 [', sprintf(fmt='%.2f', labels[2]*100), '%]')) + 
    scale_colour_manual(values=c('#59A14F', '#4E79A7', '#F28E2B', '#76B7B2'),
                        name='Study')
g$data %>% 
  write_csv('./files/source_data/ED4b.csv')
ggsave(g, filename='./figures/revision/pcoa_external.pdf',
       width = 5, height = 4, useDingbats=FALSE)


# phylum level boxplots
all.phyla <- union(rownames(motus.gtdb.lvls$phylum),
                     map(df.other.phylum, rownames) %>% unlist) 
mat.phylum <- matrix(0, nrow=length(all.phyla),
                  ncol=ncol(motus.gtdb.lvls$phylum) + sum(
                    map(df.other.phylum, ncol) %>% unlist))
rownames(mat.phylum) <- all.phyla
colnames(mat.phylum) <- c(colnames(motus.gtdb.lvls$phylum), 
                       map(df.other.phylum, colnames) %>% unlist())
mat.phylum[rownames(motus.gtdb.lvls$phylum),colnames(motus.gtdb.lvls$phylum)] <- 
  motus.gtdb.lvls$phylum
for (x in seq_along(df.other.phylum)){
  tmp <- df.other.phylum[[x]]
  mat.phylum[rownames(tmp),colnames(tmp)] <- tmp 
}

mat.phylum.rel <- prop.table(mat.phylum, 2)
mean.phylum <- rowMeans(mat.phylum.rel)
g <- mat.phylum.rel %>% as_tibble(rownames='phylum') %>% 
  pivot_longer(-phylum) %>% 
  group_by(phylum) %>% 
  mutate(mean=mean(value)) %>% 
  ungroup() %>% 
  filter(phylum!='unassigned') %>% 
  filter(mean > 0.005) %>% 
  arrange(desc(mean)) %>% 
  mutate(phylum=factor(phylum, levels=unique(phylum))) %>% 
  mutate(study=case_when(
    name %in% colnames(motus.gtdb.lvls$species)~'AWI-Gen 2',
    name %in% colnames(df.other$Schirmer)~'Schirmer',
    name %in% colnames(df.other$Franzosa)~'Franzosa',
    name %in% colnames(df.other$Yachida)~'Yachida')) %>% 
  ggplot(aes(x=phylum, y=value, fill=study)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1),
                pch=16, colour='grey90') + 
    geom_boxplot(outlier.shape = NA, colour='black') + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major.x = element_blank()) +
    xlab('') + ylab('Relative abundance') + 
    scale_fill_manual(values=c('#59A14F', '#4E79A7', '#F28E2B', '#76B7B2'),
                        name='Study')
g$data %>% 
  write_csv('./files/source_data/ED4c.csv')
ggsave(g, filename='./figures/revision/pyla_external.pdf',
       width = 7, height = 4, useDingbats=FALSE)


mat.phylum.rel %>% as_tibble(rownames='phylum') %>% 
  pivot_longer(-phylum) %>% 
  group_by(phylum) %>% 
  mutate(mean=mean(value)) %>% 
  ungroup() %>% 
  filter(phylum!='unassigned') %>% 
  filter(mean > 0.005) %>% 
  mutate(study=case_when(
    name %in% colnames(motus.gtdb.lvls$species)~'AWI-Gen 2',
    name %in% colnames(df.other$Schirmer)~'Schirmer',
    name %in% colnames(df.other$Franzosa)~'Franzosa',
    name %in% colnames(df.other$Yachida)~'Yachida')) %>% 
  group_by(phylum) %>% 
  group_map(.f=function(.x, .y){
    t <- kruskal.test(value~study, data=.x)
    tibble(p.val=t$p.val) %>% 
      mutate(.y)
  }) %>% bind_rows()

#' # A tibble: 9 × 2
#' p.val phylum              
#' <dbl> <chr>               
#' 1 2.51e-156 p__Actinobacteriota 
#' 2 6.67e-185 p__Bacteroidota     
#' 3 2.33e-276 p__Cyanobacteria    
#' 4 5.32e- 42 p__Firmicutes       
#' 5 1.44e-165 p__Firmicutes_A     
#' 6 3.29e-107 p__Firmicutes_C     
#' 7 2.69e-206 p__Proteobacteria   
#' 8 3.06e- 82 p__Spirochaetota    
#' 9 6.46e-137 p__Verrucomicrobiota

new.db <- enframe(
  colSums(motus.new.rel[str_detect(rownames(motus.new.rel), 'NEWDB'),]),
  name='Sample', value='newdb') %>% 
  mutate(study='AWI-Gen 2') %>% 
  bind_rows(map(names(df.other), .f = function(x){
    enframe(
      colSums(prop.table(df.other[[x]])[str_detect(rownames(df.other[[x]]), 'NEWDB'),]),
      name='Sample', value='newdb') %>% 
      mutate(study=x)
  }) %>% bind_rows())

g <- new.db %>% 
  ggplot(aes(x=study, y=newdb, fill=study)) + 
  geom_jitter(width = 0.1, pch=16, colour='grey90') + 
  geom_boxplot(outlier.shape = NA, colour='black') + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_blank()) +
  xlab('') + ylab('Relative abundance of added MAGs') + 
  scale_fill_manual(values=c('#59A14F', '#4E79A7', '#F28E2B', '#76B7B2'),
                    name='Study')
ggsave(g, filename='./figures/revision/new_db_external.pdf',
       width = 4, height = 3, useDingbats=FALSE)
