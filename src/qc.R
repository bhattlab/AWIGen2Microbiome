# ##############################################################################
#
## Initial QC (Buffers and empty controls)
##    also some data cleaning
#
# ##############################################################################

library('here')
library("tidyverse")
library("yaml")

plot.params <- yaml.load_file(here('files', 'params.yml'))
site.dict <- unlist(plot.params$sites)
site.colours <- unlist(plot.params$sites_colours)

# check that the colours can be distinguished
tibble(x=seq_along(site.colours), y=1, site=names(site.colours)) %>% 
  ggplot(aes(x=x, y=y, fill=site)) +
    geom_tile() + 
    theme_void() + 
    scale_fill_manual(values=site.colours)

# ##############################################################################
# Function to calculate significant differences between sites

.f_cohens_d <- function(x, y){
  d <- mean(x) - mean(y)
  
  n <- sqrt(((length(x)-1)*var(x) + 
               (length(y)-1)*var(y))/
              (length(x) + length(y) - 2))
  return(d/n)
}

.f_site_assoc <- function(tbl, col.name){
  stopifnot(is.character(col.name))
  stopifnot(col.name %in% colnames(tbl))
  if (!'site.city' %in% colnames(tbl)){
    stopifnot('Sample' %in% colnames(tbl))
    tbl <- tbl %>% mutate(site=str_remove(Sample, '[0-9]*$')) %>% 
      mutate(site.city=site.dict[site]) %>% 
      mutate(site.city=factor(site.city, levels = names(site.colours)))
  }
  f <- paste0(col.name, '~site.city')
  fit <- lm(formula = as.formula(f),
            data=tbl %>% 
              filter(!str_detect(site, '(Buffer|Zymo)')))
  res <- anova(fit)
  p.overall <- res['site.city','Pr(>F)']
  
  # site-vs-site testing
  pairwise.site.list <- t(combn(names(site.colours), 2))
  ef <- c()
  p <- c()
  for (i in seq_len(nrow(pairwise.site.list))){
    fit <- lm(formula = as.formula(f),
              data=tbl %>% 
                filter(site.city %in% pairwise.site.list[i,]))
    
    res <- anova(fit)
    p <- c(p, res['site.city', 'Pr(>F)'])
    ef <- c(ef, .f_cohens_d(tbl %>% 
                        filter(site.city == pairwise.site.list[i,1]) %>% 
                        pull(col.name),
                      tbl %>% 
                        filter(site.city == pairwise.site.list[i,2]) %>% 
                        pull(col.name)))
  }
  pairwise.site.list <- pairwise.site.list %>% 
    as.data.frame() %>% 
    as_tibble() %>% 
    mutate(cohens_d=ef, p.val=p)
  g.heat <- pairwise.site.list %>% 
    mutate(V1=factor(V1, levels = unique(V1)), 
           V2=factor(V2, levels = unique(V2))) %>% 
    mutate(annot=case_when(p.val < 0.0001~'***',
                           p.val < 0.001~'**',
                           p.val < 0.05~'*',
                           TRUE~'')) %>% 
    ggplot(aes(x=V1, y=V2, fill=cohens_d)) + 
      geom_tile() + 
      theme_minimal() + theme(panel.grid = element_blank()) + 
      xlab('') + ylab('') + 
      scale_fill_gradient2(name="Cohen's D") +
      geom_text(aes(label=annot))
  g.plot <- tbl %>% 
    filter(!str_detect(site, '(Buffer|Zymo)')) %>% 
    mutate(y=!!sym(col.name)) %>% 
    ggplot(aes(x=site.city, y=y)) + 
      geom_jitter(colour='#D3D3D3', width = 0.15) +
      geom_boxplot(aes(fill=site.city), outlier.shape = NA, alpha=0.7) + 
      scale_fill_manual(values=site.colours, guide='none') + 
      theme_bw() + theme(panel.grid.minor = element_blank(), 
                         panel.grid.major.x = element_blank()) + 
      xlab('') + ylab(col.name)
  return(list("p.value.overall"=p.overall, "pairwise.tests"=pairwise.site.list,
              "plot.ef"=g.heat, "plot"=g.plot))
  
}

# ##############################################################################
# Number of reads
read.counts <- read_tsv(here('data', 'misc', 'read_counts.tsv'), 
                        col_types = cols()) %>% 
  mutate(site=str_remove(Sample, '[0-9]*$')) %>% 
  mutate(site.city=site.dict[site]) %>% 
  mutate(site.city=factor(site.city, levels = names(site.colours)))

# any samples with a weird number of reads?
read.counts %>% 
  ggplot(aes(x=hostremoved_reads)) + 
    geom_histogram(bins=100) + 
    scale_x_log10() + xlab('No. of reads after host-removal')

# compare to Awigen1
read.counts.a1 <- read_tsv(here('data', 'misc', 'read_counts_awigen1.tsv'),
                           col_types = cols()) 

g1 <- read.counts %>% 
  mutate(study=case_when(str_detect(Sample, 'Buffer|Zymo')~'internal_controls',
                         TRUE~'awigen2')) %>% 
  bind_rows(read.counts.a1 %>% mutate(study='awigen1')) %>% 
  select(raw_reads, dedup_reads, trimmed_reads, hostremoved_reads, study) %>% 
  pivot_longer(-study) %>% 
  mutate(name=factor(name, levels = c('raw_reads', 'dedup_reads', 
                                      'trimmed_reads', 
                                      'hostremoved_reads'))) %>% 
  ggplot(aes(x=name, y=value, fill=study)) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1), 
                col='#D3D3D3') +
    geom_boxplot(outlier.shape = NA, alpha=0.7) + 
    xlab('') + ylab('Number of reads (log10)') + 
    theme_bw() + theme(panel.grid.minor = element_blank(), 
                       panel.grid.major.x = element_blank()) + 
    scale_fill_manual(values=unlist(plot.params$stage_colours), name='Study') + 
    scale_y_log10(limits=c(1e06,1.2e08)) + 
    annotation_logticks(sides='l')
ggsave(g1, filename = here('figures', 'general', 'read_counts.pdf'),
       width = 6, height = 4, useDingbats=FALSE)

# have a second look at samples with fewer than 1e07 hostremoved_reads?
hostremoved_reads_cutoff <- 1e07

g2 <- read.counts %>% 
  mutate(study=case_when(str_detect(Sample, 'Buffer|Zymo')~'internal_controls',
                         TRUE~'awigen2')) %>% 
  bind_rows(read.counts.a1 %>% mutate(study='awigen1')) %>% 
  select(dedup_frac, trimmed_frac, hostremoved_frac, study) %>% 
  pivot_longer(-study) %>% 
  mutate(name=factor(name, levels = c('dedup_frac', 'trimmed_frac', 
                                      'hostremoved_frac'))) %>% 
  ggplot(aes(x=name, y=value, fill=study)) + 
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1), 
              col='#D3D3D3') +
  geom_boxplot(outlier.shape = NA, alpha=0.7) + 
  xlab('') + ylab('Fraction of reads') + 
  theme_bw() + theme(panel.grid.minor = element_blank(), 
                     panel.grid.major.x = element_blank()) + 
  scale_fill_manual(values=unlist(plot.params$stage_colours), name='Study') 
ggsave(g2, filename = here('figures', 'general', 'read_count_fraction.pdf'),
       width = 6, height = 4, useDingbats=FALSE)

# have a second look at samples with low fraction of reads after host removal
hostremoved_frac_cutoff <- 0.55

# number of reads per site
test.hostremoved_reads <- .f_site_assoc(read.counts, 'hostremoved_reads')
ggsave(test.hostremoved_reads$plot + ylab('Number of reads after host removal'), 
       filename=here('figures', 'general', 'hostremoved_reads_per_site.pdf'),
       width = 5, height = 4, useDingbats=FALSE)
ggsave(test.hostremoved_reads$plot.ef + 
         theme(axis.text.x = element_text(angle=45, hjust=1)),
       filename=here('figures', 'general', 
                     'hostremoved_reads_per_site_ef.pdf'),
       width = 4, height = 3, useDingbats=FALSE)

# human read-counts across sites?
test.host_frac <- .f_site_assoc(
  read.counts %>% 
    mutate(host_frac = hostremoved_reads/trimmed_reads), 'host_frac')
ggsave(test.host_frac$plot + ylab('Fraction of reads after host removal'), 
       filename=here('figures', 'general', 'hostremoved_frac_per_site.pdf'),
       width = 5, height = 4, useDingbats=FALSE)
ggsave(test.host_frac$plot.ef + 
         theme(axis.text.x = element_text(angle=45, hjust=1)),
       filename=here('figures', 'general', 
                     'hostremoved_frac_per_site_ef.pdf'),
       width = 4, height = 3, useDingbats=FALSE)
# have a second look at samples with lots of human reads after trimming?
host_frac_cutoff <- 0.25
# is this redundant with the hostremoved_frac?

read.counts %>% 
  mutate(host_frac = 1 - hostremoved_reads/trimmed_reads) %>% 
  ggplot(aes(x=hostremoved_frac, y=host_frac)) + 
    geom_point() + 
    geom_hline(yintercept = host_frac_cutoff) + 
    geom_vline(xintercept = hostremoved_frac_cutoff) + 
    theme_bw() + theme(panel.grid.minor = element_blank())
# more or less, that's okay

# ##############################################################################
# read in classification data 
# motus

motus.all <- read.table(here('data', 'classification', 'motus.tsv'),
                        sep='\t',stringsAsFactors = FALSE, check.names = FALSE,
                        row.names = 1, header = TRUE, quote = '', 
                        comment.char = '', skip=2)
motus.all <- as.matrix(motus.all)

# remove samples with too few motus counts
enframe(colSums(motus.all)) %>% 
  ggplot(aes(x=value)) + 
    geom_histogram(bins=100)
enframe(colSums(motus.all), name="Sample") %>% 
  full_join(read.counts, by='Sample') %>% 
  ggplot(aes(x=value, y=hostremoved_reads)) + 
    geom_point() + 
    geom_vline(xintercept = 3000)

# everything under 3000 counts is pretty weird
motus.all <- motus.all[,colSums(motus.all) > 3000]

# put the buffer samples somewhere else
motus.buffer <- motus.all[,str_detect(colnames(motus.all), 'Buffer|Zymo')]
motus.buffer <- motus.buffer[rowMeans(motus.buffer!=0)>0,]
motus.all <- motus.all[,setdiff(colnames(motus.all), colnames(motus.buffer))]
motus.all <- motus.all[rowMeans(motus.all!=0) > 0,]
motus.rel <- prop.table(motus.all, 2)

# type of motus per sample
g <- motus.rel %>% 
  as_tibble(rownames='species') %>% 
  mutate(type=case_when(str_detect(species, 'ref_mOTU_v3')~'ref',
                        str_detect(species, 'meta_mOTU_v3')~'meta',
                        str_detect(species, 'ext_mOTU_v3')~'ext',
                        TRUE~'unassigned')) %>% 
  pivot_longer(-c(species, type)) %>% 
  group_by(name, type) %>% 
  summarise(ab=sum(value)) %>% 
  mutate(site=str_remove(name, '[0-9]{3}$')) %>% 
  mutate(site=site.dict[site]) %>% 
  mutate(site=factor(site, levels=names(site.colours))) %>%
  mutate(type=factor(type, levels=c('ref', 'meta', 'ext', 'unassigned'))) %>% 
  ggplot(aes(x=site, y=ab, fill=type)) + 
    geom_boxplot(colour='black') +
    xlab('') + ylab('Relative abundance') + 
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       panel.grid.major.x = element_blank()) + 
    scale_fill_manual(values=c('#09425A80', '#009B7680', '#E04F39', '#7F7776'), 
                      name='Type', labels=c('reference', 'metagenomic', 
                                            'environmental', 'unassigned')) + 
    ylim(c(0,1))
ggsave(g, filename=here('figures', 'general', 'motus_type.pdf'),
       width = 8, height = 5, useDingbats=FALSE)

.f_get_motus_level <- function(lvl='phylum', motus.tbl){
  motus.tax <- read_tsv(here('files', 'motus_taxonomy.tsv'), 
                        col_types = cols())
  stopifnot(lvl %in% colnames(motus.tax)[2:7])
  rownames(motus.tbl) <- str_extract(
    rownames(motus.tbl), '((ref|meta|ext)_mOTU_v3_[0-9]{5}|(unassigned))')
  new.level <- unique(motus.tax[[lvl]])
  new.mat <- matrix(0, nrow = length(new.level)+1, ncol=ncol(motus.tbl),
                    dimnames = list(c(new.level, 'unassigned'), 
                                    colnames(motus.tbl)))
  for (x in new.level){
    incl.motus <- motus.tax %>% 
      filter(eval(sym(lvl)) ==x) %>% 
      filter(motus_ID %in% rownames(motus.tbl)) %>% 
      pull(motus_ID)
    if (length(incl.motus) > 0){
      new.mat[x,] <- colSums(motus.tbl[incl.motus,,drop=FALSE])
    }
  }
  new.mat['unassigned',] <- motus.tbl['unassigned',]
  new.mat <- new.mat[rowMeans(new.mat!=0) > 0,]
  return(new.mat)
}

lvls <- c('kingdom', 'phylum', 'class', 'order', 'family', 'genus')
motus.other.levels <- map(lvls, .f = .f_get_motus_level, motus.tbl=motus.all)
names(motus.other.levels) <- lvls
lst.motus <- motus.other.levels
lst.motus[['species']] <- motus.all

# ##############################################################################
# also for GTDB taxonomy
motus.tax.gtdb <- read_tsv(here('files', 'mOTUs_3.0.0_GTDB_tax.tsv'), 
                           col_names=c('mOTU', 'domain', 'phylum', 'class', 
                                       'order', 'family', 'genus', 'species'))
.f_get_motus_level <- function(lvl='phylum', motus.tbl){
  stopifnot(lvl %in% colnames(motus.tax.gtdb)[2:8])
  rownames(motus.tbl) <- str_extract(
    rownames(motus.tbl), '((ref|meta|ext)_mOTU_v3_[0-9]{5}|(unassigned))')
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
motus.gtdb.lvls <- map(lvls, .f = .f_get_motus_level, motus.tbl=motus.all)
names(motus.gtdb.lvls) <- lvls
lst.motus.gtdb <- motus.gtdb.lvls
names(lst.motus.gtdb)[1] <- 'kingdom'

# full tax name for Luicer
motus.tmp <- motus.all
rownames(motus.tmp) <- str_extract(
  rownames(motus.tmp), '((ref|meta|ext)_mOTU_v3_[0-9]{5}|(unassigned))')
motus.tax.gtdb <- motus.tax.gtdb %>% 
  unite(col='new_name', domain, phylum, class, order, family, genus, species, 
        sep=';') 
rownames(motus.tmp) <- motus.tax.gtdb$new_name[match(rownames(motus.tmp), 
                                                     motus.tax.gtdb$mOTU)]
rownames(motus.tmp)[length(rownames(motus.tmp))] <- 'unassigned'
motus.tmp <- as_tibble(motus.tmp, rownames='species') %>% 
  pivot_longer(-species) %>% 
  group_by(species, name) %>%
  summarize(value=sum(value), .groups = 'drop') %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  as.data.frame()
rownames(motus.tmp) <- motus.tmp$species
motus.tmp$species <- NULL
lst.motus.gtdb[['species_full_name']] <- motus.tmp

# ##############################################################################
# read in classification data 
# mpa4

mpa4.all <- read.table(here('data', 'classification', 'mpa.tsv'), 
                       sep='\t', stringsAsFactors = FALSE, check.names = FALSE, 
                       row.names = 1, header = TRUE)
mpa4.all <- as.matrix(mpa4.all)
mpa4.species <- mpa4.all[grep(rownames(mpa4.all), pattern='s__', value=TRUE) %>% 
                           grep(pattern='t__', value=TRUE, invert = TRUE),]
rownames(mpa4.species) <- str_remove(rownames(mpa4.species), '^.*\\|s__')
mpa4.species <- t(t(mpa4.species)/100)

mpa4.buffer <- mpa4.species[,str_detect(colnames(mpa4.all), 'Buffer|Zymo')]
mpa4.all <- mpa4.all[,colnames(motus.all)]
mpa4.species <- mpa4.species[,colnames(motus.all)]

# other levels
mpa.other.levels <- map(lvls, .f=function(lvl){
  remove <- paste0(substr(c(lvls, 'species')[which(lvls == lvl)+1], 1,1), '__')
  keep <- paste0(substr(lvl, 1, 1), '__')
  mpa4.red <- mpa4.all[grep(rownames(mpa4.all), pattern=keep, value=TRUE) %>% 
                         grep(pattern=remove, value=TRUE, invert = TRUE),]
  rownames(mpa4.red) <- str_remove(rownames(mpa4.red), paste0('^.*\\|', keep))
  mpa4.red <- t(t(mpa4.red)/100)
})
names(mpa.other.levels) <- lvls
lst.mpa <- mpa.other.levels
lst.mpa[['species']] <- mpa4.species

# ##############################################################################
# read in classification data 
# phanta (new database)

phanta.tax <- read.table(here('data', 'classification', 'phanta', 
                              'relative_taxonomic_abundance.txt'),
                         sep='\t', header = TRUE)
phanta.rel <- read.table(here('data', 'classification', 'phanta', 
                              'relative_read_abundance.txt'),
                         sep='\t', header = TRUE)

# table wrangling
phanta.tax$Taxon_Lineage_with_IDs <- NULL
rownames(phanta.tax) <- phanta.tax$Taxon_Lineage_with_Names
phanta.tax$Taxon_Lineage_with_Names <- NULL

# split into separate tables for Eukaryotes/Bacteria/Viruses
phanta.buffer <- phanta.tax[,str_detect(colnames(phanta.tax), 'Buffer|Zymo')]
phanta.species <- phanta.tax[,setdiff(colnames(phanta.tax), 
                                      colnames(phanta.buffer))]
phanta.vir <- phanta.species[str_detect(rownames(phanta.species), 
                                        'superkingdom_Viruses'),]
phanta.mic <- phanta.species[str_detect(rownames(phanta.species), 
                                        'superkingdom_Viruses',
                                        negate=TRUE),]

# relative taxonomic abundance for non-virus taxa
phanta.mic <- prop.table(as.matrix(phanta.mic), 2)
phanta.vir <- prop.table(as.matrix(phanta.vir), 2)

# other levels
lvls <- c('superkingdom', lvls[-1])
phanta.other.levels <- map(lvls, .f=function(lvl){
  keep <- paste0(lvl, '_')
  groups <- c()
  for (x in rownames(phanta.mic)){
    tmp <- str_split(x, paste0('\\|', lvl, '_'))[[1]][2]
    tmp <- str_split(tmp, '\\|')[[1]][1]
    groups <- c(groups, tmp)
  }
  stopifnot(length(groups) == nrow(phanta.mic))
  new.mat <- matrix(0, nrow=length(unique(groups)), ncol=ncol(phanta.mic), 
                    dimnames = list(unique(groups), colnames(phanta.mic)))
  for (g in unique(groups)){
    new.mat[g,] <- colSums(phanta.mic[which(groups == g), ,drop=FALSE])
  }
  return(new.mat)})
names(phanta.other.levels) <- lvls
lst.phanta.tax <- phanta.other.levels
rownames(phanta.mic) <- str_remove(rownames(phanta.mic), '^.*\\|species_')
lst.phanta.tax[['species']] <- phanta.mic

# relative taxonomic abundance for viral taxa
phanta.vir.other.levels <- map(lvls, .f=function(lvl){
  keep <- paste0(lvl, '_')
  groups <- c()
  for (x in rownames(phanta.vir)){
    tmp <- str_split(x, paste0('\\|', lvl, '_'))[[1]][2]
    tmp <- str_split(tmp, '\\|')[[1]][1]
    if (is.na(tmp)){
      groups <- c(groups, paste0('no_', lvl))
    } else {
      groups <- c(groups, tmp)
    }
  }
  stopifnot(length(groups) == nrow(phanta.vir))
  new.mat <- matrix(0, nrow=length(unique(groups)), ncol=ncol(phanta.vir), 
                    dimnames = list(unique(groups), colnames(phanta.vir)))
  for (g in unique(groups)){
    new.mat[g,] <- colSums(phanta.vir[which(groups == g), ,drop=FALSE])
  }
  return(new.mat)})
names(phanta.vir.other.levels) <- lvls
phanta.vir.other.levels$superkingdom <- NULL
rownames(phanta.vir) <- str_remove(rownames(phanta.vir), '^.*\\|species_')
phanta.vir.other.levels[['species']] <- phanta.vir
lst.phanta.vir <- phanta.vir.other.levels

# relative read abundance for everything (?), 
#   normalized by assigned/non-assigned reads

# get number of reads file
reads.phanta <- read_tsv(here('data', 'classification', 
                              'phanta', 'total_reads.tsv'))
phanta.rel$Taxon_Lineage_with_IDs <- NULL
rownames(phanta.rel) <- phanta.rel$Taxon_Lineage_with_Names
phanta.rel$Taxon_Lineage_with_Names <- NULL

phanta.rel.species <- phanta.rel[grep(rownames(phanta.rel), 
                                      pattern='\\|species_',
                                      value=TRUE),colnames(phanta.mic)]
phanta.rel.species <- as.matrix(phanta.rel.species)
# turn into relative classified abundance
for (i in colnames(phanta.rel.species)){
  frac.assigned <- reads.phanta %>% 
    filter(Samp_Name==i) %>% pull(Assigned_Step_Three)/reads.phanta %>% 
    filter(Samp_Name==i) %>% pull(Tot_Samp_Reads)
  tmp <- phanta.rel.species[,i]
  tmp <- tmp*frac.assigned
  phanta.rel.species[,i] <- tmp
}
new.line <- 1-colSums(phanta.rel.species)
phanta.rel.species <- rbind(phanta.rel.species, new.line)
rownames(phanta.rel.species)[nrow(phanta.rel.species)] <- 'unassigned'

# fix other levels!
phanta.rel.other.levels <- map(lvls, .f=function(lvl){
  keep <- paste0(lvl, '_')
  groups <- c()
  for (x in rownames(phanta.rel.species)){
    tmp <- str_split(x, paste0('\\|', lvl, '_'))[[1]][2]
    tmp <- str_split(tmp, '\\|')[[1]][1]
    if (is.na(tmp)){
      groups <- c(groups, paste0('no_', lvl))
    } else {
      groups <- c(groups, tmp)
    }
  }
  stopifnot(length(groups) == nrow(phanta.rel.species))
  new.mat <- matrix(0, nrow=length(unique(groups)),
                    ncol=ncol(phanta.rel.species), 
                    dimnames = list(unique(groups), 
                                    colnames(phanta.rel.species)))
  for (g in unique(groups)){
    new.mat[g,] <- colSums(phanta.rel.species[which(groups == g), ,drop=FALSE])
  }
  return(new.mat)})
names(phanta.rel.other.levels) <- lvls
rownames(phanta.rel.species) <- str_remove(rownames(phanta.rel.species), 
                                           '^.*\\|species_')
phanta.rel.other.levels[['species']] <- phanta.rel.species
lst.phanta.rel <- phanta.rel.other.levels


# ##############################################################################
# save classification tables
save(lst.motus, lst.motus.gtdb, lst.mpa, lst.phanta.vir, 
     lst.phanta.tax, lst.phanta.rel,
     file = here('data', 'classification', 'all_classification_tables.RData'))

# ##############################################################################
# test the unclassified fraction for each

test.unclassified.motus <- .f_site_assoc(
  enframe(motus.rel['unassigned',], 
          name='Sample', value='unclassified'), 'unclassified')
ggsave(test.unclassified.motus$plot + 
         ylab('Fraction of unclassified abundance'), 
       filename=here('figures', 'general', 'unclassified_motus.pdf'),
       width = 5, height = 4, useDingbats=FALSE)
ggsave(test.unclassified.motus$plot.ef + 
         theme(axis.text.x = element_text(angle=45, hjust=1)),
       filename=here('figures', 'general', 
                     'unclassified_motus_ef.pdf'),
       width = 4, height = 3, useDingbats=FALSE)

test.unclassified.mpa4 <- .f_site_assoc(
  enframe(mpa4.all['UNCLASSIFIED',]/100, 
          name='Sample', value='unclassified'), 'unclassified')
ggsave(test.unclassified.mpa4$plot + 
         ylab('Fraction of unclassified abundance'), 
       filename=here('figures', 'general', 'unclassified_mpa4.pdf'),
       width = 5, height = 4, useDingbats=FALSE)
ggsave(test.unclassified.mpa4$plot.ef + 
         theme(axis.text.x = element_text(angle=45, hjust=1)),
       filename=here('figures', 'general', 
                     'unclassified_mpa4_ef.pdf'),
       width = 4, height = 3, useDingbats=FALSE)
# phanta
test.unclassified.phanta <- .f_site_assoc(
  enframe(phanta.rel.species['unassigned',], 
          name='Sample', value='unclassified'),'unclassified')
ggsave(test.unclassified.phanta$plot + 
         ylab('Fraction of unclassified abundance'), 
       filename=here('figures', 'general', 'unclassified_phanta.pdf'),
       width = 5, height = 4, useDingbats=FALSE)
ggsave(test.unclassified.phanta$plot.ef + 
         theme(axis.text.x = element_text(angle=45, hjust=1)),
       filename=here('figures', 'general', 
                     'unclassified_phanta_ef.pdf'),
       width = 4, height = 3, useDingbats=FALSE)

# ##############################################################################
# motus is very different here: could it be from the extra MAGs?
motus.rel %>% 
  as_tibble(rownames='motus') %>% 
  pivot_longer(-motus) %>% 
  mutate(type=case_when(str_detect(motus, 'ref_mOTU')~'ref',
                        str_detect(motus, 'meta_mOTU')~'meta',
                        str_detect(motus, 'ext_mOTU')~'ext',
                        TRUE~'unclassified')) %>% 
  group_by(name, type) %>% 
  summarise(ab=sum(value)) %>% 
  mutate(site=str_remove(name, '[0-9]*$')) %>% 
  mutate(site.city=site.dict[site]) %>% 
  mutate(site.city=factor(site.city, levels = names(site.colours))) %>% 
  ggplot(aes(x=site.city, y=ab, fill=type)) + 
      geom_boxplot()

# answer: yes, definitely


# ##############################################################################
# for MPA4/Phanta: fraction of eukaryotes/viruses
 
# viruses
test.vkirus.phanta <- .f_site_assoc(
  enframe(lst.phanta.rel$superkingdom['Viruses',], 
          name='Sample', value='viral_read_abundance'),'viral_read_abundance')
ggsave(test.vkirus.phanta$plot + 
         ylab('Fraction of viral read abundance'), 
       filename=here('figures', 'general', 'viruses_phanta.pdf'),
       width = 5, height = 4, useDingbats=FALSE)
ggsave(test.vkirus.phanta$plot.ef + 
         theme(axis.text.x = element_text(angle=45, hjust=1)),
       filename=here('figures', 'general', 
                     'viruses_phanta_ef.pdf'),
       width = 4, height = 3, useDingbats=FALSE)

test.euk.phanta <- .f_site_assoc(
  enframe(lst.mpa$kingdom['k__Eukaryota',], 
          name='Sample', value='unclassified'),'unclassified')
ggsave(test.euk.map$plot + 
         ylab('Fraction of Eukaryotic abundance'), 
       filename=here('figures', 'general', 'euk_mpa4.pdf'),
       width = 5, height = 4, useDingbats=FALSE)
ggsave(test.euk.map$plot.ef + 
         theme(axis.text.x = element_text(angle=45, hjust=1)),
       filename=here('figures', 'general', 
                     'euk.mpa4_ef.pdf'),
       width = 4, height = 3, useDingbats=FALSE)

# ##############################################################################
# Buffer and Zymo

# Zymo stuff
zymo.cols <- c('Escherichia coli'="#FF0000", 
               'Salmonella enterica'="#00A08A", 
               'Bacillus subtilis'="#F2AD00", 
               'Pseudomonas aeruginosa'="#F98400", 
               'Enterococcus faecalis'="#5BBCD6", 
               'Lactobacillus fermentum'='#ECCBAE',
               'Limosilactobacillus fermentum'='#ECCBAE',
               'Staphylococcus aureus'="#046C9A", 
               'Listeria monocytogenes'="#D69C4E",
               'Saccharomyces cerevisiae'='#ABDDDE', 
               'Cryptococcus neoformans'='#8C1515',
               'unassigned'='#707273',
               'unclassified'='#707273',
               'other'='#D3D3D3')

ref.motus <- tibble(
  motus=c('Escherichia coli', 'Salmonella enterica', 'Bacillus subtilis',
          'Pseudomonas aeruginosa', 'Enterococcus faecalis', 
          'Lactobacillus fermentum', 'Staphylococcus aureus', 
          'Listeria monocytogenes', 'Saccharomyces cerevisiae', 
          'Cryptococcus neoformans'), name='reference', 
  rel.ab=c(rep(0.12, 8), 0.02, 0.02))
ref.mpa <- tibble(
  species=c('Escherichia coli', 'Salmonella enterica', 'Bacillus subtilis',
          'Pseudomonas aeruginosa', 'Enterococcus faecalis', 
          'Limosilactobacillus fermentum', 'Staphylococcus aureus', 
          'Listeria monocytogenes', 'Saccharomyces cerevisiae', 
          'Cryptococcus neoformans'), name='reference', 
  rel.ab=c(rep(0.12, 8), 0.02, 0.02))

motus.zymo <- prop.table(motus.buffer,2)
motus.zymo <- motus.zymo[rowMeans(motus.zymo==0)!=1,]

g.zymo.motus <- motus.zymo %>% 
  as_tibble(rownames = 'motus') %>% 
  pivot_longer(-motus) %>% 
  mutate(motus=str_remove(motus, ' \\[.*\\]')) %>% 
  mutate(motus=case_when(motus %in% names(zymo.cols)~motus, TRUE~'other')) %>% 
  group_by(motus, name) %>% 
  summarise(rel.ab=sum(value), .groups='drop') %>% 
  bind_rows(ref.motus) %>% 
  mutate(motus=factor(motus, levels=names(zymo.cols))) %>% 
  filter(str_detect(name, 'Zymo')) %>% 
  mutate(name=factor(name, levels = c('reference', 
                                      paste0('Zymo', seq_len(20))))) %>% 
  ggplot(aes(x=name, y=rel.ab, fill=motus)) + 
    geom_bar(stat='identity') + 
    theme_bw() + 
    xlab('') + ylab('Relative abundance') + 
    scale_fill_manual(values = zymo.cols) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + 
    ggtitle('mOTUs3')
ggsave(g.zymo.motus, filename = here('figures', 'general', 'zymo_motus.pdf'),
       width = 6, height = 4, useDingbats=FALSE)

# same for Metaphlan
mpa4.zymo <- mpa4.buffer[,str_detect(colnames(mpa4.buffer), 'Zymo')]
mpa4.zymo <- mpa4.zymo[rowMeans(mpa4.zymo==0)!=1,]
g.zymo.mpa <- mpa4.zymo %>% 
  as_tibble(rownames = 'species') %>% 
  pivot_longer(-species) %>% 
  mutate(species=str_replace(species, '_', ' ')) %>% 
  mutate(species=case_when(species %in% names(zymo.cols)~species, TRUE~'other')) %>% 
  group_by(species, name) %>% 
  summarise(rel.ab=sum(value), .groups='drop') %>% 
  bind_rows(ref.mpa) %>% 
  mutate(species=factor(species, levels = names(zymo.cols))) %>% 
  mutate(name=factor(name, levels = c('reference', 
                                      paste0('Zymo', seq_len(20))))) %>% 
  ggplot(aes(x=name, y=rel.ab, fill=species)) + 
  geom_bar(stat='identity') + 
  theme_bw() + 
  xlab('') + ylab('Relative abundance') + 
  scale_fill_manual(values=zymo.cols) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + 
  ggtitle('MetaPhlAn4')
ggsave(g.zymo.mpa, filename = here('figures', 'general', 'zymo_mpa4.pdf'),
       width = 6, height = 4, useDingbats=FALSE)

# phanta (new database)
phanta.zymo <- phanta.buffer[,str_detect(colnames(phanta.buffer), 'Zymo')]
phanta.zymo <- phanta.zymo[rowMeans(phanta.zymo==0)!=1,]
g.zymo.phanta <- as_tibble(phanta.zymo, rownames='species') %>% 
  pivot_longer(-species) %>% 
  group_by(name) %>% 
  mutate(value=value/sum(value)) %>%
  mutate(species=str_remove(species, '^.*\\|species_')) %>% 
  mutate(species=str_remove(species, '^s__')) %>% 
  mutate(species=str_replace(species, '_', ' ')) %>%
  mutate(species=case_when(species %in% names(zymo.cols)~species, 
                           TRUE~'other')) %>% 
  group_by(species, name) %>% 
  summarise(rel.ab=sum(value), .groups='drop') %>% 
  bind_rows(ref.mpa) %>% 
  mutate(species=factor(species, levels = names(zymo.cols))) %>% 
  mutate(name=factor(name, levels = c('reference', 
                                      paste0('Zymo', seq_len(20))))) %>% 
  ggplot(aes(x=name, y=rel.ab, fill=species)) + 
  geom_bar(stat='identity') + 
  theme_bw() + 
  xlab('') + ylab('Relative abundance') + 
  scale_fill_manual(values=zymo.cols) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + 
  ggtitle('Phanta')
ggsave(g.zymo.phanta, filename = here('figures', 'general', 
                                      'zymo_phanta.pdf'),
       width = 6, height = 4, useDingbats=FALSE)


# kraken (from Dylan)
kraken.zymo <- read.table(here('data','zymo','kraken_species_percentage.txt'), 
                          sep='\t', header=TRUE, row.names = 1, 
                          stringsAsFactors = FALSE, check.names = FALSE)
g.zymo.kraken <- kraken.zymo %>% 
  as_tibble(rownames='species') %>% 
  pivot_longer(-species) %>% 
  filter(species!='classified at a higher level') %>%
  group_by(name) %>% 
  mutate(value=value/sum(value)) %>% 
  mutate(species=case_when(species %in% names(zymo.cols)~species, 
                           TRUE~'other')) %>%
  group_by(species, name) %>% 
  summarise(rel.ab=sum(value), .groups='drop') %>% 
  bind_rows(ref.mpa) %>% 
  mutate(species=factor(species, levels = names(zymo.cols))) %>% 
  mutate(name=factor(name, levels = c('reference', 
                                      paste0('Zymo', seq_len(20))))) %>% 
  ggplot(aes(x=name, y=rel.ab, fill=species)) + 
  geom_bar(stat='identity') + 
  theme_bw() + 
  xlab('') + ylab('Relative abundance') + 
  scale_fill_manual(values=zymo.cols) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + 
  ggtitle('Kraken standard')
ggsave(g.zymo.kraken, filename = here('figures', 'general', 'zymo_kraken.pdf'),
       width = 6, height = 4, useDingbats=FALSE)
