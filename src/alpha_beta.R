# ##############################################################################
#
## Alpha and beta diversity plots
#
# ##############################################################################

library('tidyverse')
library('here')
library("yaml")
library("vegan")

plot.params <- yaml.load_file(here('files', 'params.yml'))
site.dict <- unlist(plot.params$sites)
site.colours <- unlist(plot.params$sites_colours)

# ##############################################################################
# functions to automate alpha/beta calculation/plotting

# alpha
.f_alpha <- function(mat, l=7500, rarefy=TRUE){
  if (rarefy){
    mat.rar <- rrarefy(t(mat), l)
  } else {
    mat.rar <- t(mat)
  }
  alpha.inv <- diversity(mat.rar, index='invsimpson')
  alpha.simp <- diversity(mat.rar, index='simpson')
  alpha.shan <- diversity(mat.rar, index='shannon')
  alpha.rich <- rowSums(mat.rar != 0)
  df <- enframe(alpha.inv, name='Sample_ID', value='inv_simpson') %>% 
    full_join(enframe(alpha.simp, name='Sample_ID', value='simpson'), 
              by='Sample_ID') %>% 
    full_join(enframe(alpha.shan, name='Sample_ID', value='shannon'), 
              by='Sample_ID') %>% 
    full_join(enframe(alpha.rich, name='Sample_ID', value='richness'), 
              by='Sample_ID')  
  g <- df %>% 
    filter(!str_detect(Sample_ID, '(Buffer|Zymo)')) %>% 
    mutate(site=str_remove(Sample_ID, '[0-9]*$')) %>% 
    mutate(site.city=site.dict[site]) %>% 
    mutate(site.city=factor(site.city, levels = names(site.colours))) %>% 
    pivot_longer(cols = c(inv_simpson, simpson, shannon, richness)) %>% 
    ggplot(aes(x=site.city, y=value)) + 
    geom_jitter(width = 0.2, colour='#D3D3D3') +
    geom_boxplot(aes(fill=site.city), outlier.shape = NA) +
    scale_fill_manual(values=alpha(site.colours, 0.85), 
                      name='Site', guide='none') + 
    scale_colour_manual(values=site.colours, name='Site', guide='none') + 
    facet_grid(name~., scales = 'free') +
    xlab('') + theme_bw() + ylab("Alpha diversity") +
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major.x = element_blank())
  # correlation between measures
  c.spear <- cor(df %>% select(-Sample_ID) %>% as.matrix()) 
  # maybe pearson is better than spearman?
  return(list('alpha_df'=df, 'plot'=g, 'correlation'=c.spear))
}

# beta
.f_beta <- function(mat, dist='log-euclidean', log.n0=1e-05){
  if (dist=='log-euclidean'){
    x <- log10(mat + log.n0)
    d.mat <- vegan::vegdist(t(x), method='euclidean')
  } else if (dist=='bray'){
    d.mat <- vegan::vegdist(t(mat), method='bray')
  } else if (dist=='jaccard'){
    d.mat <- vegan::vegdist(t(mat), method='jaccard')
  } else if (dist=='dist'){
    d.mat <- as.dist(mat)
  } else {
    stop("Parameter 'dist' should be either 'log-euclidean' or 'bray'!")
  }
  # dist.mat <- as.matrix(d.mat)
  # dist.mat[lower.tri(dist.mat)] <- NA
  # diag(dist.mat) <- NA
  # df.distances <- as_tibble(dist.mat, rownames='Sample') %>% 
  #   pivot_longer(-Sample) %>% 
  #   filter(!is.na(value)) %>% 
  #   mutate(siteA=str_remove(Sample, '[0-9]*$')) %>% 
  #   filter(!str_detect(siteA, '(Buffer|Zymo)')) %>%
  #   mutate(siteB=str_remove(name, '[0-9]*$')) %>% 
  #   filter(!str_detect(siteB, '(Buffer|Zymo)')) %>% 
  #   mutate(site.cityA=site.dict[siteA]) %>% 
  #   mutate(site.cityB=site.dict[siteB])
  
  pco <- labdsv::pco(d.mat)
  df <- pco$points
  colnames(df) <- c('PCo1', 'PCo2')
  n.var <- pco$eig[1:2]/sum(pco$eig[pco$eig > 0]) * 100
  df <- as_tibble(df, rownames='Sample') %>% 
    mutate(site=str_remove(Sample, '[0-9]*$')) %>% 
    filter(!str_detect(site, '(Buffer|Zymo)')) %>% 
    mutate(site.city=site.dict[site]) %>% 
    mutate(site.city=factor(site.city, levels = names(site.colours)))
  
  g1 <- df %>% 
    ggplot(aes(x=PCo1, y=PCo2, col=site.city)) + 
    geom_point() + 
    scale_colour_manual(values=site.colours, name='Site', guide='none') + 
    theme_classic() + 
    xlab(paste0('PCo 1 [', sprintf(fmt='%.2f', n.var[1]),'%]')) + 
    ylab(paste0('PCo 2 [', sprintf(fmt='%.2f', n.var[2]),'%]'))
  g.boxes.top <- df %>% 
    ggplot(aes(x=site.city, y=PCo1, fill=site.city)) + 
    geom_boxplot() + 
    theme_classic() + 
    coord_flip() + 
    xlab('') + ylab('') + 
    theme(axis.ticks = element_blank(), axis.text=element_blank(), 
          axis.line.x=element_blank()) + 
    scale_fill_manual(values=site.colours, guide='none')
  g.boxes.side <- df %>% 
    ggplot(aes(x=site.city, y=PCo2, fill=site.city)) + 
    geom_boxplot() + 
    theme_classic() + 
    xlab('') + ylab('') +
    theme(axis.ticks = element_blank(), axis.text=element_blank(), 
          axis.line.y=element_blank()) + 
    scale_fill_manual(values=site.colours, guide='none')
  g.all <- cowplot::plot_grid(g.boxes.top, NULL, g1, g.boxes.side, 
                              rel_heights = c(0.2, 0.8),
                              rel_widths = c(0.8, 0.2), align = 'hv')
  return(list('plot'=g.all, 'coords'=df))#, 'distances'=df.distances))
}

# load the data
load(here('data', 'classification', 'all_classification_tables.RData'))

# ##############################################################################
# alpha

alpha.motus <- .f_alpha(lst.motus$species, l=7500, rarefy=TRUE)
ggsave(alpha.motus$plot, filename = here('figures', 'classification_analyses',
                                         'alpha_motus.pdf'),
       width = 5, height = 8, useDingbats=FALSE)
alpha.mpa4 <- .f_alpha(lst.mpa$species, rarefy=FALSE)
ggsave(alpha.mpa4$plot, filename = here('figures', 'classification_analyses',
                                         'alpha_mpa4.pdf'),
       width = 5, height = 8, useDingbats=FALSE)
alpha.phanta <- .f_alpha(lst.phanta.tax$species, rarefy=FALSE)
ggsave(alpha.phanta$plot, filename = here('figures', 'classification_analyses',
                                          'alpha_phanta.pdf'),
       width = 5, height = 8, useDingbats=FALSE)
alpha.phanta.vir <- .f_alpha(lst.phanta.vir$species, rarefy=FALSE)
ggsave(alpha.phanta.vir$plot, filename = here('figures', 'classification_analyses',
                                          'alpha_phanta_vir.pdf'),
       width = 5, height = 8, useDingbats=FALSE)

# ##############################################################################
# compare alpha diversity across tools

.f_pairwise_scatter <- function(df){
  stopifnot(all(c('motus', 'mpa', 'phanta', 'phanta_vir') %in% colnames(df)))
  bind_rows(df %>% transmute(x=motus, y=mpa) %>% mutate(xt='motus', yt='mpa'),
    df %>% transmute(x=motus, y=phanta) %>% mutate(xt='motus', yt='phanta'),
    df %>% transmute(x=motus, y=phanta_vir) %>% 
      mutate(xt='motus', yt='phanta_vir'),
    df %>% transmute(x=mpa, y=phanta) %>% mutate(xt='mpa', yt='phanta'),
    df %>% transmute(x=mpa, y=phanta_vir) %>% mutate(xt='mpa', yt='phanta_vir'),
    df %>% transmute(x=phanta, y=phanta_vir) %>% 
      mutate(xt='phanta', yt='phanta_vir')) %>% 
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
  bind_rows(alpha.phanta$alpha_df %>% pivot_longer(-Sample_ID) %>% 
              mutate(type='phanta')) %>% 
  bind_rows(alpha.phanta.vir$alpha_df %>% pivot_longer(-Sample_ID) %>% 
              mutate(type='phanta_vir')) 

pdf(here('figures', 'classification_analyses', 'alpha_comp.pdf'),
    width = 8, height = 8, useDingbats = FALSE)
for (i in c('inv_simpson', 'richness', 'shannon', 'simpson')){
  print(.f_pairwise_scatter(df.alpha.comp %>% 
    filter(name==i) %>% 
    pivot_wider(names_from = 'type', values_from = 'value')) +
      ggtitle(i))
}
dev.off()  

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
      x <- .f_beta(t(vegan::rrarefy(t(lst.motus[[lvl]]), 7500)), dist=d)
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
for (d in c('bray', 'jaccard', 'log-euclidean')){
  pdf(here('figures', 'classification_analyses', 
           paste0('beta_phanta_', d,'.pdf')),
      width = 6, height = 5, useDingbats = FALSE)
  for (lvl in names(lst.phanta.tax)){
    x <- .f_beta(lst.phanta.tax[[lvl]], dist=d, log.n0 = 1e-06)
    print(x$plot)
    message(lvl)
  }
  dev.off()
}
for (d in c('bray', 'jaccard', 'log-euclidean')){
  pdf(here('figures', 'classification_analyses', 
           paste0('beta_phanta_vir_', d,'.pdf')),
      width = 6, height = 5, useDingbats = FALSE)
  for (lvl in names(lst.phanta.vir)){
    x <- .f_beta(lst.phanta.vir[[lvl]], dist=d, log.n0 = 1e-06)
    print(x$plot)
    message(lvl)
  }
  dev.off()
}


# ##############################################################################
# sourmash

# ##############################################################################
# compare distances
