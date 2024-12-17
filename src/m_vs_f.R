# ##############################################################################
#
## Navrongo M vs F comparison
#
# ##############################################################################

library("tidyverse")
library("here")
library("vegan")
library("labdsv")

load(here('data', 'metadata', 'metadata_clean.RData'))
#load(here('data', 'classification', 'all_classification_tables.RData'))
load('./data/classification/new_db_tables.RData')


meta.all <- metadata.clean %>% 
  filter(general_site=='Navrongo') %>% 
  filter(!is.na(demographic_gender))

#motus.navrongo <- lst.motus.gtdb$species[,meta.all$general_sample_id]
motus.navrongo <- motus.gtdb.lvls$species[,meta.all$general_sample_id]
motus.navrongo <- motus.navrongo[rowMeans(motus.navrongo != 0) > 0.05,]

sex.colors <- c('#F28E2B', '#59A14F')

set.seed(1519)

# alpha
mat.rar <- rrarefy(t(motus.navrongo), 5000)
alpha <- map(c('invsimpson', 'simpson', 'shannon', 'richness'), 
             .f = function(x){
               if (x=='richness'){
                 enframe(rowSums(mat.rar != 0), name='metaG_id', 
                         value='alpha') %>% 
                   mutate(index=x)
               } else {
                 enframe(diversity(mat.rar, index=x), name='metaG_id', 
                         value='alpha') %>% 
                   mutate(index=x)
               }
}) %>% bind_rows() %>% 
  full_join(meta.all %>% rename(metaG_id=general_sample_id), by='metaG_id')

g1 <- alpha %>% 
  filter(index=='invsimpson') %>% 
  ggplot(aes(x=demographic_gender, y=alpha, fill=demographic_gender)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(width = 0.17, col='#D3D3D3') + 
    # facet_wrap(~index, scales = 'free_y') + 
    xlab('') + ylab('Prokaryotic diversity [Inverse Simpson]') + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major.x = element_blank()) + 
    scale_fill_manual(values=sex.colors, guide='none')
g1$data %>% 
  select(metaG_id, alpha, demographic_gender) %>% 
  write_csv('./files/source_data/ED2a.csv')
# test
map(unique(alpha$index), .f = function(x){
  t <- wilcox.test(alpha~demographic_gender, data=alpha %>% filter(index==x))
  tibble(index=x, p.val=t$p.value)}) %>% 
  bind_rows()
#' 
#' # A tibble: 4 × 2
#' index       p.val
#' <chr>       <dbl>
#' 1 invsimpson 0.0270
#' 2 simpson    0.0270
#' 3 shannon    0.0213
#' 4 richness   0.0735

# beta
pco.res <- pco(vegdist(mat.rar))
df.plot <- as.data.frame(pco.res$points) %>% 
  as_tibble(rownames = 'metaG_id') %>% 
  full_join(meta.all %>% rename(metaG_id=general_sample_id), by='metaG_id')
var.explained <- sprintf(fmt='%.2f', 100*pco.res$eig[1:2]/sum(
  pco.res$eig[pco.res$eig > 0]))
g2 <- df.plot %>% 
  ggplot(aes(x=V1, y=V2, col=demographic_gender)) + 
    geom_point() +
    scale_colour_manual(values=sex.colors, guide='none') + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) + 
    xlab(paste0('PCo 1 [ ', var.explained[1], '%]')) + 
    ylab(paste0('PCo 2 [ ', var.explained[2], '%]'))
adonis2(vegdist(mat.rar)~demographic_gender, data=meta.all)
#' adonis2(formula = vegdist(mat.rar) ~ demographic_gender, data = meta.all)
#' Df SumOfSqs      R2      F Pr(>F)   
#' demographic_gender   1    0.511 0.00963 2.2472  0.004 **
#'   Residual           231   52.570 0.99037                 
#'  Total              232   53.081 1.00000                 
#' ---
#' Signif. codes:  0 *** 0.001 ** 0.01 * 0.05 ‘.’ 0.1 ‘ ’ 1


# diff abundance
rel.motus <- prop.table(motus.navrongo, 2)
hist(log10(rel.motus), 100)
log.n0 <- 1e-04
pb <- progress::progress_bar$new(total=nrow(rel.motus))
df.enrich <- map(rownames(rel.motus), .f = function(x){
  tmp <- enframe(rel.motus[x,], name='metaG_id') %>% 
    full_join(meta.all %>% rename(metaG_id=general_sample_id), 
              by='metaG_id') %>% 
    mutate(value=log10(value + log.n0))
  t <- wilcox.test(value~demographic_gender, data=tmp)
  gfc <- tmp %>% group_by(demographic_gender) %>% 
    group_map(.f=function(.x, .y){
      enframe(quantile(.x$value, 
                       probs=seq(from=0.05, to=0.95, by=0.05))) %>% 
        mutate(.y)}) %>% 
    bind_rows() %>% 
    pivot_wider(names_from = demographic_gender, values_from = value) %>% 
    mutate(diff=male-female) %>% 
    pull(diff) %>% mean
  pb$tick()
  tibble(motus=x, p.val=t$p.value, ef=gfc)}) %>% 
  bind_rows() %>% 
  mutate(p.adj = p.adjust(p.val, method='fdr'))
g3 <- df.enrich %>% 
  mutate(label=case_when(p.adj < 0.05~motus, TRUE~'')) %>% 
  ggplot(aes(x=ef, y=-log10(p.adj))) + 
    geom_point() + 
    xlab('Effect size') + 
    ylab('-log10(q-value)') + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank()) + 
    ggrepel::geom_text_repel(aes(label=label)) + 
    ylim(0, 5) + 
    geom_hline(yintercept = c(-log10(0.05), -log10(0.01)))
g3$data %>% 
  write_csv('./files/source_data/ED2b.csv')

g4 <- as_tibble(rel.motus[df.enrich %>%
                            filter(abs(ef)> 0.5) %>%
                            pull(motus),], rownames='motus') %>%
  pivot_longer(-motus, names_to = 'metaG_id') %>%
  full_join(meta.all %>% rename(metaG_id=general_sample_id), by='metaG_id') %>%
  mutate(value=log10(value + log.n0)) %>%
  ggplot(aes(x=demographic_gender, y=value, fill=demographic_gender)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.17) +
    facet_grid(~motus) +
    theme_bw() + theme(panel.grid.minor = element_blank(),
                       panel.grid.major.x=element_blank()) +
    scale_fill_manual(values=sex.colors, guide='none') +
    xlab('') + ylab('log10(Relative abundance)')

g.all <- cowplot::plot_grid(g1, g3, ncol=2, rel_widths = c(0.4, 0.6), 
                            labels = 'auto')
ggsave(g.all, filename = here('figures', 'general', 'm_vs_f_newdb.pdf'),
       width = 8, height = 5, useDingbats=FALSE)
