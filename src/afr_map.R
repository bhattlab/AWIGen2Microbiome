# ##############################################################################
#
## General part of figure 1: cohort description
#
# ##############################################################################

library("tidyverse")
library("maps")
library("maptools")
library("yaml")
library("here")

plot.params <- yaml.load_file(here('files', 'params.yml'))
site.dict <- unlist(plot.params$sites)
site.colours <- unlist(plot.params$sites_colours)

data(wrld_simpl)
africa <- wrld_simpl[wrld_simpl$REGION==2,]

df.list <- list()
for (x in seq_len(length(africa@polygons))){
  tmp <- as.data.frame(africa@polygons[[x]]@Polygons[[1]]@coords) %>% 
    as_tibble() %>% 
    mutate(country=as.character(africa$NAME[x])) 
  tmp <- tmp %>% 
    mutate(order=seq_len(nrow(.)))
  df.list[[x]] <- tmp
}
df.plot <- df.list %>% 
  bind_rows() %>% 
  mutate(included_country=country %in% c('South Africa', 'Ghana', 
                                         'Kenya', 'Burkina Faso')) %>% 
  arrange(desc(included_country)) %>% 
  mutate(country=factor(country, levels = unique(country)))


# location of the study sites
# taken from google maps (might be pretty imprecise, but should work on a 
# map of this resolution)

study_sites <- tibble(
  site=c('Soweto', 'Agincourt', 'Dimamo', 'Nanoro', 'Navrongo', 'Nairobi'),
  V1=c(27.52, 31.26, 29.74, -2.19, -1.09, 36.88),
  V2=c(-26.16, -24.82, -23.88, 12.69, 10.88, -1.28)
)

g <- ggplot(df.plot, aes(V1, V2)) +
  geom_polygon(aes(group = country, fill=included_country), 
               colour='#EDEDED') +
  coord_map("bonne", lat0=0) + 
  scale_fill_manual(values=c('#EDEDED', '#fae0e0'), guide='none') + 
  theme_void() + 
  geom_point(data=study_sites, aes(colour=site)) + 
  ggrepel::geom_text_repel(data=study_sites, aes(label=site, colour=site)) + 
  scale_colour_manual(values=unlist(site.colours), guide='none')
ggsave(g, filename = here('figures', 'general', 'site_map.pdf'),
       width = 5, height = 6, useDingbats=FALSE)


# ##############################################################################
# Number of participants per site

df.quality <- read_tsv(here('data', 'misc', 'read_counts.tsv'),
                       col_types = cols()) %>% 
  mutate(site=str_remove(Sample, '[0-9]*$')) %>% 
  mutate(site.city=site.dict[site]) %>% 
  mutate(site.city=factor(site.city, levels = names(site.colours)))

g.n <- df.quality %>% 
  group_by(site.city) %>% 
  tally() %>% 
  filter(!is.na(site.city)) %>% 
  ggplot(aes(x=site.city, y=n, fill=site.city)) + 
    geom_bar(stat='identity') + 
    scale_fill_manual(values=unlist(site.colours), guide='none') + 
    theme_bw() + coord_flip() +
    theme(panel.grid.major.y = element_blank(), 
          panel.grid.minor = element_blank()) + 
    ylab('Number of participants') + xlab('')
ggsave(g.n, filename = here('figures', 'general', 'no_participants.pdf'),
       width = 4, height = 3, useDingbats=FALSE)
