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
study_sites %>% 
  write_csv('./files/source_data/Fig1a.csv')

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
g.n$data %>% 
  write_csv('./files/source_data/Fig1a2.csv')
ggsave(g.n, filename = here('figures', 'general', 'no_participants.pdf'),
       width = 4, height = 3, useDingbats=FALSE)


# ##############################################################################
# Number of samples per country in UHGG
## for comparison and talks

uhgg.meta <- 'http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v2.0.1/genomes-all_metadata.tsv'

uhgg.data <- read_tsv(uhgg.meta)
country.counts <- uhgg.data %>% 
  filter(Genome_type=='MAG') %>% 
  select(Country, Continent, Sample_accession) %>% 
  distinct() %>% 
  group_by(Country, Continent) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(Country=case_when(Country=='United Republic of Tanzania'~'Tanzania',
                           TRUE~Country)) %>% 
  filter(Country!='NA') %>% 
  rename("region"=Country)


library("sf")        # for manipulation of simple features objects
library("lwgeom")    # for st_transform_proj()
library("rworldmap") # for getMap()


world_sf <- st_as_sf(getMap(resolution = "low"))
world_sf <- world_sf[world_sf$continent!='Antarctica',]
world_sf$n <- 0
world_sf$n[match(country.counts$region, 
                         world_sf$NAME)] <- country.counts$n
world_sf$type_fill <- case_when(world_sf$n>10000~'>10000',
                                world_sf$n>1000~'>1000',
                                world_sf$n>100~'>100',
                                world_sf$n>=10~'>10',
                                TRUE~'0')
world_sf$type_fill <- factor(world_sf$type_fill, levels = c('>10000', '>1000', 
                                                            '>100', '>10', '0'))

crs_wintri <- "+proj=wintri +datum=WGS84 +no_defs +over"
grat_wintri <- 
  st_graticule(lat = c(-89.9, seq(-80, 80, 20), 89.9)) %>%
  st_transform_proj(crs = crs_wintri)

world_wintri <- st_transform_proj(world_sf, crs = crs_wintri)


g <- ggplot(world_wintri) + 
  geom_sf(data = grat_wintri, color = "#D3D3D3", size = 0, linewidth=0.15) +
  geom_sf(col='white', aes(fill=type_fill), linewidth=0.15) +
  coord_sf(datum = NULL) + 
  ggthemes::theme_map() +
  theme(legend.position = 'top') + 
  scale_fill_manual(values=c('#8C1515', '#007C92', '#44e3ff',  # '#8C1515',
                             #'#b9f4ff', 
                             '#CCCCCC'), name='Number of samples') 
ggsave(g, filename = './figures/general/no_samples_UHGG.pdf',
       width = 7, height = 5, useDingbats=FALSE)

