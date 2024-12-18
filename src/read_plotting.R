library(dplyr)
library(ggplot2)
library(forcats)
library(data.table)
library(ggpubr)
library(rstatix)

reads <- read.table("/Users/dylanmaghini/scg4/projects/awigen2/response/read_counts/read_counts.tsv", header = TRUE, sep="\t")
reads <- reads %>% mutate(sitecode = substr(Sample, 1, 2))

sitedf <- data.frame(sitecode = c("AG", "DM", "BF", "GH", "KY", "SW"), 
                     site = c("Agincourt", "DIMAMO", "Nanoro", "Navrongo", "Nairobi", "Soweto"))

pal <- c("#889A4D","#BE9856","#4A4C59","#82354E","#DBBF5A","#99ADD0")
names(pal) <- c("Agincourt", "Nanoro", "Nairobi", "Soweto", "DIMAMO", "Navrongo")

reads <- merge(sitedf, reads, by=c("sitecode"))

reads_long <- melt(setDT(reads), id.vars = 1:3, variable.name = "category")
reads_long$site <- as.factor(reads_long$site)
reads_long <- reads_long %>% mutate(site = fct_relevel(site, "Nanoro", "Navrongo", "DIMAMO", "Agincourt", "Soweto", "Nairobi"))

ggplot(reads_long %>% 
         filter(category %in% c("raw_reads", "dedup_reads", "trimmed_reads", "hostremoved_reads")), 
       aes(x = category, y = value, fill = site)) +
  geom_boxplot(outlier.size = 0.5, outlier.color = "grey", linewidth = 0.3) + 
  ylab("Number of Reads") + 
  scale_fill_manual(values = pal) + 
  scale_x_discrete(labels = c("Raw reads", "Deduplicated reads", "Trimmed reads", "Host-removed reads")) + 
  theme_bw() + 
  stat_kruskal_test(label = "{p.adj.format}", p.adjust.method = "hochberg", label.y = 1.0e+08) + 
  theme(axis.title.x = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_rect(color = "black", size = 1), 
        text = element_text(size = 10, color = "black"),
        axis.text = element_text(color = "black"), 
        legend.title = element_blank())


ggsave("/Users/dylanmaghini/scg4/projects/awigen2/response/read_counts/read_counts.pdf", dpi = 300, w = 7, h = 1.7)
ggsave("/Users/dylanmaghini/scg4/projects/awigen2/response/read_counts/read_counts.jpeg", dpi = 300, w = 7, h = 1.7)


write.table(reads_long %>% 
              filter(category %in% c("raw_reads", "dedup_reads", "trimmed_reads", "hostremoved_reads")), 
                     "/Users/dylanmaghini/Downloads/sourcedata/ed3a.tsv", 
                     sep = "\t", quote = FALSE, row.names = FALSE)
            