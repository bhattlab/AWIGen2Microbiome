library(dplyr)
library(gggenes)
library(ggplot2)
library(readxl)
library(scales)
library(forcats)
library(paletteer)
library(tidyr)
library(pheatmap)

phage_annot <- read_excel("/Users/dylanmaghini/scg4/projects/awigen2/phages/response/bigphages/bakta/out/allphages.xlsx",sheet=1,col_names= TRUE,col_types=NULL,skip= 0)
names(phage_annot) <- c("Phage", "Contig", "Type", "Start", "Stop", "Direction", "ID1", "ID2", "Annot", "Annot2", "AnnotCat")
phage_annot$Start <- as.numeric(phage_annot$Start)
phage_annot$Stop <- as.numeric(phage_annot$Stop)
phage_annot <- phage_annot %>% mutate(Direction = ifelse(Direction == "+", TRUE, FALSE))
phage_annot <- phage_annot %>% filter(Phage != "GH039.tsv")
phage_annot <- phage_annot %>% mutate(AnnotCat = ifelse(AnnotCat == "annotated", "hypothetical", AnnotCat))

pal <- paletteer_d("NatParksPalettes::Torres")
pal <- c("#cccccc", pal[1:4], pal[7:10], pal[5])
names(pal) <- c("hypothetical", 
                "CRISPR array", "Methyltransferase", 
                "Sporulation regulating", 
                "Recombinase", 
                "Addiction Module Toxin", 
                "Cellulase", 
                "Penicillin-binding protein", "ABC transporter", 
                "tRNA")

phage_annot <- phage_annot %>% mutate(AnnotCat = fct_relevel(AnnotCat,
                                              "CRISPR array", "Methyltransferase", 
                                              "Sporulation regulating", 
                                              "Recombinase", 
                                              "Addiction Module Toxin", 
                                              "Cellulase", 
                                              "Penicillin-binding protein", "ABC transporter", 
                                              "tRNA", "hypothetical",))

length_order <- data.frame(Phage = c("SW076.tsv", "BF107.tsv", "BF087.tsv", "AG414.tsv", "BF270.tsv", "AG493.tsv", "BF122.tsv", "BF036.tsv", "BF152.tsv"), 
                           Order = c("Phage A", "Phage B", "Phage C", "Phage D", "Phage E", "Phage F", "Phage G", "Phage H", "Phage I"))

phage_annot <- merge(phage_annot, length_order, by = "Phage")
ggplot(phage_annot, aes(xmin = Start, xmax = Stop, y = Order, forward = Direction, color = AnnotCat, fill = AnnotCat)) +
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
  facet_wrap(~ Order, ncol = 1, scales = "free_y") +
  scale_color_manual(values = pal) + 
  scale_x_continuous(labels = label_comma()) + 
  scale_fill_manual(values = pal) + 
  theme_bw() + 
  theme(strip.background = element_blank(),
          strip.text.x = element_blank(), 
        axis.title.y = element_blank(), 
        panel.border = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        #axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        panel.grid.major.y= element_blank(), 
        legend.text = element_text(size=10)
        #legend.position = "none"
        )
ggsave("/Users/dylanmaghini/scg4/projects/awigen2/phages/response/bigphages/genomemap.pdf", dpi = 300, w = 5.3, h = 2)
ggsave("/Users/dylanmaghini/scg4/projects/awigen2/phages/response/bigphages/genomemap_wider.pdf", dpi = 300, w = 7, h = 2)

write.table(phage_annot, "/Users/dylanmaghini/Downloads/sourcedata/fig3e.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)

ggplot(phage_annot, aes(x = 1, y = 1, color = AnnotCat)) + 
  geom_point() + 
  scale_color_manual(values = pal) + 
  theme_bw() + 
  theme(legend.position = "bottom")

ggsave("/Users/dylanmaghini/scg4/projects/awigen2/phages/response/bigphages/genomemap_legend.pdf", dpi = 300, w = 8, h = 4)

phage_summary <- phage_annot %>% group_by(Order, AnnotCat, .drop = FALSE) %>% 
  summarise(count = n())

heatmap_pal <- paletteer_d("RColorBrewer::BuGn")
heatmap_pal = c("#ffffff", heatmap_pal[3:9])
phage_summary_wide <- data.frame(pivot_wider(phage_summary, names_from = AnnotCat, values_from = count))
rownames(phage_summary_wide) <- phage_summary_wide$Order
phage_summary_wide <- phage_summary_wide %>% select(!Order) %>% select(!hypothetical)
phage_summary_wide <- as.matrix(phage_summary_wide)


paletteLength <- 15
temppal <- colorRampPalette(paletteer_d("RColorBrewer::BuGn"))(30)
myColor <- c("#ffffff", colorRampPalette(c("lightblue", "darkblue"))(16))
myColor <- c("#ffffff", temppal[7:23])

breaks <- c(-1:16)

pheatmap(phage_summary_wide, cluster_rows = FALSE, cluster_cols = FALSE, 
         display_numbers = TRUE, number_format = "%.0f", color = myColor, breaks=breaks, 
         border_color = NA, angle_col = "45", legend = FALSE, cellwidth = 15, cellheight = 15,
         filename = "/Users/dylanmaghini/scg4/projects/awigen2/phages/response/bigphages/geneheatmap.pdf", width = 6, height = 4)
