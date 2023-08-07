setwd("Z:/VOF/INF/JUKJ/Quasi/consDist/")
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(seqinr)
library(ape)
library(patchwork)
library(ggplotify)
library(pheatmap)
library(cowplot)

list = c(list.files(path = ".", pattern = "*_aln.fna"))
for (file in list) {
  #load fasta file (already aligned)
  myseqs <- read.alignment(file = file, format = "fasta")
  # Convert the alignment to "DNAbin" format
  testbin <- as.DNAbin(myseqs) 
  # Calculate the genetic distance matrix (pairwise = number of dissimilarity)
  testdist_pair <- as.matrix(dist.gene(testbin, method = "pairwise"))
  # Calculate the genetic distance matrix (percentage = percentage of dissimilarity)
  testdist_perc <- as.matrix(dist.gene(testbin, method = "percentage"), gap = T)
  # Decide the length of the palette
  paletteLength <- 10
  # Create color palette
  myColor <- colorRampPalette(brewer.pal(8, "Set3"))(paletteLength)
  # length(breaks) == length(paletteLength) + 1
  # use floor and ceiling to deal with even/odd length pallettelengths
  myBreaks_pair <- unique(c(seq(min(testdist_pair), 0.0001, length.out=ceiling(paletteLength/2) + 1), 
                        seq(max(testdist_pair)/paletteLength, max(testdist_pair), length.out=floor(paletteLength/2))))
  myBreaks_perc <- unique(c(seq(min(testdist_perc), 0.0001, length.out=ceiling(paletteLength/2) + 1), 
                   seq(max(testdist_perc)/paletteLength, max(testdist_perc), length.out=floor(paletteLength/2))))
  name = unlist(strsplit(sub('(^[^_]+)_(.*)$', '\\1', file), ' '))[1]
  name2 = paste("Pairwise_distance_", name, ".png", sep = "")

  # Plot the heatmap
  a =pheatmap(testdist_pair, color = myColor, breaks = myBreaks_pair, number_format = "%.0f", display_numbers = T, width = 20, height = 15, main = paste("Pairwise distance", name, "Segment",sep = " "))
  b= pheatmap(testdist_perc, color = myColor, breaks = myBreaks_perc, number_format = "%.3f", display_numbers = T, main = paste("Percentage of dissimilarity", name, "Segment", sep = " "))
  
  p1=as.ggplot(a)
  p2=as.ggplot(b)
  
  pp = plot_grid(p1,p2,  nrow = 2,
                 ncol = 1,
                 rel_widths = 20,
                 rel_heights = 30)
  #pp = ggarrange(p1, p2, nrow = 2, ncol =1)
  #pp = p1 + p2 + plot_layout(design = layout)
  ggsave(pp, filename=name2, width=15, height=15)
}

