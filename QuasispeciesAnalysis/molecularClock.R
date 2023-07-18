#Molecular clock analysis

library(ape)
library(tidyverse)
library(stringr)
library(readr)
library(ggplot2)
library(ggpmisc)

segments <- c('PB2', 'PB1', 'PA', 'HA_H3', 'NP', 'NA_N2', 'MP', 'NS')

for(segment in segments){
  treeFile <- paste(segment, "_nucleotide.nexus", sep="")
  
  #Load tree
  tree <- read.nexus(treeFile)
  
  #Compute the distance from root to tips
  N <- Ntip(tree)
  rootNode <- N + 1
  rootTipDistances = as.data.frame(dist.nodes(tree)[1:N, rootNode])
  rownames(rootTipDistances) <- tree$tip.label
  colnames(rootTipDistances)[1] <- "Root2Tip"
  
  #Extract patient data (not the root)
  patientSamples <- rootTipDistances %>%
    filter(grepl('Day',  rownames(.)))
  
  #Add days since first sample as column
  Days <- patientSamples %>%
    rownames() %>%
    str_extract("Day[0-9]+") %>%
    parse_number()
  patientSamples <- cbind(patientSamples, Days)
  
  
  prefix <- unlist(str_split(treeFile, '[.]'))[1]
  
  
  #Plot with equation
  p1 <- ggplot(data = patientSamples, aes(x=Days, y=Root2Tip)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red", formula = y ~ x) +
    stat_poly_eq(formula = y ~ x, aes (label = paste(..eq.label.., ..rr.label..,sep = "~~~")), parse = TRUE, size=3) +
    labs(title = segment,
         x = "Day since first sample",
         y = "Root to tip distance") +
    theme_classic()
  #print(p1)
  ggsave(filename = paste("../molecularClock/", prefix, "_molecularClock.png", sep=""), p1)
  