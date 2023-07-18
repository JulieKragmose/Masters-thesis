#Relative barplot showing changes between major and minor

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(cowplot)

minDepth = 50
minQuality = 30
minCount = 15
minFreq = 0.025


df <- read_table('Z:/VOF/INF/JUKJ/Quasi/run/human/variantsWithMajor.txt')

#Augment
df_augmented <- df %>%
  select(Sample, Segment, Position, Base, Count, Depth, Frequency, Average_Quality, Alleletype) %>%
  mutate(Segment = recode(Segment, 
                          A_PB2 = 'PB2',
                          A_PB1 = 'PB1',
                          A_PA = 'PA',
                          A_HA_H3 = "HA",
                          A_NP = 'NP',
                          A_NA_N2 = "NA",
                          A_MP = 'MP',
                          A_NS = 'NS')) %>%
  mutate(Alleletype = recode(Alleletype,
                             Consensus = 'Major',
                             Minority = 'Minor')) %>%
  filter(Base != '-') %>% #Remove deletions
  mutate(Base = case_when(
          Alleletype == 'Major' & Depth < minDepth ~ 'N',
          Alleletype == 'Major' & Average_Quality < minQuality ~ 'N',
          Alleletype == 'Major' & Count < minCount ~ 'N',
          Alleletype == 'Major' & Frequency < minFreq ~ 'N',
          TRUE ~ Base
    )
  )




#Add Non-synonymous and synonymous
#mutationSummaryDf <- read_table('Z:/VOF/INF/JUKJ/Quasi/run/human/mutationSummary.txt') %>%
#  select(Sample, Segment, Nucleotide_Position, Mutation_type) %>%
#  rename(Position=Nucleotide_Position)

#df_augmented <- left_join(df_augmented, mutationSummaryDf, by=c('Sample', 'Segment', 'Position'))

#df_augmented <- df_augmented %>%
#  mutate(Mutation_type = recode(Mutation_type,
#                      "Non-synonymous" = 'NS',
#                      "Synonymous" = 'S')) %>%
#  mutate(Facet_Name = paste(Position, Mutation_type))



#Plot all variants (relative)
for (i in unique(df_augmented$Segment)){
  print(i)
  #Get data belonging to segment i
  segData <- filter(df_augmented, Segment == i)
  
  #Barplot
  p <- ggplot(segData, aes(x = factor(Sample, level = c('Day1_1', 'Day1_2', 'Day3', 'Day8', 'Day14', 'Day17', 'Day21')), 
                           y = Frequency, 
                           fill = Base)) +
    geom_bar(position = "fill", stat = "identity") +
    geom_text(aes(label = paste(round(Frequency*100, digits = 2), "%", sep="")),
              #vjust = -0.5, 
              size = 2.5,
              position = position_stack(vjust = 0.5),
              colour = "white") +
    scale_fill_manual(values = c('A' = "#990000", 'G' = "#F6D04D", 'T' = "#2F3EEA", 'C' = "#FC7634", 'N' = "#DADADA")) +
    labs(title = i, x = "Sample", y = "Relative frequency") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text = element_text(size=15),
          plot.title = element_text(size=22)) +
    facet_wrap(~Position)
  
  #Save plot
  ggsave(filename = paste("C:/Users/B281386/Desktop/Quasi/Figures/Relative/quasi_relativeFreq_", i ,".png", sep = ""),
         width = 12, height = 8)
  
}



#Plot all variants (non-relative)
for (i in unique(df_augmented$Segment)){
  print(i)
  #Get data belonging to segment i
  segData <- filter(df_augmented, Segment == i)
  
  #Barplot
  p <- ggplot(segData, aes(x = factor(Sample, level = c('Day1_1', 'Day1_2', 'Day3', 'Day8', 'Day14', 'Day17', 'Day21')), 
                           y = Frequency, 
                           fill = Base)) +
    geom_bar(position = "stack", stat = "identity") +
    geom_text(aes(label = paste(round(Frequency*100, digits = 2), "%", sep="")),
              #vjust = -0.5, 
              size = 2.5,
              position = position_stack(vjust = 0.5),
              colour = "white") +
    scale_fill_manual(values = c('A' = "#990000", 'G' = "#F6D04D", 'T' = "#2F3EEA", 'C' = "#FC7634", 'N' = "#DADADA")) +
    labs(title = i, x = "Sample", y = "Relative frequency") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text = element_text(size=15),
          plot.title = element_text(size=22)) +
    facet_wrap(~Position)
  
  #Save plot
  ggsave(filename = paste("C:/Users/B281386/Desktop/Quasi/Figures/Non-Relative/quasi_Freq_", i ,".png", sep = ""),
         width = 12, height = 8)
  
}




#Plot chosen variants individually (relative)
chosenVariants <- list(list('PA',1343), list('HA', 1359), list('HA', 927), list('NA', 277), list('NS', 193))
for(v in chosenVariants){
  seg <- v[[1]]
  pos <- v[[2]]
  variantData <- filter(df_augmented, Segment == seg, Position == pos)
  
  #Plot
  p1 <- ggplot(variantData, aes(x = factor(Sample, level = c('Day1_1', 'Day1_2', 'Day3', 'Day8', 'Day14', 'Day17', 'Day21')), 
                               y = Frequency, 
                               fill = Base)) +
    geom_bar(position = "fill", stat = "identity") +
    geom_text(aes(label = paste(round(Frequency*100, digits = 2), "%", sep="")),
              #vjust = -0.5, 
              size = 2.5,
              position = position_stack(vjust = 0.5),
              colour = "white") +
    scale_fill_manual(values = c('A' = "#990000", 'G' = "#F6D04D", 'T' = "#2F3EEA", 'C' = "#FC7634", 'N' = "#DADADA")) +
    labs(title = paste(seg, pos, sep = " "), x = "Sample", y = "Relative frequency") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text = element_text(size=15),
          plot.title = element_text(size=22))
  
  #Save plot
  ggsave(filename = paste("C:/Users/B281386/Desktop/Quasi/Figures/Relative/quasi_relativeFreq_", seg, pos ,".png", sep = ""),
         width = 12, height = 8)
}



#Plot chosen variants individually (non-relative)
for(v in chosenVariants){
  seg <- v[[1]]
  pos <- v[[2]]
  variantData <- filter(df_augmented, Segment == seg, Position == pos)
  
  #Plot
  p1 <- ggplot(variantData, aes(x = factor(Sample, level = c('Day1_1', 'Day1_2', 'Day3', 'Day8', 'Day14', 'Day17', 'Day21')), 
                                y = Frequency, 
                                fill = Base)) +
    geom_bar(position = "fill", stat = "identity") +
    geom_text(aes(label = paste(round(Frequency*100, digits = 2), "%", sep="")),
              #vjust = -0.5, 
              size = 2.5,
              position = position_stack(vjust = 0.5),
              colour = "white") +
    scale_fill_manual(values = c('A' = "#990000", 'G' = "#F6D04D", 'T' = "#2F3EEA", 'C' = "#FC7634", 'N' = "#DADADA")) +
    labs(title = paste(seg, pos, sep = " "), x = "Sample", y = "Relative frequency") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.text = element_text(size=15),
          plot.title = element_text(size=22))
  
  #Save plot
  ggsave(filename = paste("C:/Users/B281386/Desktop/Quasi/Figures/Non-Relative/quasi_Freq_", seg, pos ,".png", sep = ""),
         width = 12, height = 8)
}

