#!/usr/bin/env Rscript 
args <- commandArgs(trailingOnly = TRUE)
set.seed(28123)


# read parameters
min_year <- as.numeric(args[1])
max_year <- as.numeric(args[2])
name <- args[3]
thresh_color <- as.numeric(args[4])
software_path <- args[5]
first_season <- args[6]
filter <- args[7]
width <- as.numeric(args[8])
plotFormat <- args[9]



source(paste(software_path,"/buildAD_sign_fun.R", sep=""))
buildAD(infile.prefix = name, cutoff = 0, thresh_color=thresh_color, 
        minyear=min_year, maxyear=max_year, drucke=TRUE, startingSeason=first_season,
	filter_sign=filter, plotwidth=width, format=plotFormat )
#source(paste(args[5],"R/plot_seasons.R", sep=""))
#plot_seasons(infile.prefix = args[3], cutoff = cutoff, thresh_color=.5, 
#        minyear=min_year, maxyear=max_year)


#-p -n ../../Analyse/segment1/ADPlot_human/humantree_bifurcating -t ../../Analyse/segment1/ADPlot_human/human_aa.fa -m ../../Analyse/segment1/ADPlot_human/human_aa.map -f AccTran -type PROTEIN -plot 0 -time bc
#-p -n ../../Analyse/segment1/ADPlot_human/humantree_bifurcating -t ../../Analyse/segment1/ADPlot_human/human_aa.fa -m ../../Analyse/segment1/ADPlot_human/human_aa.map -f AccTran -type PROTEIN -plot 0 -time bc
#-p -n ../../gisaid_analysis/segment1/ADPlot/humantree_bifurcating -t ../../gisaid_analysis/segment1/ADPlot/alignedAA.fa -m ../../gisaid_analysis/segment1/aa.map -f AccTran -type PROTEIN -time bw -plot 1
