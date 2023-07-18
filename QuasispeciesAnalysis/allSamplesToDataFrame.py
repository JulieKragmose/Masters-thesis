# !/usr/bin/env python3
# Usage: python3 allSamplesToDataFrame.py

# Makes a combined table of all -variants.txt and -allAlleles.txt
import os, sys


#Initialize
samples = ['Day1_1', 'Day1_2', 'Day3', 'Day8', 'Day14', 'Day17', 'Day21']
segments = ['PB2', 'PB1', 'PA', 'HA_H3', 'NP', 'NA_N2', 'MP', 'NS']


#Get a dict of just the segment and positions for significant variants 
variantDict = dict()
#labelDict = dict()
for sample in samples:
    path = '/srv/data/VOF/INF/JUKJ/Quasi/run/human/' + sample + '/tables/'
    for segment in segments:
        #Open variant file
        try:
            variantFile = open(path + 'A_' + segment + '-variants.txt', 'r')
        except:
            print('No variants file found for ' + segment + ' in ' + sample)
            continue

        #Find significant minor variants
        for line in variantFile:
            if not line.startswith('Reference_Name'): #Skip header
                splitLine = line.split('\t')               
                depth = splitLine[2]
                frequency = splitLine[8]
                
                if int(depth) >= 100 and float(frequency) >= 0.025:
                    position = splitLine[1]
                    if segment in variantDict:
                        variantDict[segment].append(position)
                    else:
                        variantDict[segment] = [position]

                    #majorBase = splitLine[3]
                    #minorBase = splitLine[4]
                    #label = majorBase + str(position) + minorBase
                    #if position in labelDict:
                    #    labelDict[position].append(label)
                    #else:
                    #    labelDict[position] = [label]
        variantFile.close()



#Initialize output file
outfile = open('variantsWithMajor.txt', 'w')
outfile.write('Sample\tSegment\tPosition\tBase\tCount\tDepth\tFrequency\tAverage_Quality\tConfidenceNotMacErr\tPairedUB\tQualityUB\tAlleletype\n')


#Read allAlleles file and get data all data for these sites
for sample in samples:
    path = '/srv/data/VOF/INF/JUKJ/Quasi/run/human/' + sample + '/tables/'
    for segment in segments:
        try:
            allAllelesFile = open(path + 'A_' + segment + '-allAlleles.txt', 'r')
        except:
            print('No allAlleles file found for ' + segment + ' in ' + sample)
            continue

        for line in allAllelesFile:
            if not line.startswith('Reference_Name'): #Skip header
                splitLine = line.split('\t')   
                position = splitLine[1]

                #Check if this is a position of interest
                if segment in variantDict and position in variantDict[segment]:
                    alleleType = splitLine[-1][:-1]

                    #If it's the major allele, then we definitely want it
                    if alleleType == 'Consensus':
                        outfile.write(sample + '\t' + line)
                    
                    #If it's a minor allele, then we need to check the depth and frequency
                    elif alleleType == 'Minority':
                        depth = int(splitLine[4])
                        frequency = float(splitLine[5])
                        if depth >= 100 and frequency >= 0.025:
                            outfile.write(sample + '\t' + line)
        allAllelesFile.close()
outfile.close()
