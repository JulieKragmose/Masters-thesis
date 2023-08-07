#!/usr/bin/env python3
#Check every indel file and determine which are real
import sys, os


#Initialize
minDepth = 50
minQuality = 30
minFrequency = 0.025
minCount = 15
samples = ['Day1_1', 'Day1_2', 'Day3', 'Day8', 'Day14', 'Day17', 'Day21']
segments = ['PB1', 'PB2', 'PA', 'HA_H3', 'NP', 'NA_N2', 'MP', 'NS']

outfile = open('/srv/data/VOF/INF/JUKJ/Quasi/run/human/indelSummary.txt','w')
outfile.write('Sample\tSegment\tType\tPosition\tLength\tCount\tDepth\tFrequency\tAverageQuality\tConfidenceNotMacError\tMutation\n')


for sample in samples:
    baseDir = '/srv/data/VOF/INF/JUKJ/Quasi/run/human/' + sample + '/tables/'
    for segment in segments:

        ##### INSERTIONS #####
        try:
            insFile = open(baseDir + 'A_' + segment + '-insertions.txt', 'r')
        except IOError as error:
            print('Warning: Cant open file, reason: ', str(error))
            continue
        
        
        for line in insFile:
            if not line.startswith('Reference_Name'):
                splitLine = line.split('\t')
                
                #Check if IRMA says its real
                status = splitLine[2]

                if status == 'TRUE':  
                    depth = splitLine[6]
                    frequency = splitLine[7]
                    quality = splitLine[8]
                    count = splitLine[5]
                    
                    #Do we think it's real?
                    if int(depth) >= minDepth and float(frequency) >= minFrequency and float(quality) >= minQuality and int(count) >= minCount:
                        position = splitLine[1]
                        mutation = splitLine[2]
                        confidence = splitLine[9]
                        length = len(mutation)
                        
                        outfile.write(sample + '\t' + segment + '\t' + 'ins' + '\t' + position + '\t' + str(length) + '\t' + count + '\t' + depth + '\t' + frequency + '\t' + quality + '\t' + confidence + '\t' + mutation + '\n')
        insFile.close()
        
        
        ##### DELETIONS #####
        try:
            delFile = open(baseDir + 'A_' + segment + '-deletions.txt', 'r')
        except IOError as error:
            print('Warning: Cant open file, reason: ', str(error))
            continue

        for line in delFile:
            if not line.startswith('Reference_Name'):
                splitLine = line.split('\t')
                
                #Check if IRMA says its real
                status = splitLine[4]
                
                if status == 'TRUE':
                    depth = splitLine[6]
                    frequency = splitLine[7]
                    count = splitLine[5]
                    
                    #Do we think it is real?
                    if int(depth) >= minDepth and float(frequency) >= minFrequency and int(count) >= minCount:
                        position = splitLine[1]
                        length = splitLine[2]
                        mutation = splitLine[3]
                        
                        outfile.write(sample + '\t' + segment + '\t' + 'del' + '\t' + position + '\t' + length + '\t' + count + '\t' + depth + '\t' + frequency + '\tNA\tNA\t' + mutation + '\n')
                        

        delFile.close()   

    outfile.write('\n')    
outfile.close()
            
