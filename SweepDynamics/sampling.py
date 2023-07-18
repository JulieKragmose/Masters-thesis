#!/usr/bin/env python3
#Usage: python3 dataInfo.py <prefix_cds.fa>
import sys, os, random
from Bio import SeqIO
from datetime import date



#Month numbers belonging to the different seasons:
Nseasons1 = ['10', '11', '12']
Nseasons2 = ['01', '02', '03']
Sseasons = ['04', '05', '06', '07', '08', '09']



#Input argument
if len(sys.argv) == 1:
    filename = input('Please enter a name of a fasta file: ')
elif len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    print('Usage: dataInfo.py <prefix_cds.fa>')
    sys.exit(1)

#Open file
try:
    infile = open(filename, 'r')
except IOError as error:
    print("Can't read file, reason: " + str(error) + "\n")
    sys.exit(1)

prefix = filename.split('.')[0].split('_')[0]




#######################################################################################
#                                   DATA GATHERING                                    #
#######################################################################################

#Make dict with headers belonging to every year
entryCount = 0
headerDict = dict()
sequenceDict = dict()
for seqEntry in SeqIO.parse(filename, "fasta"):
    entryCount += 1

    header = ' '.join((seqEntry.id).split('_'))
    sequence = str(seqEntry.seq).upper()

    d = header.split(' | ')[-1]
    splitDate = d.split('-')
    year = splitDate[0]
    month = splitDate[1]
    
    #Sequences from October, Novembre and December belong to the following year
    if month in ['10', '11', '12']:
        year = str(int(year) + 1)
    
    #Save header with its corresponding year
    if year in headerDict:
        headerDict[year].append(header)
    else:
        headerDict[year] = [header]

    #Save sequence with its header
    if header not in sequenceDict:
        sequenceDict[header] = sequence
    else:
        print('Warning! Duplicate header: ' + header)
infile.close()


#Print stats
print(str(entryCount) + ' sequences in total\n')
print('Year\t#Sequences')
for year, list in headerDict.items():
    print(year + '\t' + str(len(list)))
print('\n')


#Find year with fewest sequences
minLength = None
for headerList in headerDict.values():
    if minLength is None or len(headerList) < minLength:
        minLength = len(headerList)
print('Sampling size: ' + str(minLength))




#######################################################################################
#                                      SAMPLING                                       #
#######################################################################################

#Extract X number of random headers for each year
samplingDict = dict()
for year, headerList in headerDict.items():
    samplingList = random.choices(headerList, k = minLength)
    samplingDict[year] = samplingList


#Check the correct number of sequences has been extracted for each year
for year, dict in samplingDict.items():
    if len(dict) != minLength:
        print('Warning: ' + year + ' has ' + str(len(dict)) + ' sequences (should be ' + str(minLength) + ')')


#Add extracted identifiers to one list while also fixing header format
sampleHeaders = []
for list in samplingDict.values():
    sampleHeaders.extend(list)



#Open files
outfileName = prefix + '_sampled_cds.fa'
try:
    infile = open(filename, 'r')
    outfile = open(outfileName, 'w')
except IOError as error:
    print("Can't read file, reason: " + str(error) + "\n")
    sys.exit(1)



#Write new file with only the sampled sequences
finalCount = 0
for year, headerList in samplingDict.items():
    for header in headerList:
        finalCount += 1
        outfile.write('>' + header + '\n')
        outfile.write(sequenceDict[header] + '\n')
infile.close()
outfile.close()



#######################################################################################
#                                  FIND ROOT SEQUENCE                                 #
#######################################################################################

#Find oldest date = root date
dates = []
headerDict = {}
for headerList in samplingDict.values():
    for header in headerList:
        d = header.split('|')[-1].strip()
        splitDate = d.split('-')
        year, month, day = int(splitDate[0]), int(splitDate[1]), int(splitDate[2])
        dates.append(date(year, month, day))

        if d in headerDict:
            headerDict[d].append(header)
        else:
            headerDict[d] = [header]
rootDate = min(dates)
print('Date of root sequence: ' + str(rootDate))


#Root header
rootHeader = headerDict[str(rootDate)][0]
print('Root sequence: ' + rootHeader)


#######################################################################################
#                                  WRITE BASH SCRIPT                                  #
#######################################################################################

identifier = rootHeader.split('|')[2].strip()

bashFilename = prefix + '_sweepDynamics.sh'
bashFile = open('../' + bashFilename, 'w')

bashFile.write("#!/bin/bash\n")
bashFile.write("CONDA_PATH=$(conda info | grep -i 'base environment' | awk '{print $4}')\n")
bashFile.write("source $CONDA_PATH/etc/profile.d/conda.sh\n")
bashFile.write("\n")
bashFile.write("#Activate conda environment\n")
bashFile.write("conda activate jukj_sweepDynamics\n")
bashFile.write("\n")
bashFile.write("#Run sweep dynamics without sampling\n")
bashFile.write('bash ../docker_files/app/SDplotsPipeline/SDplot_pipeline.sh -l true -g false -i ' + prefix + '_sampled -o ' + prefix + '_sampled -r "' + identifier + '" -n 16 -p "png"')

bashFile.close()
print('Done!')