#!/usr/bin/env python3
import sys, os
from Bio import SeqIO


#------------------------------ FUNCTIONS -------------------------------#

def makeHeader(header):
    # Slit into variables
    splitHeader = header.split('|')
    segment = segmentNo[splitHeader[-1]]
    identifier = splitHeader[2]
    date = splitHeader[3]
    subtype = ''.join(splitHeader[5].split('_'))
    number = identifier.split('/')[2]

    #Make new header in correct format
    newHeader = '>' + '|'.join([number, segment, identifier, subtype, date])
    return number, newHeader, segment



def deduplicate(sequence, header, number, sequenceDict, duplicateDict):
    if sequence not in sequenceDict:
        sequenceDict[sequence] = header
        duplicateDict[sequence] = [number]
    else:
        duplicateDict[sequence].append(number)
    
    return sequenceDict, duplicateDict



def writeFasta(filename, sequenceDict):
    outfile = open(filename, 'w')
    for sequence, header in sequenceDict.items():
        outfile.write(header + '\n')
        outfile.write(sequence + '\n')
    outfile.close()





#--------------------------------- MAIN ---------------------------------#

baseDir = '/srv/data/VOF/INF/JUKJ/Quasi/phylo/data/'


#Thresholds
perN = 50


#Initialize
segmentNo = {'4': 'HA', '6': 'NA'}
HAsequences = dict()
HAduplicates = dict()
NAsequences = dict()
NAduplicates = dict()
HAentryCount = 0
NAentryCount = 0
removedCount = 0

fastaFile = baseDir + 'gisaid_epiflu_sequence.fasta'
for seqEntry in SeqIO.parse(fastaFile, "fasta"):

    # Read entry data
    oldHeader = seqEntry.id
    sequence = str(seqEntry.seq).upper()

    # Make new header in correct format
    seqID, newHeader, segment = makeHeader(oldHeader)


    
    # Check if there are too many Ns in the squence
    if (float(sequence.count("N")) / float(len(sequence)) * 100) <= perN:
        # It's good to go, save it do dicts
        if segment == 'HA':
            HAentryCount += 1
            HAsequences, HAduplicates = deduplicate(sequence, newHeader, seqID, HAsequences, HAduplicates)
        elif segment == 'NA':
            NAentryCount += 1
            NAsequences, NAduplicates = deduplicate(sequence, newHeader, seqID, NAsequences, NAduplicates)
    else:
        removedCount += 1
    """

"""
# Print stats
print(str(removedCount) + ' sequences were removed due to too many Ns\n')
print('Of the remaining sequences, there were:')

print(str(HAentryCount) + ' HA sequences in total')
print(str(len(HAsequences)) + ' unique HA sequences')

print(str(NAentryCount) + ' NA sequences in total')
print(str(len(NAsequences)) + ' unique NA sequences')


# Write unique sequences to outfile
writeFasta(baseDir + 'HA-22-23.fna', HAsequences)
writeFasta(baseDir + 'NA-22-23.fna', NAsequences)
