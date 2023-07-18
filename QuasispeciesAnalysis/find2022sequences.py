#!/usr/bin/env python3
#Find the gisaid sequences that are from 4 months before and 4 months after the last sample was taken
#(August 2022 - March 2023)
import sys, os
from Bio import SeqIO


#Input
if len(sys.argv) == 2:
    filename = sys.argv[1]
elif len(sys.argv) == 1:
    filename = input('Please write name of fasta file: ')
else:
    print("Usage: python3 find2022sequences.py <fasta file>")
    sys.exit(1)


#Initialize
late2022 = ['08','09','10','11','12']
early2023 = ['01','02','03']
patientDates = ['2022-11-08', '2022-11-10', '2022-11-15', '2022-11-21', '2022-11-24', '2022-11-28']
seenIdentifier = set()


outfilename = filename.split('.')[0] + '_filtered.fna'
outfile = open(outfilename, 'w')
for seqEntry in SeqIO.parse(filename, "fasta"):
        header = seqEntry.id
        splitHeader = header.split('|')
        date = splitHeader[-1].strip()
        dateList = date.split('-')

        year, month, day = dateList[0], dateList[1], dateList[2]
        
        #Samples from late 2022 and early 2023
        if (year == '2022' and month in late2022) or (year == '2023' and month in early2023):
            sequence = str(seqEntry.seq).upper()
            segment = splitHeader[1]
            identifier = splitHeader[2]
            ID = identifier.split('/')[2]

            #Don't include sequences with identifiers already written to outfile
            if identifier not in seenIdentifier and date not in patientDates:
                newHeader = '_'.join(['_'.join([year[-2:],month,day]), ID, segment, 'DK'])
                outfile.write('>' + newHeader + '\n')
                outfile.write(sequence + '\n')
                seenIdentifier.add(identifier)

outfile.close()  
