#!/usr/bin/env python3
import sys, os

# Reconstruct consensus sequences with deletions
# Usage: indelConsensusConstruction.py <path/to/indelSummary>
# Call from .../run/human

translationDict = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W' }




#------------------------------- FUNCTIONS -----------------------------------#
def deletionData(entry):
    '''Takes a deletion entry in list form from indelSummary and saves the data as variables'''
    sample = str(entry[0])
    segment = str(entry[1])
    startPosition = int(entry[3]) + 1 #Position in indelSummary is position before deletion. We add 1 to get the first position of the deletion.
    delLength = int(entry[4])
    mutation = str(entry[-1][:-1])
    
    return sample, segment, startPosition, delLength, mutation



def translate(nucleotideSeq):
    '''Takes a nucleotide sequence and translates to protein sequence'''
    proteinSeq = ''
    for i in range(0, len(nucleotideSeq), 3):
        codon = nucleotideSeq[i:i+3]
        if len(codon) == 3:
            try:
                aa = translationDict[codon]
            except:
                aa = 'X'
        proteinSeq += aa
    
    return proteinSeq





#---------------------------------- MAIN -------------------------------------#

#Read indelSummary as input
if len(sys.argv) == 2:
    indelFilename = sys.argv[1]
elif len(sys.argv) == 1:
    indelFilename = input('Please enter name of (or path to) indelSummary file: ')
else:
    print('Error: Too many arguments.')
    print('Usage: indelConsensusConstruction.py <path/to/indelSummary>')
    sys.exit(1)


 
baseDirectory = '/srv/data/VOF/INF/JUKJ/Quasi/run/human/'


#Open indelSummary
try:
    indelFile = open(indelFilename, 'r')
except IOError as error:
    print('Cant open indelSummary file, reason: ' + str(error))
    sys.exit(1)
    

#Go thorugh every indel
for line in indelFile:
    if not line.startswith('Sample') and not line.startswith('\n'):
        splitLine = line.split('\t')
        
        indelType = splitLine[2]
        
        
        #INSERTION
        if indelType == 'ins':
            print('This script has not yet been set up to handle insertions. Skipping.')
            continue
        
        
        #DELETION
        elif indelType == 'del':
            sample, segment, startPosition, delLength, mutation = deletionData(splitLine)
            
            #Open consensus for segment in sample
            try:
                consensusFile = open(baseDirectory + sample + '/' + segment + '_consensus.fna', 'r')
            except:
                print('ERROR: Consensus file for ' + segment + ' in ' + sample + ' could not be opened or found.')
                print('Remember to call script from .../human/run')
                print('Exiting program')
                sys.exit(1)
            
            #Read consensus to string
            seq = ''
            for line in consensusFile:
                if not line.startswith('>'):
                    seq += line[:-1]
            consensusFile.close()
            
           
            #Check that the change matches the mutation in indelSummary
            deletion = '-' * delLength
            cutOut = seq[startPosition - 6 : startPosition-1] + deletion + seq[startPosition-1 + delLength : startPosition-1 + delLength + 5]
            if cutOut != mutation:
                print('ERROR: something went wrong in ' + segment + ' position ' + startPosition + ' in ' + sample) 
                print(mutation)
                print(cutOut)
                print('Skipping')
                continue
            
            #Make the deletion in the consensus
            updatedConsensus = seq[:startPosition - 1] + seq[startPosition - 1 + delLength:]
            
            #Make a length check
            if len(updatedConsensus) != len(seq) - delLength:
                print('ERROR: something went wrong when making the deletion in ' + segment + ' position ' + startPosition + ' in ' + sample)
                print('Expected length of new sequence: ' + str(len(seq) - delLength))
                print('Actual length of new sequence: ') + str(len(updatedConsensus))
                print('Skipping')
                continue
            
            #Write to file
            outfile = open(baseDirectory + '/indelSummary/' + sample + '_' + segment + '_deletion_' + str(startPosition) + '.fna', 'w')
            outfile.write('>' + sample + '_' + segment + '_pos' + str(startPosition) + '_length:' + str(delLength) + '\n')
            for i in range(0, len(updatedConsensus), 60):
                outfile.write(updatedConsensus[i:i+60] + '\n')
            outfile.close()
            
            #Translate to protein
            proteinConsensus = translate(updatedConsensus)
            
            #Write to file
            outfile = open(baseDirectory + '/indelConsensus/' + sample + '_' + segment + '_deletion_' + str(startPosition) + '.fa', 'w')
            outfile.write('>' + sample + '_' + segment + '_basePos' + str(startPosition) + '_length:' + str(delLength) + 'b\n')
            for i in range(0, len(proteinConsensus), 60):
                outfile.write(proteinConsensus[i:i+60] + '\n')
            outfile.close()
            
            
            
        #INVALID TYPE 
        else:
            print('Unknown indel type: ' + indelType + '. Skipping.')
            continue
        
indelFile.close()
    