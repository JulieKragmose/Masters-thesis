#!/usr/bin/env python3
#Make a majority consensus sequence from allAlleles.txt
import sys, os


minDepth = 50
minQual = 30



baseDir = '/srv/data/VOF/INF/JUKJ/Quasi/run/human/'
samples = ['Day1_1', 'Day1_2', 'Day3', 'Day8', 'Day14', 'Day17', 'Day21']
segments = ['PB1', 'PB2', 'PA', 'HA_H3', 'NP', 'NA_N2', 'MP', 'NS']
stopCodons = ['TAA', 'TGA', 'TAG']


translation = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
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


for sample in samples:
    for segment in segments:
        filename = 'A_' + segment + '-allAlleles.txt'
        
        #Open file
        try:
            infile = open(baseDir + sample + '/tables/' + filename, 'r')
        except:
            print(filename + ' not found in ' + sample)
            continue
        
        #Read majority alleles into dict
        alleleDict = dict()
        for line in infile:
            splitLine = line.split('\t')
            alleleType = splitLine[-1][:-1]

            
            if alleleType == 'Consensus':
                position = splitLine[1]
                allele = splitLine[2]
                depth = splitLine[4]
                frequency = splitLine[5]
                quality = splitLine[6]
                
                #Quality check
                if int(depth) >= minDepth and float(quality) >= minQual:
                    allele = splitLine[2]
                else:
                    allele = 'N'
                
                #Add to dict
                if position in alleleDict:
                    alleleDict[position].append([allele, frequency])
                else:
                    alleleDict[position] = [[allele, frequency]]
                    
        infile.close()

        
        
        #Make consensus sequence
        seq = ''
        for position, data in alleleDict.items():
            if len(data) == 1:
                seq += data[0][0]
            else:
                print('Error: Two consensus bases found for position ' + str(position) + ' in ' + segment + ' ' + sample)
                print('Exiting program')
                sys.exit(1)
            
            
        #Length check
        if len(alleleDict) != len(seq):
            print('Error: ' + segment + ' in ' + sample + ' has the wrong length. Exiting program.')
            sys.exit(1)
            
            
        #Start and stop codon check
        if seq[:3].upper() != 'ATG':
            print('Warning: Sequence for ' + segment + ' in ' + sample + ' starts with ' + seq[:3].upper() + ' instead of ATG')
        if seq[-3:].upper() not in stopCodons:
            print('Warning: Sequence does not end with a stop codon')
            


        #Write nucleotide sequence to file
        outfile = open(baseDir + sample + '/' + segment + '_consensus.fna', 'w')
        outfile.write('>' + sample + '_' + segment + '_Consensus\n')
        for i in range(0, len(seq), 60):
            outfile.write(seq[i:i+60] + '\n')
        outfile.close()
        
        
        
        #Translate to protein
        seqAA = ''
        for i in range(0,len(seq),3):
            codon = seq[i:i+3]
            if len(codon) == 3:
                if 'N' in codon:
                    AA = 'X'
                else:
                    AA = translation[codon]
                seqAA += AA
        
        #Write protein sequence to file
        outfile = open(baseDir + sample + '/' + segment + '_consensus.fa', 'w')
        outfile.write('>' + sample + '_' + segment + '_Protein_Consensus\n')
        for i in range(0, len(seqAA), 60):
            outfile.write(seqAA[i:i+60] + '\n')
        outfile.close()
        
            