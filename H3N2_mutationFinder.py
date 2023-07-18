#!/usr/bin/env python3
# My own version of the mutation finder script
# Usage H3N2_mutationFinder.py
import sys, os, glob



#############################################################################################
#                                         INITIALIZE                                        #
#############################################################################################

# Thresholds
minFrequency = 0.025
minQuality = 30
minDepth = 100
minCount = 15

# Valid input responses
negativeResponses = ['n', 'N', 'no', 'No']
positiveResponses = ['y', 'yes', 'Y', 'Yes']

# Other stuff
segments = ['PB1', 'PB2', 'PA', 'HA_H3', 'NP', 'NA_N2', 'MP', 'NS']
geneLengths = {'PB1': 2274, 'PB1-2': 273, 'PB2': 2280, 'PA': 2151, 'PA-2': 759, 'HA_H3': 1701, 'NP': 1497, 'NA_N2': 1410, 'MP': 759, 'MP-2': 294, 'NS': 693, 'NS-2': 366}

positionDict = {'PB1':[[0,2274],[94,367]] , 
                'PB2':[[0,2274]] , 
                'PA':[[0,2151],[0,573],[574,760]] , 
                'HA_H3':[[0,1701]] , 
                'NP':[[0,1497]] , 
                'NA_N2':[[0,1410]] , 
                'MP':[[0,759],[0,26,714,982]] , 
                'NS': [[0, 693],[0,30,502,838]]}

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




#############################################################################################
#                                           FUNCTIONS                                       #
#############################################################################################

def usage():
    print('\nUsage: python3 H3N2_mutationFinder.py OPTIONAL: <sample name>')
    print('Accepted responses to input request:')
    print('Negative: n, N, no, No')
    print('Positive: y, Y, yes, yes')

    print('\nExiting program')
    sys.exit(1)



def findSamples():
    '''Searches root directory for IRMA output directories and returns a list of sample directory names'''
    sampleList = []
    pathList = []
    rootDir = os.getcwd()

    # Find everything in this folder
    print('Searcing ' + str(rootDir) + ' for sample directories')
    searchList = os.listdir(".")

    # Find the entries that are existing directories
    for entry in searchList:
        path = os.path.join(rootDir, entry)
        if os.path.isdir(path):
            pathList.append(path)
            sampleList.append(entry)

    # No directories found
    if len(pathList) == 0:
        print('No sample folders found. Exiting program.')
        sys.exit(1)

    #Confirm these are the samples the user wants
    confirmed = False
    while confirmed is False:
        sampleList, confirmed = confirmSampleList(sampleList, confirmed)

    return sampleList



def confirmSampleList(sampleList, confirmed):
    '''Takes a list of samples and gives the user the option to remove samples before confirming'''
     # Print directories found
    print(str(len(sampleList)) + ' samples were found:')
    sampleNo = 0
    for sample in sampleList:
        sampleNo += 1
        print(str(sampleNo) + '\t' + sample)

    # Get user input
    proceed = input('Do you wish to proceed with these sample directories? (y/n): ')
    
    # Make no changes and proceed with these samples
    if proceed in positiveResponses:
        confirmed = True

    # Samples are not as user wants them
    elif proceed in negativeResponses:
        # Ask whether to remove directories
        remove = input('Do you wish to remove one or more of the directories? (y/n): ')       
        if remove in positiveResponses:
            unwanted = input('Please enter the directory name(s) you wish to remove (if more than one directory, separate by space): ')
            sampleList = removeDirectories(sampleList, unwanted.split(' '))
        elif remove in negativeResponses:
            print('Exiting program')
            sys.exit(0)
        else:
            print('ERROR: Input not recognized. Please see usage output for valid input responses:')
            usage()

    # Invalid response given
    else:
        print('ERROR: Input not recognized. Please see usage output for valid input responses:')
        usage()

    return sampleList, confirmed


    
def removeDirectories(sampleList, unwantedList):
    '''Takes a list of sample directory names and a list of the names to be removed'''
    updatedSampleList = [x for x in sampleList if x not in unwantedList]
    
    return updatedSampleList



def readVariantData(dataLine):
    '''Takes a line from -variants.txt and returns the needed variables'''
    dataList = dataLine.split('\t')
    position = int(dataList[1])
    depth = int(dataList[2])
    consensusBase = str(dataList[3])
    minorityBase = str(dataList[4])
    count = int(dataList[6])
    frequency = float(dataList[8])
    quality = float(dataList[10])

    return position, depth, consensusBase, minorityBase, count, frequency, quality



def validateVariant(dataLine):
    '''Takes a line of variant data and compares to thresholds defined in top of script'''
    position, depth, consensusBase, minorityBase, count, frequency, quality = readVariantData(dataLine)

    # Check if data passes thresholds:
    if quality >= minQuality and depth >= minDepth and frequency >= minFrequency and count >= minCount:
        return 'PASS'
    else:
        return 'FAIL'



def checkConsensusLength(segment, sequence, warning_list):
    '''Check if length of majority consensus matches "official" length'''
    diff = len(sequence) - geneLengths[segment]   
    if diff < 0:
        warning_list.append('Majority sequence is ' + str(abs(diff)) + ' bases shorter than reference')
    elif diff > 0:
        warning_list.append('Majority sequence is ' + str(abs(diff)) + ' bases longer than reference')  

    return warning_list



def makeConsensus(sample, segment, warningList):
    '''Reads consensus sequence to string'''
    # Open consensus file
    try:
        consensusFile = open(sample + '/' + segment + '_consensus.fna', 'r')
    except OSError as error:
        print('File {:s} failed to open'.format(error.filename))
        sys.exit(1)

    #Read whole consensus to string
    consensusSeq = ''
    for line in consensusFile:
        if not line.startswith('>'):
            consensusSeq += line[:-1]
    consensusFile.close()

    return consensusSeq, warningList



def mutationInGene(genePositions, minorityPosition):
    '''Checks if mutation is in the gene'''

    if len(genePositions) == 2:
        if minorityPosition > genePositions[0] and minorityPosition <= genePositions[1]:
            return True
        else:
            return False
    
    elif len(genePositions) == 4:
        if (minorityPosition > genePositions[0] and minorityPosition <= genePositions[1]) or (minorityPosition > genePositions[2] and minorityPosition <= genePositions[3]):
            return True
        else:
            return False

    #Invalid number of gene positions
    else:
        print('ERROR: Invalid number of gene positions')
        print('Exiting program')
        sys.exit(1)



def compareLength(consensusSequence, minoritySequence):
    '''Check if majority and minority sequences have same length. I not, program will will stop as something must be wrong with the code'''
    if len(consensusSequence) != len(minoritySequence):
        print('ERROR: Majority and minority sequences do not have the same length')
        print('Major (' + str(len(consensusSequence)) + '):')
        print(consensusSequence)
        print('Minor (' + str(len(minoritySequence)) + '):')
        print(minoritySequence)
        print('Exiting program')
        sys.exit(1)



def checkLength(geneName, sequence, warningList):
    '''Check if length of majority consensus matches "official" length'''
    diff = len(sequence) - geneLengths[geneName]   
    if diff < 0:
        warningList.append('Majority gene sequence is ' + str(abs(diff)) + ' bases too short')
    elif diff > 0:
        warningList.append('Majority gene sequence is ' + str(abs(diff)) + ' bases too long')  

    return warningList



def translateSequence(sequence):
    '''Translate the spliced genes to amino acid to make sure it looks alright'''
    proteinSeq = ''
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
        if len(codon) == 3:
            try:
                aa = translationDict[codon]
            except:
                aa = 'X'
        proteinSeq += aa

    return proteinSeq



def checkStartCodon(sequence, warningList):
    '''Check if majority consensus starts with ATG'''
    if not sequence.startswith('ATG'):
        warningList.append('Consensus sequence starts with ' + sequence[:3] + ' instead of ATG')

    return warningList



def checkStopCodon(proteinSeq, warningList):
    '''Check if protein sequence ends with a stop codon'''
    if proteinSeq[-1] != '*':
        warningList.append('Gene does not end in stop codon. Check if protein sequence looks valid:' + '\n' + proteinSeq) 

    return warningList



def findAA(position, consensus):
    '''Takes the position of a variant site (int) and a consensus sequence (str) and returns the amino acid for the corresponding codon'''
    #Find starting position of the codon which the minority variant belongs to
    modulo = position%3
    if modulo == 0: #If pos is the 3rd base in a codon
        codon_start = position - 2
    elif modulo == 1: #If pos is 1st base in a codon
        codon_start = position
    elif modulo == 2: #If pos i 2nd base in a codon
        codon_start = position - 1

    codon = consensus[codon_start - 1 : codon_start + 2]
    #print('Codon: ' + codon)
   
    if 'N' in codon:
        return codon, 'X'
    else:
        aa = translationDict[codon]
        return codon, aa



def insertVariant(consensusSequence, variantData):
    '''Inserts variant into consensus sequence to create minority sequence'''
    dataList = readVariantData(variantData)
    position = dataList[0]
    base = dataList[3]
    minoritySequence = consensusSeq[:position-1] + base + consensusSeq[position:]

    return minoritySequence



#MAYBE NOT NEEDED
def spliceGene(sequence, genePositions):
    splicedSequence = sequence[genePositions[0]:genePositions[1]] + sequence[genePositions[2]:genePositions[3]]

    return splicedSequence



def updateVariantPosition(position, genePositions):
    '''Updates position of variant if placed in the second part of spliced gene'''

    #If gene is spliced and variant is placed in second part of spliced gene
    if len(genePositions) == 4 and genePositions[2] < position <= genePositions[3]: 
        position = position - (genePositions[2] - genePositions[1])
    #If gene is not spliced, but startposition is not 0
    elif len(genePositions) == 2 and genePositions[0] != 0:
        position = position - genePositions[0] + 1

    return position



def errorCheckSequences(geneName, minorityGene, consensusGene, warningList):
    '''Checks the majority and minority genes for issues'''
    #Make sure that minority and consensus sequence have the same length
    compareLength(consensusGene, minorityGene)

    #Check that the gene has the length it is supposed to have
    warningList = checkLength(geneName, consensusGene, warningList)

    #Start codon
    warningList = checkStartCodon(consensusGene, warningList)

    #Stop codon
    proteinGene = translateSequence(minorityGene)
    warningList = checkStopCodon(proteinGene, warningList)

    return warningList



def printWarnings(warning_list):
    '''Print warnings found for segment'''
    print(str(len(warning_list)) + ' warning(s)')
    for warning in warning_list:
        print(warning)
    print('\n\n')



def printSample(sample):
    head = '#     ' + str(sample) + '     #'
    line = '#' * len(head)
    print(line)
    print(head)
    print(line)



def printGene(segment, geneNumber):
    print('\n------------------------------------')
    if geneNumber == 1:   
        print('               ' + segment)
    elif geneNumber == 2:
        print('              ' + segment + ' prot2')



def getGene(sequence, genePositions):
    #Unspliced gene
    if len(genePositions) == 2:
        gene = sequence[genePositions[0]:genePositions[1]]
    #Spliced gene       
    elif len(genePositions) == 4:
        gene = sequence[genePositions[0]:genePositions[1]] + sequence[genePositions[2]:genePositions[3]]   

    return gene



def initializeOutfile(sample, segment, geneNumber):
    filename = None
    if geneNumber == 1:
        filename = str(sample) + '/mutations/' + str(sample) + '_mutations_' + str(segment) + '.txt'
    else:
        filename = str(sample) + '/mutations/' + str(sample) + '_mutations_' + str(segment) + '_prot2.txt'

    outfile = open(filename, 'w')
    outfile.write('\t'.join(["Sample", "Segment", "Position", "MajorCodon", "MajorAA", "MinorCodon", "MinorAA", "Mutation", '\n']))
    outfile.close()

    return filename



def writeToOutfile(filename, variantData, segment, sample, majorCodon, majorAA, minorCodon, minorAA):
    position = readVariantData(variantData)[0]
    outfile = open(filename, 'a')
    outfile.write(sample + '\t' + segment + '\t' + str(position) + '\t' + majorCodon + '\t' + majorAA + '\t' + minorCodon + '\t' + minorAA + '\t')

    if majorAA == minorAA:
        outfile.write('Synonymous\n')
    else:
        outfile.write('Non-synonymous\n')
    outfile.close()












#############################################################################################
#                                              MAIN                                         #
#############################################################################################


# No specific sample given; search for sample directories
if len(sys.argv) == 1:
    samples = findSamples()
# One specific sample given
elif len(sys.argv) == 2:
    samples = [str(sys.argv[1])]
# Error, print usage and exit
else:
    usage()



# SAMPLE
for sample in samples:   
    printSample(sample)

    #Make minority directory if it doesn't already exist
    mutDirExist = os.path.exists('./' + sample + '/mutations')
    if not mutDirExist:
        os.makedirs('./' + sample + '/mutations')


    # SEGMENT
    for segment in segments:

        # Open variant file
        try:
            variantFile = open(sample + '/tables/A_' + segment + '-variants.txt', 'r')
        except OSError as error:
            print('File {:s} failed to open'.format(error.filename))
            continue

        #Check if there are any variants in this segment
        variantEntries = []
        for line in variantFile:
            if not line.startswith('Reference_Name'):
                variantEntries.append(line)
        variantFile.close()
        if len(variantEntries) > 0:


            # GENES IN SEGMENT
            geneNumber = 0
            for genePosition in positionDict[segment]:
                warningList = []
                geneNumber += 1
                printGene(segment, geneNumber)

                # Read consensus to string
                consensusSeq, warningList = makeConsensus(sample, segment, warningList)
  

                # VARIANTS IN SEGMENT
                variantFound = False
                for variantData in variantEntries:
                    position = int(variantData.split('\t')[1]) 
                    #print('Variant position: ' + str(position))

                    # Check that the mutation is within this gene and passes threholds
                    if mutationInGene(genePosition, position) is True and validateVariant(variantData) == 'PASS':

                        #If this is the first variant found in this gene, initialize the outfile
                        if variantFound is False:
                            filename = initializeOutfile(sample, segment, geneNumber)
                            variantFound = True

                        #Insert the mutation 
                        minoritySeq = insertVariant(consensusSeq, variantData)
                        #print('Length of full sequenses: ' + str(len(consensusSeq)) + ' v. ' + str(len(minoritySeq)))
                    
                        # Update position of variant site
                        position = updateVariantPosition(position, genePosition)

                        # Extract the gene from the created sequences
                        majorityGene = getGene(consensusSeq, genePosition)
                        minorityGene = getGene(minoritySeq, genePosition)

                        #Name of gene
                        if geneNumber == 1:
                            geneName = segment
                        elif geneNumber == 2:
                            geneName = segment + '-2'

                        # Check the genes for errors
                        warningList = errorCheckSequences(geneName, minorityGene, majorityGene, warningList)
                        
                        # Find codon and amino acid for consensus and minority sequence
                        majCodon, majAA = findAA(position, majorityGene)
                        minCodon, minAA = findAA(position, minorityGene)

                        #Write to outfile
                        writeToOutfile(filename, variantData, segment, sample, majCodon, majAA, minCodon, minAA)
                        
                printWarnings(warningList)
    
                        
                

                    









        
        

