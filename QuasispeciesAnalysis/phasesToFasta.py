 # !/usr/bin/env python3
 # Turn minority variants into fasta files, while taking linked variants into account

 # Usage: python3 phasesToFasta.py
 # Remember to call the script from '.../run/human'
import sys, os, glob


#Inititalize
depthThreshold = 50
countThreshold = 15
freqThreshold = 0.025
qualThreshold = 30

segments = ['PB1', 'PB2', 'PA', 'HA_H3', 'NP', 'NA_N2', 'MP', 'NS']
lengthDict = {'PB1': 2341, 'PB2': 2341, 'PA': 2233, 'HA_H3': 1778, 'NP': 1565, 'NA_N2': 1413, 'MP': 1027, 'NS': 890}
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



#Find sample directories
sampleList = []
pathList = []
rootDir = os.getcwd()
searchList = os.listdir(".")
for entry in searchList:
    path = os.path.join(rootDir, entry)
    if os.path.isdir(path):
        pathList.append(path)
        sampleList.append(entry)

if len(pathList) == 0:
    print('No sample folder found. Exiting program.')
    sys.exit(1)

#Print samples found and aks whether to proceed
print(str(len(sampleList)) + ' samples were found:')
sampleNo = 0
for sample in sampleList:
    sampleNo += 1
    print(str(sampleNo) + '\t' + sample)



proceed = input('Do you wish to proceed? (y/n): ')



#Exit
if proceed == 'n':
    print('Exiting program')
    sys.exit(0)
elif proceed != 'n' and proceed != 'y':
    print('Invalid input. Enter "y" for yes or "n" for no')
    sys.exit(1)



#Run program
else:     
    # Go through every sample
    for path in pathList:
        sample = path.split('/')[-1]
        print('\n-------------------------------------')
        print(sample)

        #Make minority directory if it doesn't already exist
        minDirExist = os.path.exists(path + '/minority')
        if not minDirExist:
            os.makedirs(path + '/minority')

        #Go through every segment in this sample
        for segment in segments:
            phaseDict = dict()

            # Open variant file
            try:
                variantFile = open(path + '/tables/A_' + segment + '-variants.txt', 'r')
            except IOError as error:
                print('Variant file for ' + segment + ' was not found or could not be opened')
                continue
            logFile = open(path + '/minority/' + sample + '_' + segment + '_phasesToFasta.log', 'w')

            
            # Go through each minor variant
            for line in variantFile:
                splitLine = line[:-1].split('\t')
                
                # Find valid minor variants and save in dict (phased variants are saved under same key)
                if splitLine[0] != 'Reference_Name': #Skip header
                    mutation = splitLine[3] + splitLine[1] + splitLine[4]
                    depth = splitLine[2]
                    count = splitLine[6]
                    frequency = splitLine[8]
                    quality = splitLine[10]
                    phaseGroup = splitLine[-1]

                    logFile.write('\n------------------------------ ' + mutation + ' ------------------------------\n')
                    logFile.write('\nMutation\tDepth\tFrequency\tphaseGroup\tStatus\n')

                    # If minor variant is valid; save it
                    if int(depth) >= depthThreshold and int(count) >= countThreshold and float(frequency) >= freqThreshold and float(quality) >= qualThreshold:
                        logFile.write(mutation + '\t' + depth + '\t' + count + '\t' + frequency + '\t' + phaseGroup + '\t' + 'PASSED\n')
                        if phaseGroup in phaseDict:
                            phaseDict[phaseGroup].append(mutation)
                        else:
                            phaseDict[phaseGroup] = [mutation]
                    # If minor variant is invalid
                    else:
                        logFile.write(mutation + '\t' + depth + '\t' + count + '\t' + frequency + '\t' + phaseGroup + '\t' + 'FAILED\n')
            variantFile.close()
            logFile.write('\n')

                    
            # Stats on minor variants
            if len(phaseDict) > 0:
                logFile.write(str(len(phaseDict)) + ' phase group(s) found\n')
                logFile.write('\n')
            # If no valid minor variants found, remove log file and continue to next segment
            else:
                #print('No valid minor variants found for ' + segment)
                logFile.close()
                os.remove(path + '/minority/' + sample + '_' + segment + '_phasesToFasta.log')
                continue
            
            
            # read consensus (made with consensusFromAllAlleles.py) and save as string
            consensusFile = open(sample + '/' + segment + '_consensus.fna', 'r')
            consensusSeq = ''
            for line in consensusFile:
                if not line.startswith('>'):
                    consensusSeq += line[:-1]
            consensusFile.close()

            
            # Check valid start of consensus
            if consensusSeq.startswith('ATG'):
                logFile.write('Consensus sequence starts with valid codon ATG\n')   
            else:
                print("Warning: Consensus sequence for " + segment + " starts with " + consensusSeq[:3] + " instead of ATG.")
                logFile.write("Warning: Consensus sequence for " + segment + " starts with " + consensusSeq[:3] + " instead of ATG\n")


            # Length check
            if len(consensusSeq) != lengthDict[segment]:
                print('Warning: consensus sequence(=' + str(len(consensusSeq)) + ') and official length of segment (' + str(lengthDict[segment]) + ') for ' + segment + ' in ' + sample + ' differ by ' + str(abs(lengthDict[segment] - len(consensusSeq))) + ' bases')
                logFile.write('Official length of segment: ' + str(lengthDict[segment]))
                logFile.write('Consensus sequence length: ' + str(len(consensusSeq)) + '\n') 
                logFile.write('Warning: consensus sequence(=' + str(len(consensusSeq)) + ') and official length of segment (' + str(lengthDict[segment]) + ') for ' + segment + ' in ' + sample + ' differ by ' + str(abs(lengthDict[segment] - len(consensusSeq))) + ' bases')


            logFile.write('Proceeding to make fasta files for minor variants\n')
            logFile.write('\n')
            

            # Change bases at variant positions
            for phaseGroup, mutationList in phaseDict.items():
                logFile.write('Phase group ' + phaseGroup + ' contains ' + str(len(mutationList)) + ' mutation(s)\n')
                minorityFasta = open(path + '/minority/A_' + segment + '_phase' + phaseGroup + '.fna', 'w')

                minoritySeq = consensusSeq
                for mut in mutationList:
                    logFile.write(mut + '\n')
                    consensusBase = mut[0]
                    minorityBase = mut[-1]
                    position = int(mut[1:-1])

                    # Replace base at the variant position with minority
                    minoritySeq = ''.join([minoritySeq[:position- 1], minorityBase, minoritySeq[position:]])
                    logFile.write(consensusBase + ' was replaced with ' + minoritySeq[position - 1] + ' at position ' + str(position) + '\n')

                    # Length check
                    if len(minoritySeq) != len(consensusSeq):
                        print('Run for ' + segment + ' stopped due to error. See log file.')
                        logFile.write('Something went wrong when replacing ' + consensusBase + ' with ' + minorityBase + ' at position ' + str(position) + '. Consensus sequence and minority sequence are not the same length.')
                        sys.exit(1)
                logFile.write('\n')

                
                # Write variant nucleotide sequence to fasta file
                # Header
                minorityFasta.write('>A_' + segment + '_phase' + phaseGroup)
                for m in mutationList:
                    minorityFasta.write('_' + m)
                minorityFasta.write('\n')

                # Sequence
                for j in range(0, len(minoritySeq), 60):
                    minorityFasta.write(minoritySeq[j:j+60] + '\n')
                minorityFasta.close()

                
                # Translate into protein
                minoritySeqProtein = ''
                for k in range(0, len(minoritySeq), 3):
                    codon = minoritySeq[k:k+3]

                    if len(codon) == 3:
                        try:
                            aa = translationDict[codon]
                        except:
                            aa = 'X'
                    minoritySeqProtein += aa

                # Write protein sequence to fasta
                minorityFastaProtein = open(path + '/minority/A_' + segment + '_phase' + phaseGroup + '.fa', 'w')
                minorityFastaProtein.write('>A_' + segment + '_phase' + phaseGroup + '\n')
                for l in range(0, len(minoritySeqProtein), 60):
                    minorityFastaProtein.write(minoritySeqProtein[l:l+60] + '\n')
                minorityFastaProtein.close()


    logFile.close()
    print('\nphasesToFasta.py ran succesfully')
