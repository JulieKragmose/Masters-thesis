#!/usr/bin/env python3
# Run this program on IRMA output folders to detect potential contamination/coinfection
import sys, os, glob, math




# DEFINE THRESHOLDS HERE:
maxVariantPercentage = 0.0028       # Percentage of segment length expected to have variants
maxFreq = 0.02                      # Variants with frequency higher than this are considered high frequency variants
maxFreqProp = 0.10                  # Max proportion of high frequency variants out of all variants allowed
maxSecondarySize = 0.10             # Max size of secondary data compared to primary data
maxFlags = 2                        # Max number of flags raised before the whole sample is flagged as potentially contaminated




##############################################################################
#                                INITIALIZE                                  #
##############################################################################

#root = "C:/Users/julie/OneDrive/Skrivebord/Speciale/Quasispecies/"
root = os.getcwd()

segmentLengths = {'PB2': 2341, 'PB1': 2341, 'PA': 2233, 'HA_H3': 1778, 'HA_H1':  1778, 'NP': 1565, 'NA_N2': 1413, 'NA_N1': 1413, 'MP': 1027, 'NS': 890}

knownSubtypes = ['H1N1', 'H3N2', 'H1N2']




##############################################################################
#                                 FUNCTIONS                                  #
##############################################################################
def usage():
    print('\nUsage: python3 contaminationAnalysis.py')
    print('Your current working directory has to be where the IRMA output folders are (fx. .../run/human)')
    sys.exit(1)


def findDirectories(root):
    '''Returns a list of all directories found in the root directory'''
    dirList = []
    searchList = os.listdir(root)
    for entry in searchList:
        path = os.path.join(root, entry)
        if os.path.isdir(path):
            dirList.append(entry)
    return dirList



def primaryData(samplePath):
    '''Find primary genus and subtype for sample'''
    genera = set()

    #Look at the READ_COUNTS file
    readCountFile = openFile(samplePath + 'tables/READ_COUNTS.txt', 'r')
    for line in readCountFile:
        record = line.split('\t')[0]
        dataGroup = record.split('-')[0]
        
        #Find primary data
        if dataGroup == '4':
            genus = record.split('_')[0].split('-')[-1]
            genera.add(genus)

            #Find primary HA and NA subtype
            segment = '_'.join(record.split('_')[1:])
            if segment[:2] == 'HA':
                HAsubtype = segment
            elif segment[:2] == 'NA':
                NAsubtype = segment
    readCountFile.close()
    

    
    #One primary genus found 
    if len(genera) == 1:
        genus = ''.join(genera)
        subtype = HAsubtype.split('_')[1] + NAsubtype.split('_')[1]
        segments = ['PB2', 'PB1', 'PA', HAsubtype, 'NP', NAsubtype, 'MP', 'NS']
        
        #Write to sample report
        sampleReport = openFile(samplePath + sample + '_contaminationReport.txt', 'a')
        sampleReport.write('Primary data is influenza ' + genus + ' (' + subtype + ')' + '\n')
        if subtype not in knownSubtypes:
            sampleReport.write('Warning: ' + subtype + ' is an unknown subtype!\n')
        sampleReport.write('\n')
        sampleReport.close()

        return genus, segments
        
    
    #Multiple primary genera found (I don't know if it's possible but I included it anyway)
    elif len(genera) > 0:
        print('Warning: Multiple influenza genera found in primary data: ' + ' and '.join(genera) + '\n')
        print("I don't know how to do the analysis with two primary genera. Exiting program")
        usage()
    
    #No primary genus found; something is wrong
    else:
        print('Error: No genus found in primary data for ' + sample)
        print('Check if something is wrong with READ_COUNTS.txt')
        print('Exiting program')
        usage()



def openFile(filename, action):
    '''Try to open a file. Action is either "r", "w" or "a" (read, write or append)'''
    try:
        return open(filename, str(action))
    except IOError as error:
        print('Error opening file: ' + str(error))
        sys.exit(1) #Maybe it should just throw a warning and continue??



def initContaminationReport(root):
    '''Create contamination report file and write header'''
    contaminationReport = open(root + '/contaminationReport.txt', 'w')
    contaminationReport.write("{:<10} {:<10} {:^3} {:^10} {:^10} {:^10} {:^10} {:^10}".format('Sample', 'Flagged', '|', '(1)Count', '(2)Freq', '(3)Prop', '(4)Reads', '(5)Assembly') + '\n')
    contaminationReport.close()



def initSampleReport(samplePath):
    '''Create empty sample report file'''
    sampleReport = open(samplePath + sample + '_contaminationReport.txt', 'w')
    sampleReport.close()



def printSample(sample):
    '''Print name of sample to screen'''
    nameLength = len(sample)
    print('\n')
    print('#' * (nameLength + 20))
    print('#         ' + sample + '         #')
    print('#' * (nameLength + 20))
            

    
def summarizeSecondary(samplePath, primaryGenus):
    '''Writes a summary of secondary data to sample report'''
    readCountFile = openFile(samplePath + 'tables/READ_COUNTS.txt', 'r')
        
    # Initialize
    genusDict = dict()
    subtypeDict = dict()
    
    for line in readCountFile:
        splitLine = line.split('\t')
        recordSplit = splitLine[0].split('-')
        group = recordSplit[0]            
        
        # Secondary data (group 5)
        if group == '5':
            genus = recordSplit[1].split('_')[0]
            segment = recordSplit[1].split('_')[1]
            patternCount = splitLine[2] # Pattern count       
        
            # Secondary subtype
            if genus == primaryGenus:
                subtype = recordSplit[1].split('_')[-1]
                subtypeDict[subtype] = patternCount                    
            
            # Secondary genus
            else:
                if genus in genusDict:
                    genusDict[genus][segment] = patternCount
                else:
                    genusDict[genus] = {segment: patternCount}
    readCountFile.close()
    
    #Subtype read counts
    sampleReport = openFile(samplePath + sample + '_contaminationReport.txt', 'a')
    sampleReport.write('-------Subtypes of influenza ' + primaryGenus + '-------\n')
    if len(subtypeDict) > 0:
        sampleReport.write('\n')
        sampleReport.write("{:<8} {:<8}".format('Segment', 'PatternCount') + '\n')
        for subtype, readCount in subtypeDict.items():
            sampleReport.write("{:<8} {:<8}".format(subtype, readCount) + '\n') 
        sampleReport.write('\n\n')
    else:
        sampleReport.write('None\n\n')    
                            
    #Secondary genera read counts
    sampleReport.write('-----------Secondary genera----------\n\n')
    if len(genusDict) > 0:
        for genus,segmentDict in genusDict.items():
            if genus != primaryGenus:
                sampleReport.write("{:<8} {:<8} {:<8}".format('Genus', 'Segment', 'ReadCount') + '\n')
                for segment, readCount in segmentDict.items():
                    sampleReport.write("{:<8} {:<8} {:<8}".format(genus, segment, readCount) + '\n')
        sampleReport.write('\n')
    else:
        sampleReport.write('None\n\n')
    sampleReport.close()



def isSampleFolder(samplePath):
    '''Checking if a given folder is IRMA output by seeing if the directories we need exist'''
    tables = os.path.exists(samplePath + 'tables')
    secondary = os.path.exists(samplePath + 'secondary')

    #Do the folders exist? If yes, then it is an IRMA sample folder
    if tables is True and secondary is True:
        return True
    else:
        return False


#----Flag functions----#
#FLAG1
def variantCount(samplePath, segments):
    '''Counts how many minority variants are found in each segment of a given sample'''
    #Find variant files   
    variantFileList = glob.glob(samplePath + 'tables/' + primaryGenus + '*-variants.txt')
    
    if len(variantFileList) > 0:
        countDict = dict()       
        for variantFilePath in variantFileList:    
            segment = '_'.join(variantFilePath.split('/')[-1].split('-')[0].split('_')[1:])
            
            #Open variant file
            variantFile = openFile(variantFilePath, 'r')
                       
            # Count number of variants
            for line in variantFile:
                if not line.startswith('Reference_Name') and not line.strip() == '': #Skip header and empty lines
                    
                    if segment in countDict:
                        countDict[segment] += 1
                    else:
                        countDict[segment] = 1  
                        
            # Add segment with count of 0, if there were no variants in file
            if segment not in countDict:
                countDict[segment] = 0  
            variantFile.close()
            
        
    
    # Check if any segments have too many variants
    flaggedForCount = []
    for segment, count in countDict.items():
        maxCount = math.floor(segmentLengths[segment] * maxVariantPercentage)    #Define how many variants this segment is allowed to have
        if count > maxCount:
            flaggedForCount.append(segment)  

    # Write results to sample report 
    writeCount(samplePath, flaggedForCount, countDict, segments)

    # If at least one segment has too many variants; return True
    if len(flaggedForCount) > 0:
        return True, flaggedForCount, countDict
    else:
        return False, flaggedForCount, countDict
    
    
#FLAG2
def variantFrequency(samplePath, flaggedForCount, countDict, segments):
    '''Finds frequencies of mutations for segments that are flagged as having too many minority variants'''
    freqDict = dict()
    for segment in countDict.keys():
        
        # Open variant file
        variantFile = openFile(samplePath + 'tables/' + primaryGenus + '_' + segment + '-variants.txt', 'r')
        
        # # Find frequencies of all mutations in all segments
        for line in variantFile:
            if not line.startswith('Reference_Name') and not line.strip() == '': #Skip header and empty lines
                frequency = float(line.split('\t')[8])
                if segment in freqDict:
                    freqDict[segment].append(frequency)     
                else:
                    freqDict[segment] = [frequency]
        variantFile.close()

        # How many mutations are above the frequency threshold?
        flaggedForFreq = set()
        highFreqCountDict = {}
        for segment, freqList in freqDict.items():
            for f in freqList:

                # Count frequencies above threshold
                if f > maxFreq:
                    if segment in highFreqCountDict:
                        highFreqCountDict[segment] += 1
                    else:
                        highFreqCountDict[segment] = 1

                    # If this segment was previously flagged for high count of variants; flag for high frequency
                    if segment in flaggedForCount:
                        flaggedForFreq.add(segment)

    # Write results to sample report 
    writeFreq(samplePath, highFreqCountDict, segments)
    
    # If any segments were flagged for high frequency; return True
    if len(flaggedForFreq) > 0:
        return True, flaggedForFreq, freqDict, highFreqCountDict
    else:
        return False, flaggedForFreq, freqDict, highFreqCountDict
   
    
#FLAG3
def highFreqProportion(samplePath, frequencyDict, highFreqCountDict, flaggedForFreq, segments):  
    '''Calculates proportion of high frequency variants in segments previously flagged'''
    propDict = {}
    for segment, freqList in frequencyDict.items():   
        # How many variants are in this segment?
        variantCount = len(freqList)
         
        # Does the proportion exceed the threshold? (Only in flagged segments!)
        if segment in flaggedForFreq:
            proportion = round(highFreqCountDict[segment] / variantCount, 2)
            if proportion > maxFreqProp:
                propDict[segment] = proportion 
    flaggedForProp = list(propDict.keys())

    # Write results to sample report
    writeProp(samplePath, propDict, segments)

    # If any segments were flagged for high proportion; return True
    if len(flaggedForProp) > 0:
        return True, flaggedForProp, highFreqCountDict, propDict
    else:
        return False, flaggedForProp, highFreqCountDict, propDict


#FLAG4    
def secondaryVprimary(samplePath): 
    '''Calculates how big the secondary data set is compared to the primary data set (secondary read count divided by primary read count)'''  
    primPatternCount = 0
    secPatternCount = 0

    # Open read count file
    readCountFile = open(samplePath + 'tables/READ_COUNTS.txt', 'r')
    
    for line in readCountFile:
        splitLine = line.split('\t')
        record = splitLine[0]
        dataGroup = record.split('-')[0]
        
        # Sum primary read count
        if dataGroup == '4':
            primPatternCount += int(splitLine[2])
            
        # Sum secondary read count
        elif dataGroup == '5':
            secPatternCount += int(splitLine[2])           
    readCountFile.close()
    
    # How big is the secondary data compared to primary data
    size = secPatternCount / primPatternCount

    # Write to sample report
    writeReadCount(samplePath, primPatternCount, secPatternCount, size)

    # If the size of secondary data is bigger than threshold; return True
    if size > maxSecondarySize:
        return True
    else:
        return False


#FLAG5
def secondaryAssembly(samplePath):
    '''Checks for the presence of "secondary-assembly" directory'''
    sampleReport = openFile(samplePath + sample + '_contaminationReport.txt', 'a')
    sampleReport.write('\n')
    sampleReport.write('#' * 37 + '\n')
    sampleReport.write('#          Secondary Data           #\n')
    sampleReport.write('#' * 37 + '\n\n')

    # Did IRMA make a secondary assembly?
    if os.path.exists(samplePath + 'secondary_assembly') and os.path.isdir(samplePath + 'secondary_assembly'):
        sampleReport.write('(5) A secondary assembly was made! This indicates high amount of secondary data!\n\n') 
        sampleReport.close()
        summarizeSecondary(samplePath, primaryGenus)
        return True
    else:
        sampleReport.write('(5) No secondary assembly was made\n\n')    
        sampleReport.close()
        summarizeSecondary(samplePath, primaryGenus)
        return False



#----For writing output----#
def writeCount(samplePath, flaggedForCount, countDict, segments):
    sampleReport = openFile(samplePath + sample + '_contaminationReport.txt', 'a')
    sampleReport.write('\n')
    sampleReport.write('#' * 37 + '\n')
    sampleReport.write('#         Minority Variants         #\n')
    sampleReport.write('#' * 37 + '\n\n')
    
    # Write header
    sampleReport.write("{:<16}".format(''))
    for s in segments:
        if s in flaggedForCount:
            sampleReport.write("{:<8}".format(s + '*'))
        else:
            sampleReport.write("{:<8}".format(s))
    sampleReport.write('\n')
    
    #Write count thresholds
    sampleReport.write("{:<3} {:<12}".format('', 'MaxCount'))
    for s in segments:
        maxCount = math.floor(segmentLengths[s] * maxVariantPercentage)
        sampleReport.write("{:<8}".format(maxCount))
    sampleReport.write('\n')
    
    # Write counts
    sampleReport.write("{:<3} {:<12}".format('(1)', 'Count'))
    for s in segments:
        if s in countDict:
            sampleReport.write("{:<8}".format(str(countDict[s]))) 
        else:
            sampleReport.write("{:<8}".format('NA')) 
    sampleReport.write('\n')
    sampleReport.close()



def writeFreq(samplePath, highFreqCountDict, segments):
    '''Write count of high frequency variants to sample report'''
    sampleReport = openFile(samplePath + sample + '_contaminationReport.txt', 'a')

    # Write counts of variants with frequency above maxFreq
    sampleReport.write('\n')
    sampleReport.write("{:<3} {:<12}".format('(2)', 'Freq>' + str(maxFreq)))
    for s in segments:     
        if s in highFreqCountDict:
            sampleReport.write("{:<8}".format(highFreqCountDict[s]))
        else:
            sampleReport.write("{:<8}".format('0'))  
    sampleReport.write('\n')
    sampleReport.close()



def writeProp(samplePath, propDict, segments):
    '''Write proportion of high frequency variants to sample report'''
    sampleReport = openFile(samplePath + sample + '_contaminationReport.txt', 'a')
    
    sampleReport.write("{:<3} {:<12}".format('(3)', 'Proportion'))
    for s in segments:
        if s in propDict:
            sampleReport.write("{:<8}".format(propDict[s]))
        else:
            sampleReport.write("{:<8}".format('-')) 
    sampleReport.write('\n')   
    sampleReport.close()



def writeReadCount(samplePath, primPatternCount, secPatternCount, size):
    '''Write summary of read counts in sample and the size of '''   
    sampleReport = openFile(samplePath + sample + '_contaminationReport.txt', 'a')
    sampleReport.write('\n\n\n\n')
    sampleReport.write('#' * 37 + '\n')
    sampleReport.write('#        Secondary v. primary       #\n')
    sampleReport.write('#' * 37 + '\n')
    sampleReport.write('\n')

    # Size of secondary data compared to primary data
    if size >= 0.01:
        sampleReport.write('(4) Secondary data is ' + str(round(size * 100, 2)) + '% the size of primary data\n\n')
    else:
        sampleReport.write('(4) Secondary data is <1% of the size of primary data\n\n')

    # Total number of read patterns
    sampleReport.write("{:<10} {:<20}".format('', 'TotalReadPatterns') + '\n')
    sampleReport.write("{:<10} {:<20}".format('Primary', str(primPatternCount)) + '\n')
    sampleReport.write("{:<10} {:<20}".format('Secondary', str(secPatternCount)) + '\n\n\n\n')
    sampleReport.close()



def writeContaminationReport(root, flagList):
    '''Writing output to contamination report'''
    contaminationReport = open(root + '/contaminationReport.txt', 'a')

    # Does this sample have too many flags raised?
    flagCount = flagList.count(True)
    if flagCount >= maxFlags:
        flagged = 'True'
    else:
        flagged = 'False'

    # Write results
    contaminationReport.write("{:<10} {:<10} {:^3} {:^10} {:^10} {:^10} {:^10} {:^10}".format(sample, flagged, '|', flagList[0], flagList[1], flagList[2], flagList[3], flagList[4]) + '\n')
    contaminationReport.close()
    



##############################################################################
#                                   MAIN                                     #
##############################################################################
# Initialize contamination report file
initContaminationReport(root)

# Go through found directories
dirList = findDirectories(root)
for sample in dirList:
    samplePath = root + '/' + sample + '/'

    # Proceed only if this is a sample folder
    if isSampleFolder(samplePath) is True: 

        # Initialize
        flag1 = False
        flag2 = False
        flag3 = False
        flag4 = False
        flag5 = False
        printSample(sample)
        initSampleReport(samplePath)
        primaryGenus, segments = primaryData(samplePath)



        #####
        # 1 # Count minority variants
        #####
        print('--------1. Minority variant count--------')
        flag1, flaggedForCount, countDict = variantCount(samplePath, segments)
        print('\n')
        
        
        
        if flag1 is True:
            #####
            # 2 # Check frequencies in segments that were flagged as having too many minority frequencies (if any)
            #####
            print('--------2. Minority variant frequency--------')
            flag2, flaggedForFreq, freqDict, highFreqCountDict = variantFrequency(samplePath, flaggedForCount, countDict, segments)
            print('\n')
            
            
                        
            if flag2 is True:              
                #####
                # 3 # Check proportion of high frequency minority variants (if any)
                #####
                print('--------3. High frequency proportion--------')
                flag3, flaggedForProp, highFreqCountDict, propDict = highFreqProportion(samplePath, freqDict, highFreqCountDict, flaggedForFreq, segments)
                print('\n')
                



        #####
        # 4 # Calculate size of secondary data compared to primary data
        #####
        print('--------4. Secondary data size--------')
        flag4 = secondaryVprimary(samplePath)
        print('\n')



        #####
        # 5 # Look for presence of 'Secondary_assembly' folder
        #####
        print('--------5. Secondary assembly--------')
        flag5 = secondaryAssembly(samplePath)
        

        # Write results to contamination report
        writeContaminationReport(root, [flag1, flag2, flag3, flag4, flag5])
    

print('\n\nContamination report done!')

