# Contamination analysis of IRMA output

## Content

1. How to run
2. Intro to program
    1.1 Thresholds
    1.2 Flag description
3. Output
    3.1 Output file description
4. Potential issues


## 1. How to run

Navigate to the folder from where you ran IRMA; this is where all your sample output folders are. Run the analysis simply by:

```sh
python3 /path/to/contaminationAnalysis.py
```

If the program ran without errors it will say '_Contamination report done_'.


## 2. Intro to program

The general idea of this program is to look for signs of potential contamination and/or coinfection in influenza samples run through [IRMA FLU](https://wonder.cdc.gov/amd/flu/irma/). There are currently 5 criteria being checked and for each criteria a flag can either be raised or not.
- **Flag1** - Variant count
- **Flag2** - Variant frequency
- **Flag3** - High frequency variant proportion
- **Flag4** - Secondary data size
- **Flag5** - Secondary assembly

If the number of raised flags exceeds or is equal to _maxFlags_, the whole sample will be flagged in the output as potentially contaminated/coinfected. 

## 2.1 Thresholds
You can change the thresholds in the beginning of the script
```sh
maxVariantPercentage = 0.0028   #Percentage of segment length expected to have variants (0.28%)
maxFreq = 0.02                  #Frequency at which a variant is considered high frequency (2%)
maxFreqProp = 0.10              #Max proportion of high frequency variants (10%)
maxSecondarySize = 0.10         #Max size of secondary read pattern sum compared to primary read pattern sum (10%)
maxFlags = 2                    #Number of flags to be raised before the entire sample is flagged
```

## 2.2 Flag description

***

##### Flag1 - Variant count 

Counting number of minority variants in each segment. If any segment have a variant count above the threshold, flag1 will be raised (ie. set to True). This threhold is calculated using a percentage (_maxVariantPercentage_) denoting how big a part of the segment we'd expect to have variants without it being suspicious. This percentage was calculated from the HA segment; expecting maximum 5 variants divided by the length of the HA segment, which is 1778 bases:

 ```sh
 maxVariantPercentage = 5 / 1778 = 0.0028
```

The maximum number of the variants allowed is calculated for each segment by multiplying _maxVariantPercentage_ by the length of the segment and rounding down. For example, the number of allowed variants in the PB2 segment is:
 ```sh
 maxCount = 0.0028 * 2341 = 6.5548 = 6 variants allowed
```

If flag1 is not raised, the program will skip to check the criteria for flag4. If flag1 is raised, the program will pass a list of segments with too many variants, called _flaggedForCount_, and continue to check the criteria for flag2 . 

 ***
 
##### Flag2 - Variant frequency

Are there any variants in _flaggedForCount_ segments  with a frequency above _maxFreq_ (set to 2%)? If there is at least one, flag2 will be raised. Honesty time: I know this flag is not strictly necessary since only the amount of high frequency variants is important. I made this flag to help me better understand what I was doing, but in the future this flag could be removed. 

If flag2 is not raised, the program will skip to flag4. If flag2 is raised the program will pass an updated version of _flaggedForCount_ with segments also flagged for having at least one high frequency variant, called _flaggedForFreq_. It will then continue to check the criteria for flag3.

****

##### Flag3 - High frequency variant proportion

Out of all variants in a segment, is too high a proportion high frequency variants? If the proportion of high frequency variants exceeds the _maxFreqProp_ (set to 10%) flag3 will be raised. The proportion is calculated as:
 ```sh
Number of high frequency variants / Total number of variants
```
No matter if flag3 is raised or not, the program will continue on to flag4.

***

##### Flag4 - Secondary data size

Is the size of the secondary data too large compared to the primary data? First the number of read patterns is summed for both the primary and secondary data, from the READ_COUNTS.txt file (I'm using read patterns and not read counts, because IRMA does it when they calculate their observed factor. I'm sure there is a point in not including duplicate reads, but my brain hasn't understood the mathematics behind it). The size of the secondary data is calculated as:
 ```sh
Total secondary read pattern count / total primary read pattern count
```

If this number exceeds _maxSecondarySize_ (set to 10%) flag4 will be raised. For example, if the secondary data size is 5%, it means the secondary data is 5% of the size of the primary data, which is not enough to raise the flags

No matter if flag4 is raised or not, the program will continue on to flag5.

***

#####  Flag5 - Secondary assembly

Is a secondary assembly dictionary present in the IRMA output? If yes, it means there was enough secondary data to trigger an assembly (IRMA calculates an observed factor for this. You can read about it under [Residual and secondary assembly](https://wonder.cdc.gov/amd/flu/irma/) on the IRMA website), and flag5 is raised. 

This part could be expanded on in the future to also analyse the reads placed in the 'secondary folder'. This program summarises these reads in the output, but doesn't use it for flagging. 

End of program

***


## 3. Output

The program output two files:
- _contaminationReport.txt_ - Placed in your current working directory. One file summarizing results for every sample. 
- _prefix_contaminationAnalysis.txt_ - One for each sample, placed in their respective sample folder

### 3.1 Output file description

***

#### _contaminationReport.txt_

For every sample the _Flagged_ column denotes whether the sample has been flagged as potentially contaminated or coinfected, with True or False. The following 5 columns show which flags were raised, with 1 or 0 (1=flag raised, 0=flag not raised).

##### Example of output:

|Sample|Flagged|(1)Count|(2)Freq|(3)Count|(4)Reads|(5)Assembly|
|------|-------|:------:|:-----:|:------:|:------:|:---------:|
|Day1  |True   |1       |1      |1       |0       |0          |
|Day2  |False  |1       |0      |0       |0       |0          |
|Day3  |True   |1       |1      |0       |0       |0          |
NB! Ideally I don't think the last sample in this example should be flagged as potentially contaminated. There might just be one high frequency variant, but it shouldn't matter if the proportion isn't too high.

***

#### _prefix_contaminationAnalysis.txt_

Every sample will have this file containing a summary of the results found during the analysis. If a sample was flagged as being potentially contaminated or coinfected, you can go look in this file to get a better idea why. 

First the file shows the primary influenza genus and subtype, follwed by three summary parts:
1. **Minority variants** - A table showing, for each segment, the maximum number of allowed variants, the variant count, count of high frequency variants and the proportion of high frequency variants. NOTE that the frequency and proportion lines are only added if any segments have too many variants. So, if you see a file mssing these last two lines, it is not an error but just means there is no issues with the minority variants :)
Segments with at least one flag raised in this part is marked with *.
2. **Secondary v. primary** - Shows the size of the secondary data compared to the primary data, followed by a table with the summed read pattern counts. Use this to see if there is a lot of secondary data, which can be caused by contamination or coinfection.
3. **Secondary data** - A summary of the read counts for secondary subtypes and genera found in the sample.

It is shown which flag the results belong to with the number of the flag in parenthesis ((1), (2), etc...)

##### Example of output:

Primary data is influenza A (H3N2)

##### Minority Variants

|              |PB2  |PB1  |PA*  |HA_H3|NP  |NA_N2|MP  |NS  |
|--------------|:---:|:---:|:---:|:---:|:--:|:---:|:--:|:--:|
|MaxCount      |6    |6    |6    |4    |4   |3    |2   |2   |
|(1) Count     |0    |0    |8    |1    |2   |1    |0   |0   |
|              |     |     |     |     |    |     |    |    |
|(2) Freq>0.02 |0    |0    |4    |1    |0   |0    |0   |0   |
|(3) Proportion|-    |-    |0.5  |-    |-   |-    |-   |-   |

##### Secondary v. Primary
(4) Secondary data is <1% of the size of primary data

|         |TotalReadPatterns|
|---------|-----------------|
|Primary  |305058           |
|Secondary|11               |

##### Secondary Data 

(5) No secondary assembly was made

\-------Subtypes of influenza A-------
|Segment|PatternCount|
|-------|------------|
|H1     |3           |
|H8     |1           |
|N1     |6           |


\----------Secondary genera-----------
|Genus|Segment|ReadCount|
|-----|-------|---------|
|B    |HA     |1        |

***

## Potential issues
- If there is more than one genera in the primary data (for example, A_HA and B_NA) the program will exit. I don't know if this would ever happen, but I imagine it could if a sample has been heavily contaminated, so it's 50/50 influenza A and B
- Flag2 is a bit useless
- If the HA and/or NA segment have failed to be sequenced and not are present in the primary data, I think the program will just give up. 
- Not much error handling
- Adding something with phases

