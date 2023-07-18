#!/usr/bin/env Rscript
#this script performs a fisher test for all positions and seasons
#input is the prefix used in the AD plot analysis
#Call it from the Sweep Dynamics output directory (fx. HA-15-22_run)
args <- commandArgs(trailingOnly = TRUE)

# read data

#prefix <- 'HA-15-22'
prefix <- args[1]
mutations <- read.table(paste("yearOutput/significance/", prefix, ".subtreeMutationMap.per_position.txt", sep=""))
numIsolates <- read.table(paste("yearOutput/significance/", prefix, ".numIsolates.txt", sep=""))
mutnames <- read.table(paste("yearOutput/significance/", prefix, ".mutations.per_position.txt", sep=""), sep="\t")

#cutoff <- 0.5
cutoff <- args[2] 	#defines which frequencies are considered predominant
#min_diff <- 0.5		#defines how large the difference to previous season should be

if (numIsolates[1,1] == 0) {
	#remove first column when no isolates
	#alternative: maybe remove all columns when no isolates??
	mutations = mutations[,-1]
	numIsolates = numIsolates[,-1]
}

###########################
# MAKE SEASONs INTO YEARS #
###########################
minmaxseason <- read.table('yearOutput/significance/minmaxseason.txt')
firstYear <- as.numeric(minmaxseason[[1]][1])
lastYear <- as.numeric(minmaxseason[[2]][1])
firstSeason <- minmaxseason[[3]][1]

#Make a vector with all the seasons
seasons <- c()
for(i in firstYear:lastYear) {
  #First year has two seasons
  if(i == firstYear & firstSeason == 'N') {
    seasons <- append(seasons,c(i,i))
  #First year has one season
  } else if (i == firstYear & firstSeason == 'S') {
    seasons <- append(seasons,i)
  #Add two seasons for remaining years
  } else {
    seasons <- append(seasons,c(i,i))
  }
}

#Check to see if the last season should be removed
if(length(seasons) > ncol(numIsolates)) {
  seasons <- seasons[-length(seasons)]
}
if(length(seasons) != ncol(numIsolates)) {
  print('Warning: Something is wrong with the number of seasons')
}

#Add year as column names
colnames(mutations) <- seasons
colnames(numIsolates) <- seasons

#Merge columns belonging to same year
mutations <- t(rowsum(t(mutations), group = colnames(mutations), na.rm = T))
numIsolates <- t(rowsum(t(numIsolates), group = colnames(numIsolates), na.rm = T))

#Get years as vector and save in file (to be used for plotting later)
years <- colnames(mutations)
write(years, 'yearOutput/years.txt', ncol=length(years), sep="\t")
#################################




# set negative mutation counts to 0
# set mutation counts larger than number of isolates in season to number of isolates
# (both shouldn't happen with a working frequency correction)
mutations[mutations<0] <- 0
for (col in 1:ncol(mutations)) {
	mutations[mutations[,col]>numIsolates[1,col],col] <- numIsolates[1,col]
}

# initiate matrix to save pvalues
n_years <- ncol(mutations)
n_positions <- nrow(mutations)
pval <- matrix(nrow = n_positions, ncol = n_years-1)
#diffs <- matrix(nrow = n_positions, ncol = n_years-1)


for (year in 2:n_years) {
	num_thisyear = as.integer(numIsolates[year])
	num_prevyear = as.integer(numIsolates[year -1])


	for (position in 1:n_positions) {
		# get values for contigency table
		mut_thisyear <- as.integer(mutations[position, year])
		mut_prevyear <- as.integer(mutations[position, year - 1])
		
		nomut_thisyear <- num_thisyear - mut_thisyear
		nomut_prevyear <- num_prevyear - mut_prevyear
		
		# calculate fisher test and save pvalue in matrix
		data <- matrix(c(mut_thisyear, nomut_thisyear, mut_prevyear, nomut_prevyear), nrow = 2)
		pval[position, year -1] <- fisher.test(data, alternative = "greater")$p.value

		#diffs[position, season -1] <- mut_thisseason/num_thisseason - mut_prevseason/num_prevseason
	}
}

# adjust p-values for multiple testing
#pval_adjusted <- matrix(p.adjust(pval, method = "bonferroni"), nrow = n_positions)
pval_adjusted <- matrix(p.adjust(pval, method = "fdr"), nrow = n_positions)

# add names
pval_names <- cbind(t(mutnames), pval_adjusted)

# save pvalues into new file
write.table(pval_names, file=paste("yearOutput/significance/", prefix,".results.pvalues_fishertest.txt", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


#define frequency matrix instead if counts
frequency <- matrix(nrow = n_positions, ncol = n_years)
for (season in 1:n_years){
	frequency[,season] <- mutations[,season]/numIsolates[,season]
}

# save frequencies with names
freq_names <- cbind(t(mutnames), frequency)
write.table(freq_names, file=paste("yearOutput/significance/", prefix,".results.frequencies.txt", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# check if position was predominant before
prev_predom <- matrix(nrow = n_positions, ncol = n_years)
prev_predom[,1] <- 0

for (year in 2:n_years){
	if (year == 2) {
		prev_predom[,year] <- frequency[,1] > cutoff
	} else if (year > 2 && n_positions > 1) {
        	prev_predom[,year] <- rowSums(frequency[,1:year-1] > cutoff)
	} else if (year > 2 && n_positions == 1) {
		prev_predom[,year] <- sum(frequency[,1:year-1] > cutoff)
	}
}

# check which positions and years are significant and predominant 
sign_matrix <- cbind(rep(1, n_positions), pval_adjusted)
#diffs <- cbind(rep(0, n_positions), diffs)
#diffs[,2] <- 1
predom_sign <-matrix(nrow=n_positions, ncol=n_years)
predom_sign <- frequency > cutoff & sign_matrix < 0.05 & prev_predom == 0

predom_sign <- matrix(as.integer(predom_sign), nrow=n_positions)

# save results of testing
write.table(predom_sign, file=paste("yearOutput/significance/", prefix, ".significance.txt", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# overwrite numIsolates.txt and subtreeMutationMap.per_position.txt (but call it subtreeMutationMap) as we need the updated versions
write.table(numIsolates, file=paste("yearOutput/intermediateFiles/", prefix, ".numIsolates.txt", sep=""), row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(mutations, file=paste("yearOutput/intermediateFiles/", prefix, ".subtreeMutationMap.txt", sep=""), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

