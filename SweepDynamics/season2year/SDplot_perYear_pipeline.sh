#!/bin/bash
# This pipeline recalculates p-values and plots for an SD run using years instead of seasons

# Usage: bash SDplot_perYear_pipeline.sh <prefix>

SCRIPTPATH="${0%/*}"
PREFIX=$1
CUTOFF=0.5
FILTER="TRUE"
PLOTWIDTH=12
PLOTFORMAT="png"


#Make directories for new output
mkdir -p yearOutput
mkdir -p yearOutput/intermediateFiles
mkdir -p yearOutput/significance

#Move intermediate files from SD run to folder
mv $PREFIX".mutations.txt" "./yearOutput/intermediateFiles"
mv $PREFIX".subtreeMutationMap.txt" "./yearOutput/intermediateFiles"
mv $PREFIX"summary.txt" "./yearOutput/intermediateFiles"

# Move files for recalculating p-values to folder
mv minmaxseason.txt ./yearOutput/significance
cp output/$PREFIX".numIsolates.txt" "./yearOutput/significance"
cp output/$PREFIX".subtreeMutationMap.per_position.txt" "./yearOutput/significance"
cp output/$PREFIX".mutations.per_position.txt" "./yearOutput/significance"


#Find first year, last year and first season
minYear=$(cut -f1 yearOutput/significance/minmaxseason.txt)
maxYear=$(cut -f2 yearOutput/significance/minmaxseason.txt)
firstSeason=$(cut -f3 yearOutput/significance/minmaxseason.txt)

#Run modified Fisher's exact test
echo "-------- Fisher's exact test --------"
Rscript $SCRIPTPATH"/fishertest.R" $PREFIX $CUTOFF
echo " "


#Plot significant results
echo "-------- Plotting --------"
Rscript $SCRIPTPATH"/buildAD_sign.r" $minYear $maxYear $PREFIX $CUTOFF $SCRIPTPATH $firstSeason $FILTER $PLOTWIDTH $PLOTFORMAT > /dev/null
echo " "