#!/bin/bash

#Usage: subset.sh <path to Processing_canu folder>

DIR=$1


#Rename the headers of each concatenated file as the headers are to long for medaka to handle

#mkdir $DIR/subsets_renamed/

cd $DIR/subsets_renamed/

exec > subset.log 2>&1

module load seqtk/1.3-GCC-10.2.0

for file in /user_data/cko/Projects/P004/P004_J002_Mock/Processing/fastq_pass/Concatenated/*.fastq; do base=$(echo $file | awk -F '.' '{print $1}'); echo "subsetting $base..."; seqtk sample -s100 $file 4000 > "$base"_subset.fastq; done

module purge

mv /user_data/cko/Projects/P004/P004_J002_Mock/Processing/fastq_pass/Concatenated/*_subset.fastq /user_data/cko/Projects/P004/P004_J002_Mock/Processing_canu/subsets_renamed/

exec > renamed.log 2>&1

echo Start renaming 
echo
for file in *_subset.fastq; do base=$(echo $file | awk -F '_' '{print $1}'); echo "renaming $base..."; sed '/^@.*runid=/ {s/^@.*read=/@read=/g; s/_.*barcode=/_barcode=/g}' $file > "$base"_renamed.fastq; done

echo Renaming done!
echo



