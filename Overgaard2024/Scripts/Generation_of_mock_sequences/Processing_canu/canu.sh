#!/bin/bash

#Usage: canu.sh <path to Processing_canu folder>

DIR=$1



#Make assembly using canu

mkdir $DIR/canu/

cd $DIR/canu/

exec > canu.log 2>&1

cd $DIR/subsets_renamed/

module load canu/2.0-foss-2018a

echo Start assembly 
echo
for file in *_renamed.fasta; do base=$(echo $file | awk -F '_' '{print $1}'); echo "Making assembly for $base..."; canu -p "$base" -d "$base"-np genomeSize=5000 -nanopore $file; done

echo Assemblies done!
echo

module purge

mv *-np/ /$DIR/canu/
