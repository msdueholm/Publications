#!/bin/bash

#Usage: bash medaka.sh < path to processing_canu folder>

#Settings
THREADS=50
DIR=$1

mkdir Medaka/
mkdir Medaka/consensus
mkdir Medaka/mapping

## log
exec > Medaka.log 2>&1

module load medaka/1.7.2

echo Alignment
echo
for f in $DIR/canu/*.contigs.fasta; do base=$(basename $f | sed 's/.contigs.fasta//g'); echo "$base - Making medaka alignment"; mini_align -i $DIR/subsets_renamed/"$base"_renamed.fasta -r $f -p Medaka/mapping/"$base"_medaka_align -threads $THREADS; done
echo
echo
echo Medaka concensus
echo
for f in Medaka/mapping/*_medaka_align.bam; do base=$(basename $f | sed 's/_medaka_align.bam//g'); echo "$base - Making  medaka polish"; medaka consensus $f Medaka/consensus/"$base"_consensus.hdf --model r941_min_sup_g507 --threads $THREADS; done
echo
echo
echo Medaka stich
echo
for f in Medaka/consensus/*_consensus.hdf; do base=$(echo $file | awk -F '_' '{print $1}'); echo "$base - Making medaka stich"; medaka stitch $f $DIR/canu/*.contigs.fasta Medaka/"$base"_medaka_polished.fasta --threads 2; done
echo
echo
echo Medaka polish done
echo

mv Medaka.log $DIR/Medaka

module purge 
