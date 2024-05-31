#!/bin/bash

#Usage: orientation.sh <path to orientation.sh>

DIR=$1
PR2="/user_data/cko/Reference_db/PR2/v5_beta_mj/pr2_version_5_mj.fasta"

## log
exec > orientation.log 2>&1


module load VSEARCH/2.22.1-GCC-11.3.0

##Orient filtered reads using PR2 - default filtered reads 
echo Orient filtered reads using PR2
echo
vsearch --orient /user_data/cko/Projects/P004/PacBio_data/Processing/mock/dada2/filtered_default/filtered/demultiplex.bc1008--bc1008.hifi_reads.fasta -db $PR2 -fastaout filtered_default_oriented1.fasta -notmatched filtered_default_notmatched1.fasta -tabbedout filtered_default_orient1.txt

echo orientation of undetected reads 
echo
vsearch --orient filtered_default_notmatched1.fasta -db filtered_default_oriented1.fasta -fastaout filtered_default_oriented2.fasta -notmatched filtered_default_notmatched2.fasta -tabbedout filtered_default_orient2.txt

echo second orientation of undetected reads 
echo
vsearch --orient filtered_default_notmatched2.fasta -db filtered_default_oriented2.fasta -fastaout filtered_default_oriented3.fasta -notmatched filtered_default_notmatched3.fasta -tabbedout filtered_default_orient3.txt
echo
echo
echo Do manual screening before concatenation of oriented files to see of the last orientation of undetected reads is necessary

##Orient filtered reads using PR2 - strickt filtered reads 
echo Orient filtered reads using PR2
echo
vsearch --orient /user_data/cko/Projects/P004/PacBio_data/Processing/mock/dada2/filtered_strict/demultiplex.bc1008--bc1008.hifi_reads.filtered.fasta -db $PR2 -fastaout filtered_strict_oriented1.fasta -notmatched filtered_strict_notmatched1.fasta -tabbedout filtered_strict_orient1.txt

echo orientation of undetected reads 
echo
vsearch --orient filtered_strict_notmatched1.fasta -db filtered_strict_oriented1.fasta -fastaout filtered_strict_oriented2.fasta -notmatched filtered_strict_notmatched2.fasta -tabbedout filtered_strict_orient2.txt

echo second orientation of undetected reads 
echo
vsearch --orient filtered_strict_notmatched2.fasta -db filtered_strict_oriented2.fasta -fastaout filtered_strict_oriented3.fasta -notmatched filtered_strict_notmatched3.fasta -tabbedout filtered_strict_orient3.txt
echo
echo
echo Do manual screening before concatenation of oriented files to see of the last orientation of undetected reads is necessary


#cat them together and but them into folders 
#echo Concatenating oriented files
#echo
#cat consensus_oriented1.fasta consensus_oriented2.fasta > consensus_oriented.fasta
#echo
#echo Done concatenating oriented files 

module purge