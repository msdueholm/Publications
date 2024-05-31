#!/bin/bash

#Usage: orientation.sh 

THREADS=22
PR2="/user_data/cko/Reference_db/PR2/v5_beta_mj/usearch_udb/pr2_version_5_mj.udb"


exec > orientation.log 2>&1

##Orient filtered reads using PR2
echo Orientation of concensus sequences using PR2
echo
usearch11 -orient /user_data/cko/Projects/P004/P004_J002_Mock/Processing_canu/Medaka/consensus.fasta -db $PR2 -threads $THREADS -fastaout consensus_oriented1.fasta -notmatched consensus_notmatched1.fasta -tabbedout consensus_orient1.txt

echo Orientation of undetected reads 
echo
usearch11 -orient consensus_notmatched1.fasta -db consensus_oriented1.fasta -threads $THREADS -fastaout consensus_oriented2.fasta -notmatched consensus_notmatched2.fasta -tabbedout consensus_orient2.txt

echo Second orientation of undetected reads 
echo
usearch11 -orient consensus_notmatched2.fasta -db consensus_oriented2.fasta -threads $THREADS -fastaout consensus_oriented3.fasta -notmatched consensus_notmatched3.fasta -tabbedout consensus_orient3.txt

echo Do manual screening before concatenation of oriented files to see of the last orientation of undetected reads is necessary

