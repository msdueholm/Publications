#!/bin/bash

#Usage: filtering.sh <path to consensus file>

DIR=$1
THREADS=65
PR2="/srv/PHN/users/cko/Reference_db/PR2/v5_beta_mj/usearch_udb/pr2_version_5_mj.udb"

## log
exec > orientation.log 2>&1

##Orient filtered reads using PR2
echo Orient filtered reads using PR2
echo
usearch11 -orient /srv/PHN/users/cko/Projects/P004/Undiluted_CKO_P004_J017/nanopore_pipeline/consensus_raconx3_medakax2.fa -db $PR2 -threads $THREADS -fastaout consensus_oriented1.fasta -notmatched consensus_notmatched1.fasta -tabbedout consensus_orient1.txt

echo orientation of undetected reads 
echo
usearch11 -orient consensus_notmatched1.fasta -db consensus_oriented1.fasta -threads $THREADS -fastaout consensus_oriented2.fasta -notmatched consensus_notmatched2.fasta -tabbedout consensus_orient2.txt

echo second orientation of undetected reads 
echo
usearch11 -orient consensus_notmatched2.fasta -db consensus_oriented2.fasta -threads $THREADS -fastaout consensus_oriented3.fasta -notmatched consensus_notmatched3.fasta -tabbedout consensus_orient3.txt
echo
echo
echo Do manual screening before concatenation of oriented files to see of the last orientation of undetected reads is necessary

#cat them together and but them into folders 
echo Concatenating oriented files
echo
cat consensus_oriented1.fasta consensus_oriented2.fasta consensus_oriented3.fasta > consensus_oriented.fasta
echo
echo Done concatenating oriented files 
