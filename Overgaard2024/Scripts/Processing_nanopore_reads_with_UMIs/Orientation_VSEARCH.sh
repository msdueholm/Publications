#!/bin/bash

#Usage: orientation.sh <path to orientation.sh>

DIR=$1
PR2="/user_data/cko/Reference_db/PR2/v5_beta_mj/pr2_version_5_mj.fasta"

## log
exec > orientation.log 2>&1


module load VSEARCH/2.22.1-GCC-11.3.0

##Orient filtered reads using PR2 - undil filtered reads 
echo Orient filtered reads using PR2
echo
vsearch --orient /user_data/cko/Projects/P004/Undiluted_CKO_P004_J017/nanopore_pipeline/consensus_raconx3_medakax2.fa -db $PR2 -fastaout undil_oriented1.fasta -notmatched undil_notmatched1.fasta -tabbedout undil_orient1.txt

echo orientation of undetected reads 
echo
vsearch --orient undil_notmatched1.fasta -db undil_oriented1.fasta -fastaout undil_oriented2.fasta -notmatched undil_notmatched2.fasta -tabbedout undil_orient2.txt

echo second orientation of undetected reads 
echo
vsearch --orient undil_notmatched2.fasta -db undil_oriented2.fasta -fastaout undil_oriented3.fasta -notmatched undil_notmatched3.fasta -tabbedout undil_orient3.txt
echo
echo third orientation of undetected reads 
echo
vsearch --orient undil_notmatched3.fasta -db undil_oriented3.fasta -fastaout undil_oriented4.fasta -notmatched undil_notmatched4.fasta -tabbedout undil_orient4.txt
echo
echo
echoecho fourth orientation of undetected reads 
echo
vsearch --orient undil_notmatched4.fasta -db undil_oriented4.fasta -fastaout undil_oriented5.fasta -notmatched undil_notmatched5.fasta -tabbedout undil_orient5.txt
echo
echo Do manual screening before concatenation of oriented files to see of the last orientation of undetected reads is necessary


module purge