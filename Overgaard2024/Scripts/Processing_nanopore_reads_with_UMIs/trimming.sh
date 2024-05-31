#!/bin/bash

## log
exec > trimming.log 2>&1

echo Trim sequence specific primer sequences using Cutadapt
echo
module load cutadapt/3.4-GCCcore-10.2.0-Python-3.8.6

echo cutadapt loaded
echo
echo Trimming seqs 
cutadapt -b file:primers.fa consensus_oriented.fasta -j 65 > consensus_trimmed.fasta
echo
