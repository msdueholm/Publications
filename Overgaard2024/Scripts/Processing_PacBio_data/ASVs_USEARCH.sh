#!/bin/bash


### log
exec > denoising.log 2>&1


echo Dereplicating 
usearch11 -fastx_uniques /user_data/cko/Projects/P004/PacBio_data/Processing/mock/Orient_vsearch/filtered_default_oriented.fasta -fastaout uniques.fa -sizeout -relabel uniq

echo
echo Denoising minsize8
usearch11 -unoise3 uniques.fa -zotus zotus_minsize8.fa

echo
echo Denoising minsize7
usearch11 -unoise3 uniques.fa -minsize 7 -zotus zotus_minsize7.fa

echo
echo Denoising minsize6
usearch11 -unoise3 uniques.fa -minsize 6 -zotus zotus_minsize6.fa

echo
echo Denoising minsize5
usearch11 -unoise3 uniques.fa -minsize 5 -zotus zotus_minsize5.fa

echo
echo Denoising minsize4
usearch11 -unoise3 uniques.fa -minsize 4 -zotus zotus_minsize4.fa

echo
echo Denoising minsize3
usearch11 -unoise3 uniques.fa -minsize 3 -zotus zotus_minsize3.fa

echo
echo Denoising minsize2
usearch11 -unoise3 uniques.fa -minsize 2 -zotus zotus_minsize2.fa

module purge 

echo Denoising "DONE! 
