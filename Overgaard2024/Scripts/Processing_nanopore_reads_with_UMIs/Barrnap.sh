#!/bin/bash

# usage: Barrnap.sh <fasta file containing concensus sequences with full path>

VAR=$1
DIR=$(dirname "${VAR}")
THREADS=65

### log
exec > $DIR/Barrnap.log 2>&1

### change to right directory
cd $DIR

### Only keep sequences with one copy of 18S and 28S
echo
echo =====================================
echo Extract 18S and 28S genes 
echo =====================================
echo

module load Barrnap/0.9-foss-2018a

barrnap --threads "${THREADS}" --reject 0.4 --kingdom euk zotus_minsize2.fa > barrnap.gff

module purge 

cat barrnap.gff | grep "+" | grep "18S" | cut -f 1 | sort | uniq -u > barrnap.18S.list
cat barrnap.gff | grep "+" | grep "28S" | cut -f 1 | sort | uniq -u > barrnap.28S.list

comm -12 barrnap.18S.list barrnap.28S.list > barrnap.comm.list

module load Mothur/1.41.0-foss-2018a-Python-2.7.14

echo "Extracting seqs with 18S and 28S ..."
mothur -q "#get.seqs(accnos=barrnap.comm.list, fasta=zotus_minsize2.fa)"


echo
echo =====================================
echo Extract 18S and 28S
echo =====================================
echo

### Finally, extract 18S and 28S from each sequence

perl /user_data/cko/Projects/P003/Mahwash_scripts/Original_scripts/extract_18S-28S.pl barrnap.gff zotus_minsize2.pick.fa 18S.fasta 28S.fasta

echo Sanity check..
echo
echo Number of 18S sequences: $(cat 18S.fasta | grep -c ">")
echo Number of 28S sequences: $(cat 28S.fasta | grep -c ">")


echo DONE!


