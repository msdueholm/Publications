#!/bin/bash

# usage: clustering99.sh <fasta file containing 18S>

VAR=$1
BASE=$(basename "${VAR}") 
SAMPLE="${BASE%.*}"
THREADS=45

### log
exec > clustering99.log 2>&1

### Start by preclustering all filtered sequences!
echo
echo =====================================
echo "Preclustering all sequences (99% id)"
echo =====================================
echo
module load VSEARCH/2.18.0-GCC-10.2.0

mkdir preclusters
vsearch --cluster_fast $BASE --id 0.99 --threads "${THREADS}" --clusters preclusters/precluster.c- --uc precluster.uc

module purge 
### Separate preclusters > 2 sequences
mkdir preclusters/large_preclusters
cd preclusters
for i in precluster.c-*; do number=$(grep -c ">" $i); if (($number > 2)); then mv $i large_preclusters/; fi; done
cd large_preclusters

### construct count table for chimera detection downstream
for i in precluster.c-*; do x=$(echo $i | cut -d '.' -f2); count=$(grep -c ">" $i); echo -e ""$x"_conseq\t"$count"" >> $SAMPLE.count_table; done
mv $SAMPLE.count_table ./..

### Align large preclusters with mafft
module load MAFFT/7.470-foss-2020b-with-extensions
module load SeqKit/2.0.0

echo "Aligning large preclusters ..."
echo
for i in precluster.c-*; do if [ $(grep -c ">" $i) -lt 10000 ]; then mafft --quiet --thread "${THREADS}" $i > $i.aligned.fasta; else mafft --quiet --thread "${THREADS}" <(seqkit head -n 1000 $i) > $i.aligned.fasta; fi; done

module purge

### Generate majority rule consensus sequences
module load Mothur/1.41.0-foss-2018a-Python-2.7.14

for i in *.aligned.fasta; do mothur -q "#consensus.seqs(fasta=$i, cutoff=51)"; done
rm mothur*

### Change header name so that it is informative
for i in *.aligned.cons.fasta; do x=$(echo $i | cut -d '.' -f2); sed -i -E "s/(>)(.*)/\1${x}_\2/" $i; done

### mv consensus sequences fasta files back to the parent folder
mv *.aligned.cons.fasta ./..
cd ..
cat precluster.c-* > $SAMPLE.preclusters.fasta

### add small-preclusters to count table
grep ">m" $SAMPLE.preclusters.fasta | while read line; do header=$(echo $line | sed -E 's/>(.*)/\1/'); echo -e "$header\t1" >> $SAMPLE.count_table; done
echo -e "Representative_Sequence\ttotal" | cat - $SAMPLE.count_table > temp && mv temp $SAMPLE.count_table
mv $SAMPLE.count_table ./..
mv $SAMPLE.preclusters.fasta ./..
cd ..


### degap sequences
mothur "#degap.seqs(fasta=$SAMPLE.preclusters.fasta, processors="${THREADS}")"

echo Total number of preclusters: $(cat $SAMPLE.preclusters.ng.fasta | grep -c "^>")

module purge

echo script DONE!


