#!/bin/bash

# usage: taxonomy_round1.sh 

DB="/user_data/cko/Reference_db/PR2/v5_beta_mj/pr2_mj.db" 
FASTA="/user_data/cko/Reference_db/PR2/v5_beta_mj/pr2_version_5_mj.fasta"
EUKREF="/user_data/cko/Reference_db/PR2/v5_beta_mj/pr2.main_groups.fasta" # subset of PR2 containing the main groups (124 seqs)
THREADS=40


### log
exec > taxonomy_round1.log 2>&1

## First fix the headers of the fasta file. Replace '/' with '_' for ease downstream
#cat 18S_Askov_umi.fasta | tr '/' '_' > 18S_Askov_umi.fixed.fasta

## Blast against PR2 database
module load BLAST+/2.13.0-gompi-2022a

blastn -query 18S_uniques.fasta -db "$DB" -num_threads $THREADS -out 18S_Askov_umi_vs_pr2.blastn -outfmt '6 std salltitles'

#module purge
## Get top 50 hits for each sequence

mkdir blasthits
for header in $(cut -f1 18S_Askov_umi_vs_pr2.blastn | sort -u); do grep -w -m 50 "$header" 18S_Askov_umi_vs_pr2.blastn | cut -f13 | sort -u > blasthits/"$header".top50hits.list; done

# Get the sequences for each of these hits
module load SeqKit/2.0.0

cd blasthits/
for i in *.top50hits.list; do cat $FASTA | seqkit grep -f "$i" > "$i".fasta; done
for i in *.fasta; do cat $i | seqkit rmdup -n -o "$i".clean; done

# Add the original query sequence to these sequences
for i in *.clean; do header=$(echo $i | awk -F '.' '{print $1}'); cat ../18S_uniques.fasta | seqkit grep -p "$header" >> "$i"; done

module purge

## Calculate ML distance and extract top 2 seqs
### align sequences
module load MAFFT/7.490-GCC-10.2.0-with-extensions

for i in *.clean; do base=$(echo $i | awk -F '.' '{print $1}'); mafft --thread $THREADS --adjustdirection $i > "$base".mafft.fasta; done

module purge

### calculate ML distance
module load RAxML/8.2.12-foss-2020b-pthreads-sse3

for i in *.mafft.fasta; do base=$(echo $i | awk -F '.' '{print $1}'); raxmlHPC-PTHREADS-SSE3 -T $THREADS -f x -p 12345 -s $i -m GTRGAMMA -n "$base".out; done
rm RAxML_info.* RAxML_parsimonyTree.*

#module purge 

# Extract top 2 hits 
for i in RAxML_distances.*; do header=$(echo $i | awk -F '.' '{print $2}'); cat $i | sort -Vk2 | grep -m 2 "$header" > "$header".top2.tsv; done

# Move to parent directory
cat *top2.tsv | cut -f1 -d ' ' >> ../top2hits.list
cd ..

# How many references?
echo Number of references: $(cat top2hits.list | sort | uniq | wc -l)
cat top2hits.list | sort | uniq > top2hits.uniq.list

# Extract these ref seqs
module load SeqKit/2.0.0

cat $FASTA | seqkit grep -f top2hits.uniq.list > top2hits.fasta

# Cat everything together
cat top2hits.fasta 18S_uniques.fasta $EUKREF > 18S_Askov_umi.ref.fasta

# Remove duplicates
cat 18S_Askov_umi.ref.fasta | seqkit rmdup -n -o 18S_Askov_umi.ref.clean.fasta

module purge

## Align
module load MAFFT/7.490-GCC-10.2.0-with-extensions

mafft --retree 2 --maxiterate 1000 --thread $THREADS --reorder --adjustdirection 18S_Askov_umi.ref.clean.fasta > 18S_Askov_umi.ref.mafft.fasta

module purge

## Check the alignment manually in Aliview

## Trim
module load TrimAl/1.4.1-foss-2020b

trimal -in 18S_Askov_umi.ref.mafft.fasta -out 18S_Askov_umi.ref.mafft.trimal.fasta -gt 0.01 -st 0.001

module purge

## Run tree with RAxML (v8 so can use GTRCAT model). An example command is shown for one tree/sample. Calculate SH-like support values.
module load RAxML/8.2.12-foss-2020b-pthreads-sse3

raxmlHPC-PTHREADS-SSE3 -T 3 -m GTRCAT -p 73927 -N 20 -s 18S_Askov_umi.ref.mafft.trimal.fasta -n 18S_Askov_umi.tre #- makes the trees 
raxmlHPC-PTHREADS-SSE3 -f J -p 12347 -m GTRGAMMA -s 18S_Askov_umi.ref.mafft.trimal.fasta -t RAxML_bestTree.18S_Askov_umi.tre -n 18S_Askov_umi # Compute SH­-like support values on the best tree
cat RAxML_fastTreeSH_Support.18S_Askov_umi | sed -E 's/\)\:([^\[]+)\[([0-9]+)\]/\)\2\:\1/g'> RAxML_fastTreeSH_Support.18S_Askov_umi.fixed # fixed the [] which contains the SH-like values. Figtree will not open files containing []

## Open  RAxML_fastTreeSH_Support.Cologne.fixed in figtree and type SH when it askes for type of value 

## Examine the tree manually (figtree) and mark nucleomorph sequences (sister to red algea) (green, HSV:120°, 100%, 100%), mislabelled reference sequences (blue, HSV:240°, 100%, 100%), and any OTU sequences that look like artefacts (ridiculously long branch for example) (magenta, HSV:300°, 100%, 100%).
#Mahwash is rooting the tree at Opistkonta 
# fasta Reorder.py - this file reorder the sequences in according to a tree so that more similar sequences are closer. This way is can be easier to see 

module purge
echo 
echo
echo script ended
echo 
