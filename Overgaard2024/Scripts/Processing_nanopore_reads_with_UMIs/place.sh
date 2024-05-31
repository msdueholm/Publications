#!/bin/bash

TAX=/user_data/cko/Reference_db/PR2/v5_beta_mj/pr2_version_5_mj.tax

### log
exec > place.log 2>&1

## Strategy 2 taxonomy - prune away the OTUs and phylogenetically place them back on the reference-only tree. Get taxonomy and confidence scores for each taxonomic rank.

## split alignment into OTUs and references - this is important for running EPA-ng later!
module load SeqKit/2.0.0

for i in 18S_Askov_UMI.ref.mafft.trimal.fasta; do base=$(echo $i | awk -F '.' '{print $1}'); cat $i | seqkit grep -r -p ASV > 18S_Askov_UMI.otus.fasta; done
for i in 18S_Askov_UMI.ref.mafft.trimal.fasta; do base=$(echo $i | awk -F '.' '{print $1}'); cat $i | seqkit grep -r -p ASV -v > 18S_Askov_UMI.ref.fasta; cat "$base".ref.fasta | grep ">" | tr -d '>' > 18S_Askov_UMI.ref.list; done 

module purge

## prune tree - use a python script for this (uses the ete toolkit)
module load ETE/3.1.2-foss-2020b #Axomomma
#module load ETE/3.1.2-foss-2020b-Python-3.8.6 #i768

for i in 18S_Askov_UMI.RAxML_fastTreeSH_Support.18S_Askov_UMI.fixed; do base=$(echo $i | awk -F '.' '{print $1}'); python /user_data/cko/Projects/P003/Mahwash_scripts/Original_scripts/prune.py $i 18S_Askov_UMI.ref.list 18S_Askov_UMI.pruned.tre; done

module purge

## get model details for use in EPA and optimise branch lengths
module load RAxML-NG/1.0.3-GCC-10.2.0

for i in 18S_Askov_UMI.pruned.tre; do base=$(echo $i | awk -F '.' '{print $1}'); raxml-ng --evaluate --tree $i --msa 18S_Askov_UMI.ref.fasta --model GTR+G --threads 4 --prefix 18S_Askov_UMI.pruned; done

module purge

## run EPA-ng
module load EPA-NG/0.3.8-GCC-10.2.0

for i in *pruned.raxml.bestTree; do base=$(echo $i | awk -F '.' '{print $1}'); epa-ng --tree $i --ref-msa 18S_Askov_UMI.ref.fasta --query 18S_Askov_UMI.otus.fasta --model 18S_Askov_UMI.pruned.raxml.bestModel --threads 10 --no-heur --no-pre-mask --verbose; mv epa_info.log 18S_Askov_UMI.epa_info.log; mv epa_result.jplace 18S_Askov_UMI.epa_result.jplace; done

module purge

## set up taxonomy file

### fasta headers
for i in 18S_Askov_UMI.ref.fasta; do base=$(echo $i | awk -F '.' '{print $1}'); cat $i | grep ">" | tr -d '>' | sed -E 's/(.*@)(.*)/\2\t\1\2/' | sort -k1,1 > 18S_Askov_UMI.ref.headers; done

### list of accessions
for i in 18S_Askov_UMI.ref.fasta; do base=$(echo $i | awk -F '.' '{print $1}'); cat $i | grep ">" | sed -E 's/(.*@)(.*)/\2/' > 18S_Askov_UMI.accessions; done

### grep for taxonomy

for i in 18S_Askov_UMI.accessions; do base=$(echo $i | awk -F '.' '{print $1}'); grep -w -f $i ${TAX} | sed -E 's/(.*)\.[0-9]+\..*(\tEukaryota.*)/\1\2/' | sed -E 's/;$//' | sort -u | sort -k1,1 > 18S_Askov_UMI.ref.tax; done

## join the files based on the common accession number
for i in 18S_Askov_UMI.ref.headers; do base=$(echo $i | awk -F '.' '{print $1}'); join $i 18S_Askov_UMI.ref.tax | tr ' ' '\t' | cut -f 2,3 > 18S_Askov_UMI.taxonomy.txt; done


## assign taxonomy with gappa
module load Gappa/0.7.1-GCC-10.2.0

for i in 18S_Askov_UMI.ref.fasta; do base=$(echo $i | awk -F '.' '{print $1}'); grep "Opisthokonta" $i | tr -d '>' > 18S_Askov_UMI.outgroup; done

for i in 18S_Askov_UMI.epa_result.jplace; do base=$(echo $i | awk -F '.' '{print $1}'); gappa examine assign --threads 10 --jplace-path $i --sativa --taxon-file 18S_Askov_UMI.taxonomy.txt --root-outgroup 18S_Askov_UMI.outgroup --per-query-results --sativa --file-prefix 18S_Askov_UMI.assign_; done

module purge 

## Strategy 1 taxonomy - get taxonomy of OTUs based on their position on the tree. See "Long-read metabarcoding of the eukaryotic rDNA operon to phylogenetically and taxonomically resolve environmental diversity" for details.
## Uses the genesis app "partial-tree-taxassign"
module load partial_tree_tax_assign/20211015-GCC-10.2.0

for i in *taxonomy.txt; do base=$(echo $i | awk -F '.' '{print $1}'); partial-tree-taxassign 18S_Askov_UMI.RAxML_fastTreeSH_Support.18S_Askov_UMI.fixed $i 18S_Askov_UMI.outgroup > 18S_Askov_UMI.assignment.tsv; done

module purge

## Combine taxonomy from both strategies.
module load Array-Utils/0.5-foss-2020b

for i in *assignment.tsv; do base=$(echo $i | awk -F '.' '{print $1}'); perl /user_data/cko/Projects/P003/Mahwash_scripts/Original_scripts/taxonomy_merge.pl $i 18S_Askov_UMI.assign_sativa.tsv 0.5 18S_Askov_UMI.taxonomy_merged.tsv; done

module purge

### I recommend scanning through the table manually. If the two taxonomy strategies give highly conflicting results, you can discard the sequence after checking it manually, or re-label it manually if you are confident that the sequence is not a chimera.
