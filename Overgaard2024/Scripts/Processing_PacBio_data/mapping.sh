#!/bin/bash


### log
exec > mapping.log 2>&1

echo vsearch oriented data
echo 1
usearch11 -search_global /user_data/cko/Projects/P004/PacBio_data/Processing/mock/Orient_vsearch/filtered_strict_oriented.fasta -db /user_data/cko/Projects/P004/P004_J002_Mock/EUK_mock_ref_260123.fa -id 0 -maxrejects 0 -maxaccepts 0 -strand plus -top_hit_only -blast6out filtered_oriented_strict_vsearch_tophitonly_search_mockdb.b6 -threads 60

echo
echo 2
usearch11 -usearch_global /user_data/cko/Projects/P004/PacBio_data/Processing/mock/Orient_vsearch/filtered_default_oriented.fasta -db /user_data/cko/Projects/P004/P004_J002_Mock/EUK_mock_ref_260123.fa -id 0 -maxrejects 0 -maxaccepts 0 -strand plus -top_hit_only -blast6out filtered_oriented_default_vsearch_tophitonly_mockdb.b6 -threads 60

echo
echo 3
usearch11 -search_global /user_data/cko/Projects/P004/PacBio_data/Processing/mock/Orient_vsearch/filtered_default_oriented.fasta -db /user_data/cko/Projects/P004/P004_J002_Mock/EUK_mock_ref_260123.fa -id 0 -maxrejects 0 -maxaccepts 0 -strand plus -top_hit_only -blast6out filtered_oriented_default_vsearch_tophitonly_search_mockdb.b6 -threads 60


echo Usearch oriented data
echo
echo 4
usearch11 -usearch_global /user_data/cko/Projects/P004/PacBio_data/Processing/mock/Orient_usearch/filtered_strict_oriented.fasta -db /user_data/cko/Projects/P004/P004_J002_Mock/EUK_mock_ref_260123.fa -id 0 -maxrejects 0 -maxaccepts 0 -strand plus -top_hit_only -blast6out filtered_oriented_strict_usearch_tophitonly_mockdb.b6 -threads 60

echo
echo 5
usearch11 -search_global /user_data/cko/Projects/P004/PacBio_data/Processing/mock/Orient_usearch/filtered_strict_oriented.fasta -db /user_data/cko/Projects/P004/P004_J002_Mock/EUK_mock_ref_260123.fa -id 0 -maxrejects 0 -maxaccepts 0 -strand plus -top_hit_only -blast6out filtered_oriented_strict_usearch_tophitonly_search_mockdb.b6 -threads 60

echo
echo 6
usearch11 -usearch_global /user_data/cko/Projects/P004/PacBio_data/Processing/mock/Orient_usearch/filtered_default_oriented.fasta -db /user_data/cko/Projects/P004/P004_J002_Mock/EUK_mock_ref_260123.fa -id 0 -maxrejects 0 -maxaccepts 0 -strand plus -top_hit_only -blast6out filtered_oriented_default_usearch_tophitonly_mockdb.b6 -threads 60

echo
echo 7
usearch11 -usearch_global /user_data/cko/Projects/P004/PacBio_data/Processing/mock/Orient_usearch/filtered_default_oriented.fasta -db /user_data/cko/Projects/P004/P004_J002_Mock/EUK_mock_ref_260123.fa -id 0 -maxrejects 0 -maxaccepts 0 -strand plus -top_hit_only -blast6out filtered_oriented_default_usearch_tophitonly_search_mockdb.b6 -threads 60

echo DADA2 trimmed data 
echo
echo 
usearch11 -usearch_global /user_data/cko/Projects/P004/PacBio_data/Processing/mock/dada2/filtered_default/filtered/demultiplex.bc1008--bc1008.hifi_reads.fasta -db /user_data/cko/Projects/P004/P004_J002_Mock/EUK_mock_ref_260123.fa -id 0 -maxaccepts 0 -maxrejects 0 -strand plus -top_hit_only -blast6out dada2_default_tophitonly_mockdb.b6 -threads 60

echo
echo 9
usearch11 -search_global /user_data/cko/Projects/P004/PacBio_data/Processing/mock/dada2/filtered_default/filtered/demultiplex.bc1008--bc1008.hifi_reads.fasta -db /user_data/cko/Projects/P004/P004_J002_Mock/EUK_mock_ref_260123.fa -id 0 -maxaccepts 0 -maxrejects 0 -strand plus -top_hit_only -blast6out dada2_default_tophitonly_search_mockdb.b6 -threads 60

echo
echo 10
usearch11 -usearch_global /user_data/cko/Projects/P004/PacBio_data/Processing/mock/dada2/filtered_strict/demultiplex.bc1008--bc1008.hifi_reads.filtered.fasta -db /user_data/cko/Projects/P004/P004_J002_Mock/EUK_mock_ref_260123.fa -id 0 -maxaccepts 0 -maxrejects 0 -strand plus -top_hit_only -blast6out dada2_strict_tophitonly_mockdb.b6 -threads 60

echo
echo 11
usearch11 -search_global /user_data/cko/Projects/P004/PacBio_data/Processing/mock/dada2/filtered_strict/demultiplex.bc1008--bc1008.hifi_reads.filtered.fasta -db /user_data/cko/Projects/P004/P004_J002_Mock/EUK_mock_ref_260123.fa -id 0 -maxaccepts 0 -maxrejects 0 -strand plus -top_hit_only -blast6out dada2_strict_tophitonly_search_mockdb.b6 -threads 60

echo Raw data 
echo
echo 12
usearch11 -usearch_global /user_data/cko/Projects/P004/PacBio_data/Processing/mock/demultiplex.bc1008--bc1008.hifi_reads.fastq -db /user_data/cko/Projects/P004/P004_J002_Mock/EUK_mock_ref_260123.fa -id 0 -maxrejects 0 -maxaccepts 0 -strand both -top_hit_only -blast6out raw_data_tophitonly_mockdb.b6 -threads 60

echo
echo 13
usearch11 -search_global /user_data/cko/Projects/P004/PacBio_data/Processing/mock/demultiplex.bc1008--bc1008.hifi_reads.fastq -db /user_data/cko/Projects/P004/P004_J002_Mock/EUK_mock_ref_260123.fa -id 0 -maxrejects 0 -maxaccepts 0 -strand both -top_hit_only -blast6out raw_data_tophitonly_search_mockdb.b6 -threads 60


echo script done