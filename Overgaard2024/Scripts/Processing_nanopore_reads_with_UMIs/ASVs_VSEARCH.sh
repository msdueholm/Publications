#!/bin/bash


# log
exec > denoising2.log 2>&1

module load VSEARCH/2.22.1-GCC-11.3.0


echo Dereplicate
vsearch --derep_fulllength /user_data/cko/Projects/P004/Undiluted_CKO_P004_J017/Orientation/vsearch/undil_oriented.fasta --output uniques.fa --sizeout --relabel uniq  

echo
echo Denosing minsize8 
vsearch -cluster_unoise uniques.fa --centroids zotus_minsize8.fa 

echo
echo uchime3
vsearch --uchime3_denovo zotus_minsize8.fa --nonchimeras zotus_minsize8_nonchimeras_uchime3.fa 

echo
echo uchime2
vsearch --uchime2_denovo zotus_minsize8.fa --nonchimeras zotus_minsize8_nonchimeras_uchime2.fa 

echo
echo uchime
vsearch --uchime_denovo zotus_minsize8.fa --nonchimeras zotus_minsize8_nonchimeras_uchime.fa 

echo
echo Denoising minsize7
vsearch -cluster_unoise uniques.fa --minsize 7 --centroids zotus_minsize7.fa 

echo
echo uchime3
vsearch --uchime3_denovo zotus_minsize7.fa --nonchimeras zotus_minsize7_nonchimeras_uchime3.fa 

echo
echo uchie2
vsearch --uchime2_denovo zotus_minsize7.fa --nonchimeras zotus_minsize7_nonchimeras_uchime2.fa 

echo
echo uchime
vsearch --uchime_denovo zotus_minsize7.fa --nonchimeras zotus_minsize7_nonchimeras_uchime.fa 

echo
echo Denoising minsize6
vsearch -cluster_unoise uniques.fa --minsize 6 --centroids zotus_minsize6.fa 

echo
echo uchime3
vsearch --uchime3_denovo zotus_minsize6.fa --nonchimeras zotus_minsize6_nonchimeras_uchime3.fa 

echo
echo uchime2
vsearch --uchime2_denovo zotus_minsize6.fa --nonchimeras zotus_minsize6_nonchimeras_uchime2.fa 

echo
echo uchime
vsearch --uchime_denovo zotus_minsize6.fa --nonchimeras zotus_minsize6_nonchimeras_uchime.fa 

echo
echo Denoising minsize5
vsearch -cluster_unoise uniques.fa --minsize 5 --centroids zotus_minsize5.fa 

echo
echo uchime3
vsearch --uchime3_denovo zotus_minsize5.fa --nonchimeras zotus_minsize5_nonchimeras_uchime3.fa 

echo
echo uchime2
vsearch --uchime2_denovo zotus_minsize5.fa --nonchimeras zotus_minsize5_nonchimeras_uchime2.fa 

echo
echo uchime
vsearch --uchime_denovo zotus_minsize5.fa --nonchimeras zotus_minsize5_nonchimeras_uchime.fa 

echo
echo Denoising minsize4
vsearch -cluster_unoise uniques.fa --minsize 4 --centroids zotus_minsize4.fa 

echo
echo uchime3
vsearch --uchime3_denovo zotus_minsize4.fa --nonchimeras zotus_minsize4_nonchimeras_uchime3.fa 

echo
echo uchime2
vsearch --uchime2_denovo zotus_minsize4.fa --nonchimeras zotus_minsize4_nonchimeras_uchime2.fa 

echo
echo uchime
vsearch --uchime_denovo zotus_minsize4.fa --nonchimeras zotus_minsize4_nonchimeras_uchime.fa 

echo 
echo Denosing minsize3
vsearch -cluster_unoise uniques.fa --minsize 3 --centroids zotus_minsize3.fa 
 
echo
echo uchime3
vsearch --uchime3_denovo zotus_minsize3.fa --nonchimeras zotus_minsize3_nonchimeras_uchime3.fa 

echo
echo uchime2
vsearch --uchime2_denovo zotus_minsize3.fa --nonchimeras zotus_minsize3_nonchimeras_uchime2.fa 

echo
echo uchime
vsearch --uchime_denovo zotus_minsize3.fa --nonchimeras zotus_minsize3_nonchimeras_uchime.fa 

echo
echo Denosing minsize2
vsearch -cluster_unoise uniques.fa --minsize 2 --centroids zotus_minsize2.fa 

echo
echo uchime3
vsearch --uchime3_denovo zotus_minsize2.fa --nonchimeras zotus_minsize2_nonchimeras_uchime3.fa 

echo
echo uchime2
vsearch --uchime2_denovo zotus_minsize2.fa --nonchimeras zotus_minsize2_nonchimeras_uchime2.fa 

echo 
echo uchime
vsearch --uchime_denovo zotus_minsize2.fa --nonchimeras zotus_minsize2_nonchimeras_uchime.fa 