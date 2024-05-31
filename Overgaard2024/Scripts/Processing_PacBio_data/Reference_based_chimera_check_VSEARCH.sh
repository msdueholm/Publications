#!/bin/bash


### log
exec > chimer_check_ref.log 2>&1

module load VSEARCH/2.22.1-GCC-11.3.0

echo
echo uchime_ref_default_unfiltered_minsize2
vsearch --uchime_ref /user_data/cko/Projects/P004/PacBio_data/Processing/mock/ASV_usearch/unfiltered/filter_default/zotus_minsize2.fa --chimeras chimera_unfiltered_minsize2.fa  --nonchimeras nonchimera_unfiltered_minsize2.fa --db /user_data/cko/Projects/P004/P004_J002_Mock/Processing_canu/trimmed/EUK_mock_ref_210123.fa --threads 10 
echo
echo
echo uchime_ref_default_unfiltered_minsize3
vsearch --uchime_ref /user_data/cko/Projects/P004/PacBio_data/Processing/mock/ASV_usearch/unfiltered/filter_default/zotus_minsize3.fa --chimeras chimera_unfiltered_minsize3.fa  --nonchimeras nonchimera_unfiltered_minsize3.fa --db /user_data/cko/Projects/P004/P004_J002_Mock/Processing_canu/trimmed/EUK_mock_ref_210123.fa --threads 10
echo
echo
echo uchime_ref_default_unfiltered_minsize4
vsearch --uchime_ref /user_data/cko/Projects/P004/PacBio_data/Processing/mock/ASV_usearch/unfiltered/filter_default/zotus_minsize4.fa --chimeras chimera_unfiltered_minsize4.fa  --nonchimeras nonchimera_unfiltered_minsize4.fa --db /user_data/cko/Projects/P004/P004_J002_Mock/Processing_canu/trimmed/EUK_mock_ref_210123.fa --threads 10
echo
echo
echo uchime_ref_default_unfiltered_minsize5
vsearch --uchime_ref /user_data/cko/Projects/P004/PacBio_data/Processing/mock/ASV_usearch/unfiltered/filter_default/zotus_minsize5.fa --chimeras chimera_unfiltered_minsize5.fa  --nonchimeras nonchimera_unfiltered_minsize5.fa --db /user_data/cko/Projects/P004/P004_J002_Mock/Processing_canu/trimmed/EUK_mock_ref_210123.fa --threads 10
echo
echo
echo uchime_ref_default_unfiltered_minsize6
vsearch --uchime_ref /user_data/cko/Projects/P004/PacBio_data/Processing/mock/ASV_usearch/unfiltered/filter_default/zotus_minsize6.fa --chimeras chimera_unfiltered_minsize6.fa  --nonchimeras nonchimera_unfiltered_minsize6.fa --db /user_data/cko/Projects/P004/P004_J002_Mock/Processing_canu/trimmed/EUK_mock_ref_210123.fa --threads 10
echo
echo
echo uchime_ref_default_unfiltered_minsize7
vsearch --uchime_ref /user_data/cko/Projects/P004/PacBio_data/Processing/mock/ASV_usearch/unfiltered/filter_default/zotus_minsize7.fa --chimeras chimera_unfiltered_minsize7.fa  --nonchimeras nonchimera_unfiltered_minsize7.fa --db /user_data/cko/Projects/P004/P004_J002_Mock/Processing_canu/trimmed/EUK_mock_ref_210123.fa --threads 10
echo
echo
echo uchime_ref_default_unfiltered_minsize8
vsearch --uchime_ref /user_data/cko/Projects/P004/PacBio_data/Processing/mock/ASV_usearch/unfiltered/filter_default/zotus_minsize8.fa --chimeras chimera_unfiltered_minsize8.fa  --nonchimeras nonchimera_unfiltered_minsize8.fa --db /user_data/cko/Projects/P004/P004_J002_Mock/Processing_canu/trimmed/EUK_mock_ref_210123.fa --threads 10
echo

echo
echo uchime_ref_selfid_unfiltered_minsize2
vsearch --uchime_ref /user_data/cko/Projects/P004/PacBio_data/Processing/mock/ASV_usearch/unfiltered/filter_default/zotus_minsize2.fa --chimeras chimera_unfiltered_minsize2_selfid.fa  --nonchimeras nonchimera_unfiltered_minsize2_selfid.fa --selfid --db /user_data/cko/Projects/P004/P004_J002_Mock/Processing_canu/trimmed/EUK_mock_ref_210123.fa --threads 10 
echo
echo
echo uchime_ref_selfid_unfiltered_minsize3
vsearch --uchime_ref /user_data/cko/Projects/P004/PacBio_data/Processing/mock/ASV_usearch/unfiltered/filter_default/zotus_minsize3.fa --chimeras chimera_unfiltered_minsize3_selfid.fa --nonchimeras nonchimera_unfiltered_minsize3_selfid.fa --selfid --db /user_data/cko/Projects/P004/P004_J002_Mock/Processing_canu/trimmed/EUK_mock_ref_210123.fa --threads 10
echo
echo
echo uchime_ref_selfid_unfiltered_minsize4
vsearch --uchime_ref /user_data/cko/Projects/P004/PacBio_data/Processing/mock/ASV_usearch/unfiltered/filter_default/zotus_minsize4.fa --chimeras chimera_unfiltered_minsize4_selfid.fa  --nonchimeras nonchimera_unfiltered_minsize4_selfid.fa --selfid --db /user_data/cko/Projects/P004/P004_J002_Mock/Processing_canu/trimmed/EUK_mock_ref_210123.fa --threads 10
echo
echo
echo uchime_ref_selfid_unfiltered_minsize5
vsearch --uchime_ref /user_data/cko/Projects/P004/PacBio_data/Processing/mock/ASV_usearch/unfiltered/filter_default/zotus_minsize5.fa --chimeras chimera_unfiltered_minsize5_selfid.fa  --nonchimeras nonchimera_unfiltered_minsize5_selfid.fa --selfid --db /user_data/cko/Projects/P004/P004_J002_Mock/Processing_canu/trimmed/EUK_mock_ref_210123.fa --threads 10
echo
echo
echo uchime_ref_selfid_unfiltered_minsize6
vsearch --uchime_ref /user_data/cko/Projects/P004/PacBio_data/Processing/mock/ASV_usearch/unfiltered/filter_default/zotus_minsize6.fa --chimeras chimera_unfiltered_minsize6_selfid.fa  --nonchimeras nonchimera_unfiltered_minsize6_selfid.fa --selfid --db /user_data/cko/Projects/P004/P004_J002_Mock/Processing_canu/trimmed/EUK_mock_ref_210123.fa --threads 10
echo
echo
echo uchime_ref_selfid_unfiltered_minsize7
vsearch --uchime_ref /user_data/cko/Projects/P004/PacBio_data/Processing/mock/ASV_usearch/unfiltered/filter_default/zotus_minsize7.fa --chimeras chimera_unfiltered_minsize7_selfid.fa  --nonchimeras nonchimera_unfiltered_minsize7_selfid.fa --selfid --db /user_data/cko/Projects/P004/P004_J002_Mock/Processing_canu/trimmed/EUK_mock_ref_210123.fa --threads 10
echo
echo
echo uchime_ref_selfid_unfiltered_minsize8
vsearch --uchime_ref /user_data/cko/Projects/P004/PacBio_data/Processing/mock/ASV_usearch/unfiltered/filter_default/zotus_minsize8.fa --chimeras chimera_unfiltered_minsize8_selfid.fa  --nonchimeras nonchimera_unfiltered_minsize8_selfid.fa --selfid --db /user_data/cko/Projects/P004/P004_J002_Mock/Processing_canu/trimmed/EUK_mock_ref_210123.fa --threads 10
echo