#!/bin/bash

longread_umi nanopore_pipeline -d rawreads.fastq -v 2 -o nanopore_pipeline -s 90 -e 73 -m 3000 -M 6000 -t 60 -f TTTCTGTTGGTGCTGATATTGC -F GGCAAGTCTGGTGCCAG -r ACTTGCCTGTCGCTCTATCTTC -R GACGAGGCATTTGGCTACCTT -c 3 -p 2 -u nanopore_pipeline/umi_binning -t 60 -T 2