#!/usr/bin/env bash
input=${1}
output=${2}

samtools view "$input" |\
  parallel --pipe --progress --block 100M \
  awk \'{print \$1,\$16}\' \
  > "$output"