#!/bin/bash

python ../phaser/phaser.py \
    --vcf data/NA06986.vcf.gz \
    --bam data/NA06986.2.M_111215_4.bam \
    --paired_end 1 \
    --mapq 255 \
    --baseq 10 \
    --sample NA06986 \
    --blacklist references/hg19_hla.bed \
    --haplo_count_blacklist references/hg19_haplo_count_blacklist.bed \
    --threads 4 \
    --o results/phaser_test_case
