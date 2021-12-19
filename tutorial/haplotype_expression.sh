#!/bin/bash

python ../../phaser/phaser_gene_ae/phaser_gene_ae.py \
    --haplotypic_counts results/phaser_test_case.haplotypic_counts.txt \
    --features references/gencode.v19.GRCh37.genes.bed \
    --o results/phaser_test_case_gene_ae.txt
