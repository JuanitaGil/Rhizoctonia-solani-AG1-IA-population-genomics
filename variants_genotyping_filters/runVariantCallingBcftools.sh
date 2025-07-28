#!/bin/bash

 # input files
BAMS=/path/to/bam_files;
REFERENCE=/path/to/ref.fna;

# output variable
OUT_PREFIX=

 # call variants
bcftools mpileup --threads 20 -a DP,AD,INFO/AD -f ${REFERENCE} ${BAMS}/*bam | bcftools call --threads 20 -m --ploidy 2 --output-type z -o ${OUT_PREFIX}.vcf.gz > bcftools.log 2>&1
