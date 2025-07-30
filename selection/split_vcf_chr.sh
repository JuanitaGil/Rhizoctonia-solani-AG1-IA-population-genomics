#!/bin/bash

for i in $(cat chromosomes.txt)
do
vcftools --gzvcf pop.vcf --chr ${i} --out pop_${i} --recode --recode-INFO-all
done;
