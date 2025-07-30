#!/bin/bash

#Run for each chromosome
CHR=chr1;
VCF=/path/to/targetPop_${CHR}_phased.vcf.gz;
REFVCF=/path/to/refPop_${CHR}_phased.vcf.gz;
MAP=/path/to/targetPop_${CHR}_phased.maps;
SELSCAN=/path/to/selscan

${SELSCAN} --xpehh --vcf ${VCF} --vcf-ref ${REFVCF} --map ${MAP} --out pop_${CHR}_phased --threads 25 --cutoff 0.20 --pmap > ${CHR}_selscan.log 2>&1
