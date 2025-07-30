#!/bin/bash

beagle=/path/to/beagle.22Jul22.46e.jar;
t=16;
pop=AR;

for chr in $(cat chromosomes.txt)
do
java -Xmx16g -jar ${beagle} gt=${pop}_${chr}.vcf out=${pop}_${chr}_phased impute=false nthreads=${t} > beagle_${pop}_${chr}.log 2>&1
done;
