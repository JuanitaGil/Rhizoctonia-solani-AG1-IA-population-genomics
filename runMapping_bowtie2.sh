#!/bin/bash
p=$1;
x=1000;
threads=16

 # input files
f1=${p}_1.fastq.gz;
f2=${p}_2.gz;
REFERENCE=ref.fna;

 # software variables. Write paths only if you can not install the programs or can not use installed versions
BOWTIE2=/path/to/bowtie2-2.3.4.1-linux-x86_64/bowtie2;
SAMTOOLS=/path/to/samtools
JAVA="/path/to/java -d64 -XX:MaxHeapSize=1g";

 # jars for java packages
PICARD=/path/to/picard.jar;

 # map the reads and sort the alignment
mkdir ${p}_tmpdir;
bowtie2 --rg-id ${p} --rg SM:${p} --rg PL:ILLUMINA -p ${threads} -X ${x} -k 3 -t -x ${REFERENCE} -1 ${f1} -2 ${f2} 2> ${p}_bowtie2.log | java -Xmx8g -jar ${PICARD} SortSam MAX_RECORDS_IN_RAM=1000000 SO=coordinate CREATE_INDEX=true TMP_DIR=${p}_tmpdir I=/dev/stdin O=${p}_bowtie2_sorted.bam > ${p}_bowtie2_sort.log 2>&1;
rm -rf ${p}_tmpdir;
