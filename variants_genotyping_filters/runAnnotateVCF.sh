#!/bin/bash

 # input files
VCF=/path/to/filtered.vcf.gz;
REFERENCE=/path/to/ref.fna;
GFF=/path/to/ref.gff3;

 # software variables. Write paths only if you can not install the programs or can not use installed versions
#BOWTIE2=/path/to/bowtie2-2.3.4.1-linux-x86_64/bowtie2;
#SAMTOOLS=/path/to/samtools
#JAVA="/path/to/java -d64 -XX:MaxHeapSize=1g";

 # jars for java packages
NGSEP=/path/to/NGSEPcore_4.3.1.jar;

 # call variants
java -Xmx32g -jar ${NGSEP} VCFAnnotate -i ${VCF} -r ${REFERENCE} -t ${GFF} -o output_annotated.vcf > out.log 2>&1;
