#!/bin/bash

#folder and file extension
folder="/mnt/DataSSD/Rhizoctonia_pop_AG1/wdirectory/variants/HG81_reference/bcftools/list_of_genes"
extension=".txt"
ref=/mnt/DataSSD/Rhizoctonia_pop_AG1/wdirectory/references/AG1-IA_HG81/Rsolani_AG1-IA_HG81.fna

for file in $folder/*$extension; do
echo "Processing file: $file"
python3 vcf2fastaCDS_rhizoctonia.py --listvcf ${file} --path ./CDS/ --ref ${ref} --out ./CDS_fastas/ > vcf2fastaCDS.log 2>&1
done;
