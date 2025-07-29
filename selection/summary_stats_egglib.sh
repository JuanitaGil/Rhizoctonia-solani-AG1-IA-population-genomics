#!/bin/bash

#folder and file extension
folder="path/to/CDS"
extension=".fasta"

for file in $folder/*$extension; do
echo "Processing file: $file"

# Extracting the filename without extension
filename=$(basename "$file")
filename_no_ext="${filename%.*}"

# Output file name
output_file="${filename_no_ext}_egglib.stats"

python3 summary_stats_egglib.py --path ${file} --min_size 10 --min_length 400 --max_missing 0.50 --pathout $output_file
done;
