#!/bin/bash

# Input folder containing the fastq.gz files
input_folder="/data/leuven/350/vsc35091/THESIS/DATA/000_fastq"

# Output folder for trimmed files
output_folder="/data/leuven/350/vsc35091/THESIS/DATA/trimmed_fastq"

# Create the output folder if it doesn't exist
mkdir -p "$output_folder"

# Iterate over each fastq.gz file in the input folder
for file in "$input_folder"/*.fastq.gz; do
    # Extract the file name (excluding the path and extension)
    filename=$(basename "$file")
    filename_no_ext="${filename%.*}"

    # Run bbduk.sh for each file
    bbduk.sh in="$file" out="$output_folder/$filename_no_ext.trimmed.fastq.gz" qtrim=rl trimq=10 minlength=50
done
