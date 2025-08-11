#!/bin/bash

# Source TSV file
source_tsv="/data/leuven/350/vsc35091/THESIS/DATA/sample_IDs_fixed.tsv"

# Destination TSV file
output_tsv="unmatched_filepaths.tsv"

# Directory containing fastq.gz files
fastq_dir="/data/leuven/350/vsc35091/THESIS/DATA/fastq"

# Create the header in the output TSV file
echo -e "absolute-filepath" > "$output_tsv"

# Loop through the fastq.gz files in the directory
for fastq_file in "$fastq_dir"/*.fastq.gz; do
  if [ -f "$fastq_file" ]; then
    # Extract the barcode from the filename
    barcode=$(basename "$fastq_file" | awk -F'.' '{print $1}')
    
    # Check if the barcode is not present in the source TSV file
    if ! grep -q -w "$barcode" "$source_tsv"; then
      # Get the absolute filepath of the unmatched file
      absolute_filepath=$(realpath "$fastq_file")
      
      # Append the filepath to the output TSV file
      echo -e "$absolute_filepath" >> "$output_tsv"
    fi
  fi
done

echo "Unmatched Filepaths TSV file created: $output_tsv"
