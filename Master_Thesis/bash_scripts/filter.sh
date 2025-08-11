#!/bin/bash

# Source directory containing fastq files
source_dir="/Users/camillocolleluori/Desktop/KUL/Thesis/last_sequencing/fastq_results"

# TSV file with IDs in the first column
tsv_file="/Users/camillocolleluori/Desktop/KUL/Thesis/last_sequencing/samples_run4.tsv"

# Destination directory
destination_dir="/Volumes/MaxOne/DATA/fastq/run4"


# Loop through each line in the TSV file
while IFS=$'\t' read -r id rest; do
  # Check if the ID is not empty (skip header lines)
  if [ -n "$id" ]; then
    # Construct the source file path
    source_file="$source_dir/$id.fastq.gz"
    
    # Check if the source file exists
    if [ -e "$source_file" ]; then
      # Move the file to the destination directory
      mv "$source_file" "$destination_dir/"
      echo "Moved: $source_file to $destination_dir/"
    else
      echo "File not found: $source_file"
    fi
  fi
done < "$tsv_file"

echo "File movement complete."

