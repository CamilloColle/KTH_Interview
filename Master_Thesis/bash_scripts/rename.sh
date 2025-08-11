#!/bin/bash

# Input folder containing files
input_folder="/data/leuven/350/vsc35091/THESIS/DATA/fastq/run4"

# Change to the input folder
cd "$input_folder" || exit

# Rename files matching the pattern X--Y.fastq.gz to X----Y.fastq.gz
for file in *--*.fastq.gz; do
  new_name=$(echo "$file" | sed 's/--/----/')
  mv "$file" "$new_name"
  echo "Renamed: $file to $new_name"
done

# Print a message
echo "File renaming complete."
