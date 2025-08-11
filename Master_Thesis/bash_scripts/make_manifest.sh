#!/bin/bash

# Directory containing the fastq.gz files
fastq_dir="/data/leuven/350/vsc35091/THESIS/DATA/000_fastq"

# Metadata file
metadata_file="/data/leuven/350/vsc35091/THESIS/DATA/input_files/sample_IDs_fixed.tsv"

# Output TSV file
output_file="manifest33V2_new.tsv"

# Create the header in the output TSV file
echo -e "sample-id\tabsolute-filepath" > "$output_file"

# Loop through the fastq.gz files in the directory
for fastq_file in "$fastq_dir"/*.fastq.gz; do
  if [ -f "$fastq_file" ]; then
    # Extract the barcode from the filename
    barcode=$(echo "$fastq_file" | awk -F'.' '{print $1}' | awk -F'/' '{print $9}')
    
    
    # Look up the sample-id in the metadata file based on the barcode
    sample_id=$(grep "$barcode" "$metadata_file" | cut -f 1)
    echo "$barcode, $sample_id"
    #echo "$sample_id"
    # Get the absolute filepath of the fastq.gz file
    absolute_filepath=$(realpath "$fastq_file")
    
    # Append the data to the output TSV file
    #echo -e "$sample_id\t$absolute_filepath"
    
    echo -e "$sample_id\t$absolute_filepath" >> "$output_file"
  fi
done

echo "TSV file created: $output_file"

