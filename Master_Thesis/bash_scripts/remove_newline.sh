#!/bin/bash

# Input TSV file
input_tsv="/data/leuven/350/vsc35091/THESIS/DATA/sample_IDs_categorical.tsv"

# Output TSV file
output_tsv="/data/leuven/350/vsc35091/THESIS/DATA/sample_IDs_categorical2.tsv"

# Create the header in the output TSV file
head -n 1 "$input_tsv" > "$output_tsv"

# Process the rest of the lines, removing "\n" from every other line starting from the second line
awk -F'\t' 'NR>1 && NR%2==0 {gsub(/\\n/, "", $0)} {print}' "$input_tsv" >> "$output_tsv"

echo "TSV file with removed newline characters created: $output_tsv"
