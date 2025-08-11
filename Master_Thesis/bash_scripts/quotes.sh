#!/bin/bash

# Input TSV file
input_tsv="/data/leuven/350/vsc35091/THESIS/DATA/sample_IDs_fixed.tsv"

# Output TSV file
output_tsv="/data/leuven/350/vsc35091/THESIS/DATA/sample_IDs_categorical.tsv"

# Create the header in the output TSV file
head -n 1 "$input_tsv" > "$output_tsv"

# Process the rest of the lines, enclosing the trimmed third column in double quotes
awk -F'\t' '{OFS="\t"; gsub("\n", "", $3); $3 = "\"" $3 "\""; $0=$1 FS $2 FS $3} NR>1' "$input_tsv" >> "$output_tsv"

echo "TSV file with corrected quotation marks created: $output_tsv"