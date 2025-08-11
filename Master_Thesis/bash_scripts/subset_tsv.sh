#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 file1.tsv file2.tsv"
    exit 1
fi

file1="$1"
file2="$2"

# Check if the files exist
if [ ! -f "$file1" ] || [ ! -f "$file2" ]; then
    echo "Error: Input files not found."
    exit 1
fi

# Use awk to compare and filter lines
awk -F'\t' 'NR==FNR { subjects[$1]; next } $1 in subjects' "$file1" "$file2" > temp_file

# Replace the original file2 with the temp_file
mv temp_file "$file2"

echo "Comparison and filtering complete."
