#!/bin/bash

# Function to classify fasta as nucl or prot
classify_fasta() {
    local fasta_file="$1"
    local line
    line=$(grep -v '^>' "$fasta_file" | head -n 1)
    if [[ "$line" =~ ^[ACGTUNacgtun]+$ ]]; then
        echo "nucl"
    else
        echo "prot"
    fi
}

create_id2() {
    # Dummy id2 function, replace with your logic if needed
    echo "$1"
}

# Check for arguments
if [ $# -ne 2 ]; then
    echo "Usage: $0 <sample_sheet> <output_file>"
    exit 1
fi

input_csv="$1"
output_csv="$2"

# Read the header, append new columns, and write to output
header=$(head -n 1 "$input_csv")
echo "${header},type,id2" > "$output_csv"

# Process each row, including the last line even if no newline at EOF
tail -n +2 "$input_csv" | while IFS=, read -r id h5ad fasta annotation || [ -n "$id$h5ad$fasta$annotation" ]; do
    type=$(classify_fasta "$fasta")
    id2=$(create_id2 "$id")
    echo "${id},${h5ad},${fasta},${annotation},${type},${id2}" >> "$output_csv"
done

