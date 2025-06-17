#!/bin/bash
#==========================================================
# Author: Ryan Sonderman
# Date: 2025-06-17
# Version: 1.0.0
# Description: This script adds a 2-char id and classifies 
#              fasta files for a sample sheet
# Dependencies: None
# Usage: ./update_sample_sheet.sh sample_sheet
#==========================================================

#==========================================================
# Variable Initialization
#==========================================================
input_csv=""
output_csv="out.csv"
hashed_id_len=2
declare -A used_ids

# Parse command-line arguments
if [ $# -ne 1 ]; then
    echo "Usage: $0 <sample_sheet>"
    exit 1
fi
input_csv="$1"
if [ ! -f "$input_csv" ]; then
    echo "Error: File '$input_csv' not found or not readable." >&2
    exit 1
fi

#==========================================================
# Helper Functions
#==========================================================

# Classifies a fasta file as either 'nucl' (nucleotide) or 'prot' (protein) based on the content
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

generate_id() {
    local valid_chars=({a..z}) 
    local key_length=$1
    local key=""
    for ((i=0; i<key_length; i++)); do
        key+="${valid_chars[RANDOM % ${#valid_chars[@]}]}" # Randomly select a character
    done
    echo "$key"
}

#==========================================================
# Main Logic
#==========================================================

# Add the new headers
header=$(head -n 1 "$input_csv")
echo "${header},type,id2" > "$output_csv"

# Process each row in the csv
tail -n +2 "$input_csv" | while IFS=, read -r id h5ad fasta annotation; do
    type=$(classify_fasta "$fasta")
    
    # Generate a reproducible key by hashing a combination of id and fasta
    hash=$(echo -n "${id}${fasta}" | sha256sum | awk '{print $1}' | tr '[:upper:]' '[:lower:]')
    # Take the first two characters of the hash, ensuring they are valid lowercase alphabetic characters
    key=$(echo "$hash" | tr -cd 'a-z' | cut -c1-"$hashed_id_len")

    # Ensure unique key generation with max 100 attempts
    attempt=0
    max_attempts=100
    while [[ -n "${used_ids[$key]}" || ${#key} -lt $hashed_id_len ]]; do
        # If key is not long enough or already seen, shift and regenerate
        hash=$(echo "$hash" | tr 'a-z' 'b-za')
        key=$(echo "$hash" | tr -cd 'a-z' | cut -c1-"$hashed_id_len")
        attempt=$((attempt + 1))
        if [[ $attempt -ge $max_attempts ]]; then
            echo "Exceeded attempts to generate unique id2 for ${id} (${fasta})" >&2
            exit 1
        fi
    done

    used_ids[$key]=1    # Store the unique key in the seen_keys array
    id2="$key"          # Assign the unique key to id2

    echo "${id},${h5ad},${fasta},${annotation},${type},${id2}" >> "$output_csv"
done < <(tail -n +2 "$input_csv")

