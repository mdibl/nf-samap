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

# Function to create id key with given length
generate_id() {
    local valid_chars=({a..z}) # Define usable characters
    local key_length=$1 # Default key length
    local key=""

    for ((i=0; i<key_length; i++)); do
        key+="${valid_chars[RANDOM % ${#valid_chars[@]}]}" # Randomly select a character
    done

    echo "$key"
}


# Check for arguments
if [ $# -ne 1 ]; then
    echo "Usage: $0 <sample_sheet>"
    exit 1
fi

# Input and output file paths
input_csv="$1"
timestamp=$(date '+%Y%m%d_%H%M%S')
output_csv="${input_csv%.csv}_${timestamp}.csv"

# Read the header, append new columns, and write to output
header=$(head -n 1 "$input_csv")
echo "${header},type,id2" > "$output_csv"

# Params for generating id2
declare -A seen_keys
desired_length=2 # Set desired key length

# Process each row, including the last line even if no newline at EOF
tail -n +2 "$input_csv" | while IFS=, read -r id h5ad fasta annotation; do
    type=$(classify_fasta "$fasta")
    
    # Generate a reproducible key by hashing a combination of id and fasta
    hash=$(echo -n "${id}${fasta}" | sha256sum | awk '{print $1}' | tr '[:upper:]' '[:lower:]')

    # Take the first two characters of the hash, ensuring they are valid alphabetic characters
    key=$(echo "$hash" | cut -c1-"$desired_length")

    # Ensure unique key generation with max 100 attempts
    attempt=0
    max_attempts=100
    while [[ -n "${seen_keys[$key]}" ]]; do
        key=$(echo "$key" | tr 'a-z' 'b-za')  # Modify the key if it exists (shifting characters)
        attempt=$((attempt + 1))
        if [[ $attempt -ge $max_attempts ]]; then
            echo "Exceeded attempts to generate unique id2 for ${id} (${fasta})" >&2
            exit 1
        fi
    done

    seen_keys[$key]=1  # Store the unique key in the seen_keys array
    id2="$key"         # Assign the unique key to id2

    echo "${id},${h5ad},${fasta},${annotation},${type},${id2}" >> "$output_csv"
done < <(tail -n +2 "$input_csv")

