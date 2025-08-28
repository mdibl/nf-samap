#!/bin/bash
#==========================================================
# Author: Ryan Sonderman
# Date: 2025-06-17
# Version: 1.1.0
# Description: This script adds a 2-char id and classifies 
#              fasta files for a sample sheet
# Dependencies: None
# Usage: ./update_sample_sheet.sh sample_sheet outfile
#==========================================================

# Standardized logging
log() {
    local level="$1"
    local message="$2"
    local timestamp
    timestamp=$(date +"%Y-%m-%d %H:%M:%S.%3N")
    echo "$timestamp [$level]: $message"
}


#==========================================================
# Variable Initialization
#==========================================================
input_csv=""
output_csv=""
hashed_id_len=2
declare -A used_ids

log "INFO" "Beginning execution of script"

# Parse command-line arguments
if [ $# -ne 2 ]; then
    log "ERROR" "Usage: $0 <sample_sheet> <outfile>"
    exit 1
fi
input_csv="$1"
if [ ! -f "$input_csv" ]; then
    log "ERROR" "FIle '$input_csv' not found or not readable"
    exit 2
fi
log "INFO" "Found sample sheet: $input_csv"
output_csv="$2"
log "INFO" "Output will be stored in: $output_csv"

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
log "INFO" "Adding new headers for 'type' and 'id2'"
header=$(head -n 1 "$input_csv")
echo "${header},type,id2" > "$output_csv"

# Process each row in the csv
log "INFO" "Updating sample sheet with type and id2 values"
while IFS=, read -r id h5ad fasta annotation; do
    log "INFO" "> Updating entry $id - ($fasta)"
    log "INFO" "  > Classifying fasta"
    type=$(classify_fasta "$fasta")
    log "INFO" "  > Fasta classified as '$type'"
    
    # Generate a reproducible key by hashing a combination of id and fasta
    log "INFO" "  > Generating hash"
    hash=$(echo -n "${id}${fasta}" | sha256sum | awk '{print $1}' | tr '[:upper:]' '[:lower:]')
    # Take the first two characters of the hash, ensuring they are valid lowercase alphabetic characters
    key=$(echo "$hash" | tr -cd 'a-z' | cut -c1-"$hashed_id_len")
    log "INFO" "  > Trying '$key'"

    # Ensure unique key generation with max 100 attempts
    attempt=0
    max_attempts=100
    while [[ -n "${used_ids[$key]}" || ${#key} -lt $hashed_id_len ]]; do
        log "WARNING" "  > $key not valid"
        # If key is not long enough or already seen, shift and regenerate
        hash=$(echo "$hash" | tr 'a-z' 'b-za')
        key=$(echo "$hash" | tr -cd 'a-z' | cut -c1-"$hashed_id_len")
        log "INFO" "  > Trying '$key'"
        attempt=$((attempt + 1))
        if [[ $attempt -ge $max_attempts ]]; then
            log "ERROR" "Exceeded attemps to generate unique id2"
            exit 4
        fi
    done

    log "INFO" "  > '$key' is valid"
    used_ids[$key]=1    # Store the unique key in the seen_keys array
    id2="$key"          # Assign the unique key to id2

    echo "${id},${h5ad},${fasta},${annotation},${type},${id2}" >> "$output_csv"
done < <(tail -n +2 "$input_csv")

log "INFO" "Script complete, see $output_csv"
