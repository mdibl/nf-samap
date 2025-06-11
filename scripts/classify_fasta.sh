#!/usr/bin/env bash

# Check for argument
if [ $# -ne 1 ]; then
    echo "Usage: $0 <fasta_file>"
    exit 1
fi
fasta_file="$1"

# Grab first non-header sequence line
line=$(grep -v '^>' "$fasta_file" | head -n 1)

# If it contains only DNA letters, it's probably nucleotide
if [[ "$line" =~ ^[ACGTUNacgtun]+$ ]]; then
    echo "nucl"
else
    echo "prot"
fi