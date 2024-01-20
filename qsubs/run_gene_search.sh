#!/bin/bash

# Check if gene list file is provided
if [ -z "$1" ]; then
    echo "Usage: $0 [path to gene list file]"
    exit 1
fi

# File containing list of genes
GENE_LIST_FILE=$1

# Loop through each line in the file
while IFS= read -r gene; do
    ./gene_search.sh "$gene"
done < "$GENE_LIST_FILE"

