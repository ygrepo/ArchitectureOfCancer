#!/bin/bash

# Check if gene list file is provided
if [ -z "$1" ]; then
    echo "Usage: $0 [path to gene list file]"
    exit 1
fi
 
# File containing list of genes
GENE_LIST_FILE=$1

# Check if the gene list file exists
if [ ! -f "$GENE_LIST_FILE" ]; then
    echo "Gene list file not found: $GENE_LIST_FILE"
    exit 1
fi

# Loop through each line in the file
while IFS= read -r GENE_NAME; do
    # Create a directory for each gene, if it doesn't exist
    mkdir -p "data/genes/${GENE_NAME}"
done < "$GENE_LIST_FILE"