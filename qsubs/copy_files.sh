#!/bin/bash

# Destination directory where you want to copy the files
destination_dir=$1


# Check if the file list exists
if [ -f "BRCA1/BRCA1_annotated_BreastCancer_files_search_1.txt" ]; then
  # Read each file path from the file_list.txt and copy them to the destination directory
  while IFS= read -r file; do
    # Check if the file exists before copying
    if [ -e "$file" ]; then
      # Use the "basename" command to extract the filename without the path
      filename=$(basename "$file")
      # Copy the file to the destination directory
      cp "$file" "$destination_dir/$filename"
      echo "Copied $file to $destination_dir/$filename"
    else
      echo "File $file not found."
    fi
  done < "BRCA1/BRCA1_annotated_BreastCancer_files_search_1.txt"
else
  echo "File list file (file_list.txt) not found."
fi
