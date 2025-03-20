#!/bin/bash

# Define the files to keep
KEEP_FILES=("ana.dat" "ccdcif_info.log" "clean.sh")

# Remove all files except the ones in KEEP_FILES
for file in *; do
    if [[ ! " ${KEEP_FILES[@]} " =~ " $file " ]]; then
        rm -f "$file"
        echo "Deleted: $file"
    fi
done

echo "Cleanup complete. Only ${KEEP_FILES[@]} remain."

