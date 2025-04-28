#!/bin/bash

# Original script file
original_script="test.sh"

# Loop to create modified copies
for i in {1..50}; do
    # Define output filename
    output_script="runpfas${i}.sh"

    # Replace
    sed "s/test.log/pfas${i}.log/g" "$original_script" > "$output_script"

    # Make the new script executable (optional)
    chmod +x "$output_script"

    echo "Created $output_script"
done

echo "All scripts have been generated."

