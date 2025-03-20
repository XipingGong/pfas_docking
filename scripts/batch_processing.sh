#!/bin/bash

# Function to display help message
usage() {
    echo "Usage: $0 script.py|script.sh input.txt [--options]"
    echo "Runs the given Python or Bash script for each line in the input text file."
    echo "Lines starting with '#' are treated as comments and printed as is."
    echo "Additional options will be passed to the script."
    exit 1
}

# Check if -h option is passed
if [[ "$1" == "-h" ]]; then
    usage
fi

# Check for minimum number of arguments
if [[ $# -lt 2 ]]; then
    usage
fi

SCRIPT_FILE="$1"
INPUT_FILE="$2"
shift 2 # Shift arguments to access additional options
EXTRA_OPTIONS="$@"

# Check if the script file exists
if [[ ! -f "$SCRIPT_FILE" ]]; then
    echo "Error: Script file '$SCRIPT_FILE' not found."
    exit 1
fi

# Check if the input file exists
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file '$INPUT_FILE' not found."
    exit 1
fi

# Determine script type based on extension
SCRIPT_TYPE=""
if [[ "$SCRIPT_FILE" == *.py ]]; then
    SCRIPT_TYPE="python"
elif [[ "$SCRIPT_FILE" == *.sh ]]; then
    SCRIPT_TYPE="bash"
else
    echo "Error: Unsupported script type. Please provide a Python (.py) or Bash (.sh) script."
    exit 1
fi

# Ensure Python or Bash is installed
if [[ "$SCRIPT_TYPE" == "python" && ! $(command -v python) ]]; then
    echo "Error: Python is not installed or not in PATH."
    exit 1
elif [[ "$SCRIPT_TYPE" == "bash" && ! $(command -v bash) ]]; then
    echo "Error: Bash is not installed or not in PATH."
    exit 1
fi

# Start total execution timer
TOTAL_START_TIME=$(date +%s.%N)

# Read the input file line by line
while IFS= read -r line; do
    # If line starts with '#', print it as a comment
    if [[ "$line" =~ ^# ]]; then
        echo "$line"
    else
        # Construct command based on script type
        if [[ "$SCRIPT_TYPE" == "python" ]]; then
            COMMAND="python \"$SCRIPT_FILE\" $line $EXTRA_OPTIONS"
        else
            COMMAND="bash \"$SCRIPT_FILE\" $line $EXTRA_OPTIONS"
        fi

        echo "Executing: $COMMAND"
        eval "$COMMAND"
        echo ""
    fi
done < "$INPUT_FILE"

# End total execution timer
TOTAL_END_TIME=$(date +%s.%N)
TOTAL_ELAPSED_TIME=$(echo "$TOTAL_END_TIME - $TOTAL_START_TIME" | bc)
echo "Total execution time: $TOTAL_ELAPSED_TIME seconds"

