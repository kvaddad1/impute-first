#!/bin/bash

# Usage: ./create_SDF.sh <path_to_reference_sequence>

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path_to_reference_sequence>"
    exit 1
fi

REFERENCE=$1
OUTPUT_DIR="${REFERENCE}_SDF"

# Path to RTG Tools
RTG_PATH="RTG.jar"

# Run RTG format to create SDF
java -jar $RTG_PATH format -o "$OUTPUT_DIR" "$REFERENCE"

echo "SDF created at $OUTPUT_DIR"

