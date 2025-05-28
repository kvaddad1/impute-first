#!/bin/bash

# Create the target directory if it doesn't exist
mkdir -p ./exampleData

# Zenodo record ID
ZENODO_ID="14902947"

# Use the Zenodo REST API to get the record metadata
echo "Fetching metadata for Zenodo record $ZENODO_ID..."
RESPONSE=$(curl -s "https://zenodo.org/api/records/$ZENODO_ID")

# Extract download links for all files (using jq if available, otherwise with grep)
if command -v jq > /dev/null; then
    # Using jq for more robust parsing
    echo "Using jq to parse JSON response..."
    FILES=$(echo "$RESPONSE" | jq -r '.files[] | "\(.links.self) \(.key)"')
    
    # Download each file
    echo "$FILES" | while read -r URL FILENAME; do
        echo "Downloading $FILENAME to ./exampleData/"
        curl -L -o "./exampleData/$FILENAME" "$URL"
    done
else
    # Fallback to grep/sed if jq is not available (less robust)
    echo "jq not found, using grep/sed (basic parsing)..."
    FILES=$(echo "$RESPONSE" | grep -o '"download": "[^"]*' | sed 's/"download": "//')
    
    for URL in $FILES; do
        FILENAME=$(basename "$URL")
        echo "Downloading $FILENAME to ./exampleData/"
        curl -L -o "./exampleData/$FILENAME" "$URL"
    done
fi

echo "Download complete! Files are available in ./exampleData/"
