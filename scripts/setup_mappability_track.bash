#!/bin/bash

# Setup script to download and process mappability tracks
# Usage: bash setup_mappability_track.bash hg19 [output_directory]

# Get the directory of this script
SCRIPT_DIR=$(dirname "$(realpath "$0")")

# Validate input
if [ -z "$1" ]; then
    echo "Error: No genome provided. Usage: bash setup_mappability_track.bash hg19 [output_directory]"
    exit 1
fi

GENOME=$1
OUTPUT_DIR=${2:-"$SCRIPT_DIR/../resources"}  # Default to 'resources' folder in the repo root's directory

if [[ "$GENOME" != "hg19" && "$GENOME" != "GRCh37" ]]; then
    echo "Error: Only 'hg19' or 'GRCh37' are currently implemented."
    echo "Please raise an issue on [github link] if you need another genome."
    exit 1
fi

# Check for bigWigToBedGraph
if ! command -v bigWigToBedGraph &> /dev/null; then
    echo "Error: 'bigWigToBedGraph' is not installed or not found in your PATH."
    echo "Please ensure you have activated the correct environment or created one using 'requirements.txt'."
    exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR" || exit

# Download file
URL="http://hgdownload.cse.ucsc.edu/gbdb/hg19/bbi/wgEncodeCrgMapabilityAlign50mer.bw"
FILE="wgEncodeCrgMapabilityAlign50mer.bw"

echo "Downloading reference file for $GENOME into directory: $OUTPUT_DIR"
curl -O "$URL" || { echo "Download failed!"; exit 1; }

# Convert to BedGraph
echo "Converting to BedGraph. This may take a few minutes..."
bigWigToBedGraph "$FILE" "${GENOME}_mappability.bedGraph" || { echo "Conversion failed!"; exit 1; }

# Delete original file
echo "Cleaning up..."
rm "$FILE"

echo "Setup complete! Processed file saved to: $OUTPUT_DIR/${GENOME}_mappability.bedGraph"
