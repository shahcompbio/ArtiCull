#!/bin/bash

# Setup script to download and process mappability tracks
# Usage: bash setup_mappability_track.bash [output_directory]

# Get the directory of this script
SCRIPT_DIR=$(dirname "$(realpath "$0")")

GENOME=hg19
OUTPUT_DIR=${1:-"$SCRIPT_DIR/../resources"}  # Default to 'resources' folder in the repo root's directory

# Check for bigWigToBedGraph
if ! command -v bigWigToBedGraph &> /dev/null; then
    echo "Error: 'bigWigToBedGraph' is not installed or not found in your PATH."
    echo "Please ensure you have activated the correct environment or created one using 'requirements.yml'."
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
