#!/bin/bash

# Argument1:  BlobDir directory 
# Argument2:  ID to name output files (ToLID or GCA accession)
# Argument3:  coverage data available (string: "TRUE" or "FALSE")
# Argument4:  args = other parameters required to run the module

BLOBDIR=$1
ID=$2
COV=$3
ARGS=$4

# process depending on availability of coverage data
if [[ "$COV" == "TRUE" ]]; then
    blobtools view --view blob --param plotShape=circle --param largeFonts=true --format png --out $BLOBDIR "$ID" $ARGS
    blobtools view --view blob --param plotShape=hex --param largeFonts=true --format png --out $BLOBDIR "$ID" $ARGS
    blobtools view --view blob --param plotShape=square --param largeFonts=true --format png --out $BLOBDIR "$ID" $ARGS
    blobtools view --view blob --param plotShape=kite --param largeFonts=true --format png --out $BLOBDIR "$ID" $ARGS
    blobtools view --view cumulative --param largeFonts=true --format png --out $BLOBDIR "$ID" $ARGS
    blobtools view --view snail --param largeFonts=true --format png --out $BLOBDIR "$ID" $ARGS
    mv $BLOBDIR/*.png ./
elif [[ "$COV" == "FALSE" ]]; then
    blobtools view --view cumulative --param largeFonts=true --format png --out $BLOBDIR "$ID" $ARGS
    blobtools view --view snail --param largeFonts=true --format png --out $BLOBDIR "$ID" $ARGS
    mv $BLOBDIR/*.png ./
else
    echo "ERROR: $COV is not a valid value, please use TRUE or FALSE"
    exit 1
fi
