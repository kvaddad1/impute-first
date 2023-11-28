#!/bin/bash

# Usage: prepare_mapfiles.sh <linkage_map_dir>
# Example: prepare_mapfiles.sh linkage_maps
# Note: This script expects PLINK format files.

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <linkage_map_dir>"
    exit 1
fi

linkage_map_dir=$1

for i in $(seq 1 1 22); do
    echo "pos chr cM" > GLIMPSE_plink.chr${i}.GRCh38.map
    awk '{print $4,$1,$3}' ${linkage_map_dir}/plink.chr${i}.GRCh38.map >> GLIMPSE_plink.chr${i}.GRCh38.map
done

echo "Done prepping linkage map files for GLIMPSE"

