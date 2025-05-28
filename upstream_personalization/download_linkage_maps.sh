#!/bin/bash
set -x
# Create directories
mkdir -p linkage_maps/{beagle,glimpse}

# Download and extract the linkage maps
echo "Downloading linkage maps..."
wget "https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip" -O linkage_maps/plink.GRCh38.map.zip
unzip -o linkage_maps/plink.GRCh38.map.zip -d linkage_maps
rm linkage_maps/plink.GRCh38.map.zip

# Process maps for BEAGLE (adding 'chr' prefix)
echo "Processing maps for BEAGLE..."
for chr in {1..22}; do
    input_file="linkage_maps/plink.chr${chr}.GRCh38.map"
    output_file="linkage_maps/beagle/beagle_chr${chr}_GRCh38.map"
    
    if [ -f "$input_file" ]; then
        sed 's/^/chr/' "$input_file" > "$output_file"
        echo "Processed BEAGLE map for chromosome ${chr}"
    else
        echo "Warning: Input file $input_file not found"
    fi
done

# Process maps for GLIMPSE
echo "Processing maps for GLIMPSE..."
for chr in {1..22}; do
    input_file="linkage_maps/plink.chr${chr}.GRCh38.map"
    output_file="linkage_maps/glimpse/glimpse_chr${chr}_GRCh38.map"
    
    if [ -f "$input_file" ]; then
        # Create header
        echo "pos	chr	cM" > "$output_file"
        # Process the file and add 'chr' prefix to the chromosome field
        #awk '{print $4 "\t" "chr" $1 "\t" $3}' "$input_file" >> "$output_file"
        awk '{print $4 "\t" $1 "\t" $3}' "$input_file" >> "$output_file"
        echo "Processed GLIMPSE map for chromosome ${chr}"
    else
        echo "Warning: Input file $input_file not found"
    fi
done

echo "Done processing linkage maps for both BEAGLE and GLIMPSE"

