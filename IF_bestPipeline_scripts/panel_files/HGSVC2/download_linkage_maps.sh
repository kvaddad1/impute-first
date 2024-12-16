set -x  
mkdir -p linkage_maps
wget "https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip" -O linkage_maps/plink.GRCh38.map.zip
unzip -o linkage_maps/plink.GRCh38.map.zip -d linkage_maps
rm linkage_maps/plink.GRCh38.map.zip
echo "done downloading linkage_maps"
