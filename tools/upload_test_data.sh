#!/usr/bin/env bash
#
# Package Fortran reference data into tarballs and upload everything to S3.
#
# Usage: ./tools/upload_test_data.sh
#
# Requires s3cmd configured with ~/.s3cfg-gadopt for the DigitalOcean Spaces
# gadopt bucket (syd1 region).

set -euo pipefail

GEODYN="/Users/sghelichkhani/Workplace/tomographic_filtering/srts/geodyn"
DATA_DIR="$(cd "$(dirname "$0")/.." && pwd)/src/srts/data"
S3_BUCKET="s3://gadopt/srts"
S3CMD="s3cmd -c $HOME/.s3cfg-gadopt"
WORKDIR=$(mktemp -d)

trap 'rm -rf "$WORKDIR"' EXIT

echo "=== Building tier_1.tar.gz ==="
TIER1="$WORKDIR/tier1"
mkdir -p "$TIER1/rawfiles" "$TIER1/pwrfiles" "$TIER1/comparefiles" "$TIER1/outputfiles"

# .sph files
cp "$GEODYN/inpm.S40.examplefile.dvs.repar.sph" "$TIER1/"
cp "$GEODYN/oupm.S40.examplefile.dvs.filt.sph" "$TIER1/"

# .raw files at 7 depths, both repar and filt
for depth in 0100 0500 0675 1000 1500 2000 2500; do
    cp "$GEODYN/rawfiles/inpm.S40.examplefile.dvs.repar.${depth}.raw" "$TIER1/rawfiles/"
    cp "$GEODYN/rawfiles/oupm.S40.examplefile.dvs.filt.${depth}.raw" "$TIER1/rawfiles/"
done

# Power spectrum files
cp "$GEODYN/pwrfiles/inpm.S40.examplefile.dvs.repar.pwr.dat" "$TIER1/pwrfiles/"
cp "$GEODYN/pwrfiles/inpm.S40.examplefile.dvs.repar.pwr.deg.dat" "$TIER1/pwrfiles/"
cp "$GEODYN/pwrfiles/oupm.S40.examplefile.dvs.filt.pwr.dat" "$TIER1/pwrfiles/"
cp "$GEODYN/pwrfiles/oupm.S40.examplefile.dvs.filt.pwr.deg.dat" "$TIER1/pwrfiles/"

# Correlation files
cp "$GEODYN/comparefiles/corr.S40RTS..inpm.S40.examplefile.dvs.repar.corr.dat" "$TIER1/comparefiles/"
cp "$GEODYN/comparefiles/corr.S40RTS..inpm.S40.examplefile.dvs.repar.corr.deg.dat" "$TIER1/comparefiles/"
cp "$GEODYN/comparefiles/corr.S40RTS..oupm.S40.examplefile.dvs.filt.corr.dat" "$TIER1/comparefiles/"
cp "$GEODYN/comparefiles/corr.S40RTS..oupm.S40.examplefile.dvs.filt.corr.deg.dat" "$TIER1/comparefiles/"

# Point evaluation output files
cp "$GEODYN/outputfiles/inpm.S40.examplefile.dvs.repar.0677.dat" "$TIER1/outputfiles/"
cp "$GEODYN/outputfiles/inpm.S40.examplefile.dvs.repar.1038.dat" "$TIER1/outputfiles/"
cp "$GEODYN/outputfiles/inpm.S40.examplefile.dvs.repar.1806.dat" "$TIER1/outputfiles/"
cp "$GEODYN/outputfiles/oupm.S40.examplefile.dvs.filt.0677.dat" "$TIER1/outputfiles/"
cp "$GEODYN/outputfiles/oupm.S40.examplefile.dvs.filt.1038.dat" "$TIER1/outputfiles/"

tar -czf "$WORKDIR/tier_1.tar.gz" -C "$TIER1" .
echo "tier_1.tar.gz: $(du -h "$WORKDIR/tier_1.tar.gz" | cut -f1)"


echo "=== Building tier_2.tar.gz ==="
TIER2="$WORKDIR/tier2"
mkdir -p "$TIER2/examplemodel"
cp "$GEODYN/examplemodel/"* "$TIER2/examplemodel/"

tar -czf "$WORKDIR/tier_2.tar.gz" -C "$TIER2" .
echo "tier_2.tar.gz: $(du -h "$WORKDIR/tier_2.tar.gz" | cut -f1)"


echo "=== Uploading to $S3_BUCKET ==="

# Tarballs
$S3CMD put --acl-public "$WORKDIR/tier_1.tar.gz" "$S3_BUCKET/tier_1.tar.gz"
$S3CMD put --acl-public "$WORKDIR/tier_2.tar.gz" "$S3_BUCKET/tier_2.tar.gz"

# HDF5 model data files
for f in S12RTS.h5 S20RTS.h5 S40RTS.h5; do
    if [ -f "$DATA_DIR/$f" ]; then
        echo "Uploading $f..."
        $S3CMD put --acl-public "$DATA_DIR/$f" "$S3_BUCKET/$f"
    else
        echo "WARNING: $DATA_DIR/$f not found, skipping"
    fi
done

echo "=== Done ==="
echo "Verify with: $S3CMD ls $S3_BUCKET/"
