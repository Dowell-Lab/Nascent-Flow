#!/bin/sh

# Extract read statistics from fastqc output, providing a stopping point if we
# don't have enough.

# Set options
set -e
set -o nounset
set -o errexit
set -o pipefail

SRR=C413_1_S3_R1_001
GC=$(unzip -c "$(find . -name "$SRR"_fastqc.zip)" "$SRR"_fastqc/fastqc_data.txt \
         | grep "%GC" | grep -o "[0-9]*")
SEQ=$(unzip -c "$(find . -name "$SRR"_fastqc.zip)" "$SRR"_fastqc/fastqc_data.txt | \
          grep "Total Sequences" | \
          grep -o "[0-9]*")
DEDUP=$(unzip -c "$(find . -name "$SRR"_fastqc.zip)" "$SRR"_fastqc/fastqc_data.txt | \
            grep "#Total Deduplicated Percentage" | \
            grep -o "[0-9,.]*")

echo -e "SRR\t%GC\tTotal_Sequences\t%Total_Deduplicated"
echo -e "$SRR"\t"$GC"\t"$SEQ"\t"$DEDUP"

if [ "$GC" -lt 3000000 ]
then
		echo "Insufficient read depth"
		exit 1
fi
