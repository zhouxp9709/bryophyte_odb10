#!/bin/bash

# Check if a file path is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <path_to_HOG.tsv>"
    exit 1
fi

file_path=$1

# Process the file using awk
awk '
BEGIN {
    FS = OFS = "\t"
    Ppatens_col = 0
    Sphfalx_col = 0
}

NR == 1 {
    # Find the column numbers for species Ppatens and Sphfalx
    for (i = 1; i <= NF; i++) {
        if ($i == "Ppatens") {
            Ppatens_col = i
        }
        if ($i == "Sphfalx") {
            Sphfalx_col = i
        }
    }
}

NR > 1 {
    valid = 1
    total_count = 0  # Used to store the total number of gene copies
    for (i = 2; i <= NF; i++) {  # Check the gene copy number for each species
        if (i == Ppatens_col || i == Sphfalx_col) {
            if ($i == 0 || $i > 4) {
                valid = 0
                break
            }
        } else {
            if ($i == 0 || $i > 2) {
                valid = 0
                break
            }
        }
        total_count += $i
    }
    if (valid == 1 && total_count <= 33)
        print $1  # Print the HOG that meets the criteria
}
' $file_path