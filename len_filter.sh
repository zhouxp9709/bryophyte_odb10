#!/bin/bash

# Setup environment and variables
set -e
PEP_DIR="/vol2/BUSCO/bryophyte/HOG/02.copy.pep"
OUTPUT_DIR="/vol2/BUSCO/bryophyte/HOG/03.len.pep"
FILES=($(ls ${PEP_DIR}/*.fa))

# Load the necessary environment
export PATH=/public/home/miniconda3/envs/mamba/bin/:\$PATH

# Function to filter sequences based on relative standard deviation (RSD) of sequence lengths
filter_sequences() {
    local file="\$1"
    local lengths_file="\${file}_lengths.txt"
    local output_file="${OUTPUT_DIR}/\$(basename "\$file")"
    
    # Generate index using samtools faidx
    samtools faidx "\$file"
    
    # Extract sequence length information
    awk '{print \$2}' "\${file}.fai" > "\$lengths_file"

    # Read lengths into an array
    lengths=(\$(cat "\$lengths_file"))

    # Calculate mean and standard deviation
    sum=0
    sumsq=0
    count=\${#lengths[@]}
    for length in "\${lengths[@]}"; do
        sum=\$((sum + length))
        sumsq=\$((sumsq + length * length))
    done
    mean=\$(echo "scale=2; \$sum / \$count" | bc)
    stddev=\$(echo "scale=2; sqrt(\$sumsq / \$count - \$mean * \$mean)" | bc)

    # Calculate relative standard deviation (RSD)
    if (( \$(echo "\$mean > 0" | bc -l) )); then
        rsd=\$(echo "scale=2; \$stddev / \$mean * 100" | bc)
    else
        rsd=0
    fi

    # Set RSD threshold
    rsd_threshold=20  # Conservative threshold, allowing 20% RSD

    # Check if RSD exceeds the threshold
    if (( \$(echo "\$rsd > \$rsd_threshold" | bc -l) )); then
        echo "File \$file has high relative sequence length standard deviation: \$rsd%"
        return
    fi

    # Copy file that meets the criteria to the output directory
    cp "\$file" "\$output_file"
}

# Process each file in the PEP_DIR
for file in "\${FILES[@]}"; do
    filter_sequences "\$file"
done
