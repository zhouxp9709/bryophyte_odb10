
#!/bin/bash

# Set environment and variables
set -e
PEP_DIR="/vol2/BUSCO/bryophyte/HOG/02.copy.pep"
SPLIT_DIR="/vol2/BUSCO/bryophyte/HOG/03.len.filter.work"
OUTPUT_DIR="/vol2/BUSCO/bryophyte/HOG/04.len.pep"
FILES=($(ls ${PEP_DIR}/*.fa))
TOTAL_FILES=${#FILES[@]}
FILES_PER_DIR=$((TOTAL_FILES / 2))
COUNTER=0
DIR_INDEX=1

# Create target directories and split files
for FILE in "${FILES[@]}"; do
    DIR_PATH="${SPLIT_DIR}/dir_${DIR_INDEX}"
    mkdir -p "${DIR_PATH}"
    
    FILE_BASENAME=$(basename "$FILE")
    cp "$FILE" "${DIR_PATH}/"

    let COUNTER=COUNTER+1
    if [ $((COUNTER % FILES_PER_DIR)) -eq 0 ] && [ $DIR_INDEX -lt 2 ]; then
        let DIR_INDEX=DIR_INDEX+1
    fi
done

# Create work.sh script in each directory
for DIR in ${SPLIT_DIR}/dir_*; do
    cat > ${DIR}/work.sh <<EOF
#!/bin/bash
export PATH=/public/home/miniconda3/envs/mamba/bin/:\$PATH

cd ${DIR}

# Use samtools faidx to compute sequence lengths
filter_sequences() {
    local file="\$1"
    local lengths_file="\${file}_lengths.txt"
    local output_file="${OUTPUT_DIR}/\$(basename "\$file")"
    
    # Generate index using samtools faidx
    samtools faidx "\$file"
    
    # Extract sequence length information
    awk '{print \$2}' "\${file}.fai" > "\$lengths_file"

    # Read length information
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
    rsd_threshold=20  # Conservative setting, allows 20% RSD

    # Check if RSD exceeds threshold
    if (( \$(echo "\$rsd > \$rsd_threshold" | bc -l) )); then
        echo "File \$file has high relative sequence length standard deviation: \$rsd%"
        return
    fi

    # Copy files that meet criteria to the output directory
    cp "\$file" "\$output_file"
}

export -f filter_sequences

# Process each file
for file in *.fa; do
    filter_sequences "\$file"
done
EOF
    chmod +x ${DIR}/work.sh
done

# Submit jobs
for DIR in ${SPLIT_DIR}/dir_*; do
    qsub -cwd -pe smp 1 ${DIR}/work.sh
done
