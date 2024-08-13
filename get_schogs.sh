#!/bin/bash

# Setup environment and variables
set -e
PEP_DIR="/vol2/BUSCO/bryophyte/HOG/00.pep"
SPLIT_DIR="/vol2/BUSCO/bryophyte/HOG/01.copy.filter.work"
OUTPUT_DIR="/vol2/BUSCO/bryophyte/HOG/02.copy.pep"
FILES=($(ls ${PEP_DIR}/*.fa))
TOTAL_FILES=${#FILES[@]}
FILES_PER_DIR=$((TOTAL_FILES / 30))
COUNTER=0
DIR_INDEX=1

# Create target directories and split files
for FILE in "${FILES[@]}"; do
    DIR_PATH="${SPLIT_DIR}/dir_${DIR_INDEX}"
    mkdir -p "${DIR_PATH}"
    
    FILE_BASENAME=$(basename "$FILE")
    cp "$FILE" "${DIR_PATH}/"

    let COUNTER=COUNTER+1
    if [ $((COUNTER % FILES_PER_DIR)) -eq 0 ] && [ $DIR_INDEX -lt 30 ]; then
        let DIR_INDEX=DIR_INDEX+1
    fi
done

# Create work.sh script in each directory
for DIR in ${SPLIT_DIR}/dir_*; do
    cat > ${DIR}/work.sh <<EOF
#!/bin/bash
export PATH=/public/home/miniconda3/envs/mamba/bin/:\$PATH

cd ${DIR}
process_hog() {
    local file="\$1"
    local dir=\$(dirname "\$file")
    local base=\$(basename "\$file" .fa)
    local output_file="${OUTPUT_DIR}/\${base}_best.fa"
    local tmp_aligned="\${dir}/\${base}_aligned.fa"

    # Perform global alignment using MAFFT
    mafft --auto --thread 6 "\$file" > "\$tmp_aligned"

    # Calculate similarity matrix using ClustalW
    clustalw -INFILE="\$tmp_aligned" -OUTPUT=FASTA

    # Extract gene IDs for each species
    grep "^>" "\$tmp_aligned" | sed 's/^>//' | cut -d '|' -f 1 | sort | uniq | while read -r species; do
        # Calculate average similarity score for each gene and save to temp file
        local tmp_scores=\$(mktemp)
        grep "^>\$species" "\$tmp_aligned" | while read -r gene; do
            gene_id=\$(echo "\$gene" | sed 's/^>//')
            sim_score=\$(grep -A 1 "\$gene" "\$tmp_aligned" | tail -n 1 | awk '{ total += length; } END { print total/NR; }')
            echo "\$gene_id \$sim_score" >> "\$tmp_scores"
        done

        # Select the gene with the highest similarity score
        best_gene=\$(sort -k2,2nr "\$tmp_scores" | head -n 1 | awk '{print \$1}')
        echo "Best \$species gene for \$file: \$best_gene"

        # Extract the best gene sequence from the original file and append to output
        awk -v best_gene="\$best_gene" '
            BEGIN { output_best = 0 }
            /^>/ {
                if (\$0 == ">" best_gene) {
                    output_best = 1
                } else {
                    output_best = 0
                }
            }
            { if (output_best) print \$0 }
        ' "\$file" >> "\$output_file"

        # Clean up temp file
        rm "\$tmp_scores"
    done

    # Clean up temporary alignment file
    rm "\$tmp_aligned"
}

export -f process_hog

# Process each .fa file in the directory
for file in *.fa; do
    process_hog "\$file"
done
EOF
    chmod +x ${DIR}/work.sh
done

# Submit jobs
for DIR in ${SPLIT_DIR}/dir_*; do
    qsub -cwd -pe smp 6 ${DIR}/work.sh
done