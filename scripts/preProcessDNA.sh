#!/bin/bash

# Usage: ./preprocessing_script.sh input_forward input_reverse output_forward output_reverse adapters

# Check if all required arguments are provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 input_forward input_reverse output_forward output_reverse adapters"
    exit 1
fi

# Input arguments
input_forward="$1"
input_reverse="$2"
output_forward="$3"
output_reverse="$4"
adapters="$5"

# Run bbduk for quality trimming and adapter removal
bbduk.sh in1="$input_forward" in2="$input_reverse" out1="$output_forward" out2="$output_reverse" \
  qtrim=rl trimq=20 adapters="$adapters"

# Generate quality control statistics
stats_file="preprocessing_stats.txt"
bbduk.sh in1="$input_forward" in2="$input_reverse" out1="$output_forward" out2="$output_reverse" \
  qtrim=rl trimq=20 adapters="$adapters" stats="$stats_file"

# Display summary statistics
cat "$stats_file"

# Run FastQC on the preprocessed data
fastqc "$output_forward" "$output_reverse"

# Optional: Remove intermediate files
#rm "$input_forward" "$input_reverse"
