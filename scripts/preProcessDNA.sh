#!/bin/bash

# Set input and output file paths
input_forward="input_forward.fastq"
input_reverse="input_reverse.fastq"
output_forward="cleaned_forward.fastq"
output_reverse="cleaned_reverse.fastq"

# Adapter sequences file
adapters="adapters.fa"

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
