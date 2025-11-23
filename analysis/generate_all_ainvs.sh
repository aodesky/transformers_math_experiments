#!/bin/bash
# Generate ainvs files for all search_output_*.txt files in elliptic_gpu_test

# Find the GPU test directory
GPU_DIR="output/elliptic_big_run/a6jmms42lu"

if [ ! -d "../$GPU_DIR" ]; then
    echo "Error: Directory ../$GPU_DIR not found"
    exit 1
fi

echo "Finding all search_output_*.txt files in ../$GPU_DIR"
echo

# Find all search_output_N.txt files (excluding tokenized versions)
search_files=$(find "../$GPU_DIR" -name "search_output_*.txt" -type f | grep -v tokenized | sort -V)

if [ -z "$search_files" ]; then
    echo "No search_output_*.txt files found"
    exit 1
fi

# Count files
num_files=$(echo "$search_files" | wc -l)
echo "Found $num_files files to process"
echo

# Process each file
count=0
for file in $search_files; do
    count=$((count + 1))
    echo "[$count/$num_files] Processing $file"
    sage -python generate_ainvs_dataset.py "$file"
    echo
done

echo "All done!"
echo
echo "Generated files:"
find "../$GPU_DIR" -name "ainvs_of_search_output_*.txt" | sort -V
