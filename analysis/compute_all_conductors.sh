#!/bin/bash
# Compute conductors for all ainvs_of_search_output_*.txt files

# Find the GPU run directory
GPU_DIR="output/elliptic_big_run/a6jmms42lu"

if [ ! -d "../$GPU_DIR" ]; then
    echo "Error: Directory ../$GPU_DIR not found"
    exit 1
fi

echo "Finding all ainvs_of_search_output_*.txt files in ../$GPU_DIR"
echo

# Find all ainvs_of_search_output_N.txt files
ainvs_files=$(find "../$GPU_DIR" -name "ainvs_of_search_output_*.txt" -type f | sort -V)

if [ -z "$ainvs_files" ]; then
    echo "No ainvs_of_search_output_*.txt files found"
    exit 1
fi

# Count files
num_files=$(echo "$ainvs_files" | wc -l)
echo "Found $num_files files to process"
echo

# Process each file
count=0
for file in $ainvs_files; do
    count=$((count + 1))
    echo "[$count/$num_files] Processing $file"
    sage -python compute_conductors_parallel.sage.py "$file"
    echo
done

echo "All done!"
echo
echo "Generated files:"
find "../$GPU_DIR" -name "conductors_of_search_output_*.txt" | sort -V
