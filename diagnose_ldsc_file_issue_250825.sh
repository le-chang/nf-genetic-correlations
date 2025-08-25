#!/bin/bash

echo "=== Diagnosing LDSC file staging issue ==="
echo ""

# Check the failing work directory
WORK_DIR="/lustre06/project/6054517/lchang24/github/nf-genetic-correlations_fork/nf-genetic-correlations/work/72/ecd5df90733b4014d24db17c2c420b"

echo "1. Contents of failing work directory:"
ls -la "$WORK_DIR"
echo ""

echo "2. Looking for .sumstats.gz files in work directory:"
find "$WORK_DIR" -name "*.sumstats.gz" -ls 2>/dev/null || echo "No .sumstats.gz files found in work directory"
echo ""

echo "3. Contents of results/munged/ (where files actually exist):"
ls -la results/munged/
echo ""

echo "4. Check the command that was executed:"
if [ -f "$WORK_DIR/.command.sh" ]; then
    echo "Command script contents:"
    cat "$WORK_DIR/.command.sh"
else
    echo ".command.sh not found in work directory"
fi
echo ""

echo "5. Check for any symlinks or file references:"
find "$WORK_DIR" -type l -ls 2>/dev/null || echo "No symlinks found"
echo ""

echo "6. Check Nextflow's staging of inputs:"
if [ -f "$WORK_DIR/.command.begin" ]; then
    echo "Command begin log:"
    cat "$WORK_DIR/.command.begin"
else
    echo ".command.begin not found"
fi
echo ""

echo "7. Check if there are any input file references:"
find "$WORK_DIR" -name "*.txt" -o -name "*.tsv" -o -name "*.gz" | head -10
