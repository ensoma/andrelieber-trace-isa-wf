#!/usr/bin/env python

import sys
from Bio import SeqIO

# Check if the command-line argument is provided
if len(sys.argv) != 2:
    print("Usage: python script.py input.genbank")
    sys.exit(1)

# Get the input GenBank file name from the command-line argument
input_genbank = sys.argv[1]

# Open the GenBank file and parse it
with open(input_genbank, 'r') as handle:
    for record in SeqIO.parse(handle, 'genbank'):
        for feature in record.features:
            # Extract the start, end, and strand information
            start = feature.location.start
            end = feature.location.end
            strand = '*'
            
            # Attempt to extract the feature label, use 'unknown' if not available
            name = feature.qualifiers.get('label', ['unknown'])[0]
            
            # Print the BED format: chrom start end name score strand
            print(f"{record.id}\t{start}\t{end}\t{name}\t0\t{strand}")