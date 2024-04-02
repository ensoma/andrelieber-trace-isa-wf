#!/usr/bin/env python3

import pysam
import argparse
import sys
from collections import Counter

def main():
    # Main parser.
    parser = argparse.ArgumentParser()

    # Sub parsers.
    subparsers = parser.add_subparsers(help='sub-command help', dest='subcommand')

    # Softclip filtering subcommand.
    parser_softclip = subparsers.add_parser('softclip_filter', help='Filter reads based on softclipping')
    parser_softclip.add_argument('-i', '--input', help='BAM file to parse')
    parser_softclip.add_argument('-o', '--output', help='SAM file to write', default="-")
    parser_softclip.add_argument('-s', '--softclip', help='Softclip threshold', default=5)

    # Insertion site counting subcommand.
    parser_is = subparsers.add_parser('is_counter', help='Count insertion sites')
    parser_is.add_argument('-i', '--input', help='BAM file to parse')
    parser_is.add_argument('-o', '--output', help='BED file to write', default=sys.stdout, type=argparse.FileType('w'))
    parser_is.add_argument('-g', '--group_range', help='Insertion sites within this range (bp) will be grouped', default=5, type=int)

    # Parse the arguments.
    args = parser.parse_args()

    # Call the appropriate function.
    if args.subcommand == 'softclip_filter':
        softclip_filter(args.input, args.output, args.softclip)
    elif args.subcommand == 'is_counter':
        is_counter(args.input, args.output, args.group_range)

# Iterate through the BAM file.
# Discard reads where the R1 5' softclip is greater than the threshold.
# pysam cigartuple returns a list of tuples where (operation, length).
# 4 is the code for softclipping.
def softclip_filter(input_bam, output_sam, softclip_threshold):
    # Open the input and output files.
    infile = pysam.AlignmentFile(input_bam, 'rb')
    outfile = pysam.AlignmentFile(output_sam, 'w', template=infile)

    # Iterate through the BAM file.
    for line in infile:
        if line.is_read2:
            # Write the read2 to the output file.
            outfile.write(line)
        if line.is_read1 and line.is_reverse:
            # Retrieve the tuple at the end of the list (end of the cigar string).
            cigartuple = line.cigartuples[-1]
            # Write the read1 to the output file if the softclip is less than the threshold.
            if cigartuple[0] != 4 or cigartuple[1] < int(softclip_threshold):
                outfile.write(line)
        if line.is_read1 and not line.is_reverse:
            # Retrieve the tuple at the beginning of the list (beginning of the cigar string).
            cigartuple = line.cigartuples[0]
            # Write the read1 to the output file if the softclip is less than the threshold.
            if cigartuple[0] != 4 or cigartuple[1] < int(softclip_threshold):
                outfile.write(line)

# Iterate through the BAM file and retrieve the ISs.
# Soft clipped bases are excluded in the returned ranges.
def is_counter(input_bam, output_bed, group_range):
    # Open the BAM file.
    bam_data = pysam.AlignmentFile(input_bam, 'rb')

    # Iterate through the BAM file and retrieve the insertion sites.
    is_counter = Counter()
    for line in bam_data:
        if line.is_read1:
            if line.is_reverse:
                start = line.reference_end - 1
                end = line.reference_end
                strand = '-'
            else:
                start = line.reference_start
                end = line.reference_start + 1
                strand = '+'
            seq = line.reference_name
            is_counter.update({(seq, start, end, strand)})
    bam_data.close()

    is_ranges = [
        (seq, start, end, strand, count) for (seq, start, end, strand), count in is_counter.items()
    ]

    # For ranges within 5 bp of each other, keep the one with the highest count.
    # If there are multiple ranges with the same count, keep the first one.
    # The score of this range should be the sum of the counts of the ranges that were grouped.
    def group_ranges(ranges):
        # Sort the list by seqname, strand, and then start position
        # This makes it so that ranges withing 5 bp of each other are next to each other.
        sorted_ranges = sorted(ranges, key=lambda x: (x[0], x[3], x[1]))
        grouped_ranges = []  # This will hold the result
        current_group = []  # This holds the tuples of the current group
        for r in sorted_ranges:
            if not current_group:
                # If the current group is empty, this is the start of a new group
                current_group.append(r)
            else:
                # Check if the current tuple should be grouped with the current group.
                # The seqname and strand must be the same.
                # The start position must be within 5 bp of the end position of the last tuple in the current group.
                last_tuple = current_group[-1]
                if r[0] == last_tuple[0] and r[3] == last_tuple[3] and r[1] <= last_tuple[2] + group_range:
                    current_group.append(r)
                else:
                    # If the next item is not withing 5 bp of the last item in the current group,
                    # The IS location will be one with the highest score.
                    # If the scores are the same, the first one will be kept.
                    # The final score will be the sum of the scores of the ranges that were grouped. 
                    max_score_tuple = max(current_group, key=lambda x: x[4])
                    total_score = sum(x[4] for x in current_group)
                    grouped_ranges.append((max_score_tuple[0], max_score_tuple[1], max_score_tuple[2], max_score_tuple[3], total_score))
                    current_group = [r]
        # Need to process the last group.
        if current_group:
            max_score_tuple = max(current_group, key=lambda x: x[4])
            total_score = sum(x[4] for x in current_group)
            grouped_ranges.append((max_score_tuple[0], max_score_tuple[1], max_score_tuple[2], max_score_tuple[3], total_score))
        return grouped_ranges

    grouped_ranges = group_ranges(is_ranges)

    # Format the output as a bed file.
    bed = []
    for is_range in grouped_ranges:
        seqnames = is_range[0]
        start = is_range[1]
        end = is_range[2]
        score = is_range[4]
        strand = is_range[3]
        bed.append((seqnames, start, end, '.', score, strand))

    # Output the bed file.
    for line in bed:
        output_bed.write('\t'.join([str(x) for x in line]) + '\n')
    
    # If the input was a file, close it.
    if output_bed != sys.stdout:
        output_bed.close()

if __name__ == '__main__':
    main()