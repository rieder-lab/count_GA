# may need to install some of these packages using anaconda or pip

import os # standard package for interfacing with operating system
import sys # same
import gzip # package for reading/writign gzip files
import regex # package for doing
import logging # package for logging runtime info
import csv # package for writign csv/tsv files
import argparse # package for passing command line arguments to python scripts
from Bio import SeqIO # package for easily parsing through fastq files. You may need to install biopython before you can load this package. See: https://biopython.org/wiki/Download

### This is a function that parses command line arguements. This is not necessary. You could simply hard code these variable, but doing it this way allows you to change the varibales more easily.
def process_args():
    parser = argparse.ArgumentParser(description='Finds H3H4p reads and counts GA repeats.')
    parser.add_argument('-i', '--input_fastq', required=True, type=str, help='Input fastq file - gzipped.')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output path and file prefix.')
    parser.add_argument('-l', '--left_anchor', required=False, type=str, default="TAGCAATCGT", help='Left (5prime) anchor sequence to search for. [default=TAGCAATCGT]')
    parser.add_argument('-r', '--right_anchor', required=False, type=str, default="CATTTCATTTGACGAGC", help='Right (3prime) anchor sequence to search for. [default=CATTTCATTTGACGAGC]')
    parser.add_argument('-s', '--repeat_seed', required=False, type=str, default="AGAGAG", help='Seed sequence to identify GA repeat stretch. [default=AGAGAG]')
    parser.add_argument('-m', '--mismatches', required=False, type=int, default=0, help='How many mismatches to allow in searchign for anchor sequences? [default=0]')
    return parser.parse_args()

### This is the main code of the script

# parse command line arguements
args = process_args()

#configure the logging output
logging.basicConfig(format='[%(filename)s] %(asctime)s %(levelname)s: %(message)s', datefmt='%I:%M:%S', level=logging.INFO)

# log run info
logging.info("### Preparing to count GA repeats in " + args.input_fastq +
             "\n### Using left anchor " + args.left_anchor + ", right anchor " + args.right_anchor +
             ",\n### GA repeat seed " + args.repeat_seed + ", and allowing " + str(args.mismatches) + " mismatches...")

# format query strings to include mismatches if desired
left_query = args.left_anchor + "{s<=" + str(args.mismatches) + "}"
right_query = args.right_anchor + "{s<=" + str(args.mismatches) + "}"

# define output files
left_only_out = gzip.open(args.output + "_left_only.fastq.gz", "wt")
right_only_out = gzip.open(args.output + "_right_only.fastq.gz", "wt")
dual_match_out = gzip.open(args.output + "_dual_match.fastq.gz", "wt")
no_match_out = gzip.open(args.output + "_no_match.fastq.gz", "wt")
strange_match_out = gzip.open(args.output + "_strange_match.fastq.gz", "wt")

# set counters
total_count = 0
left_match_count = 0
left_match_rc_count = 0
right_match_count = 0
right_match_rc_count = 0
dual_match_count = 0
dual_match_rc_count = 0
no_match_count = 0
strange_match_count = 0

# iterate through fastq file read by read to find anchor match(es)
logging.info("### Iterating through fastq file...")
with gzip.open(args.input_fastq, "rt") as handle:
    for read in SeqIO.parse(handle, "fastq"):

        total_count += 1

        # set match trackers to False
        left_match = False
        left_match_rc = False
        right_match = False
        right_match_rc = False

        ### Search for matches to anchor sequences
        # search sequence for left anchor
        if regex.search(left_query, str(read.seq)):
            left_match = True
        # search reverse complmenet of sequence for left anchor
        if regex.search(left_query, str(read.seq.reverse_complement())):
            left_match_rc = True
        # search sequence for right anchor
        if regex.search(right_query, str(read.seq)):
            right_match = True
        # search reverse complmenet of sequence for right anchor
        if regex.search(right_query, str(read.seq.reverse_complement())):
            right_match_rc = True

        ### classify sequences according to matches and write to seperate files
        # no matches
        if(not(left_match or left_match_rc or right_match or right_match_rc)):
            no_match_count += 1
            no_match_out.write(read.format("fastq"))
            continue
        # unexpected configurations like matching both forward and reverse complement
        if((left_match and left_match_rc) or (right_match and right_match_rc) or (left_match and right_match_rc) or (right_match and left_match_rc)):
            strange_match_count += 1
            strange_match_out.write(read.format("fastq"))
            continue
        # match both left and right anchors in same orientation (checking to make sure they are in the correct order)
        if left_match and right_match:
            if regex.search(left_query, str(read.seq)).span()[0] < regex.search(right_query, str(read.seq)).span()[0]:
                dual_match_count += 1
                dual_match_out.write(read.format("fastq"))
            else:
                strange_match_count += 1
                strange_match_out.write(read.format("fastq"))
            continue
        elif left_match_rc and right_match_rc:
            if regex.search(left_query, str(read.seq.reverse_complement())).span()[0] < regex.search(right_query, str(read.seq.reverse_complement())).span()[0]:
                dual_match_rc_count += 1
                dual_match_out.write(read.format("fastq"))
            else:
                strange_match_count += 1
                strange_match_out.write(read.format("fastq"))
            continue
        # match left anchor only
        if left_match or left_match_rc:
            if left_match:
                left_match_count += 1
            if left_match_rc:
                left_match_rc_count += 1
            left_only_out.write(read.format("fastq"))
            continue
        # match right anchor only
        if right_match or right_match_rc:
            if right_match:
                right_match_count += 1
            if right_match_rc:
                right_match_rc_count += 1
            right_only_out.write(read.format("fastq"))
            continue

# close file handles
handle.close()
left_only_out.close()
right_only_out.close()
dual_match_out.close()
no_match_out.close()
strange_match_out.close()

# log counts
logging.info("### Total seqs examined:\t" + str(total_count))
logging.info("### No anchor matches:\t" + str(no_match_count))
logging.info("### Strange anchor match configurations:\t" + str(strange_match_count))
logging.info("### Left only anchor matches:\t" + str(left_match_count))
logging.info("### Left only anchor matches(rc):\t" + str(left_match_rc_count))
logging.info("### Right only anchor matches:\t" + str(right_match_count))
logging.info("### Right only anchor matches(rc):\t" + str(right_match_rc_count))
logging.info("### Dual anchor matches:\t" + str(dual_match_count))
logging.info("### Dual anchor matches(rc):\t" + str(dual_match_rc_count))

### This block parses through the dual match fastq file created above to find matches to the ga seed and then count repeats

# set counters
total_count = 0
ga_found_count = 0
no_ga_count = 0

#open tsv output file to record repeat lengths
out_file = open(args.output + "_ga_counts.tsv", 'wt')
tsv_out = csv.writer(out_file, delimiter='\t')
tsv_out.writerow(['GA_rep_len', 'GA_seq', 'full_seq', 'seq_name'])

logging.info("### Searching dual anchor match seqs for GA repeats...")

# parse through dual match fastq to count GA repeats
with gzip.open(args.output + "_dual_match.fastq.gz", "rt") as handle:
    for read in SeqIO.parse(handle, "fastq"):

        total_count += 1

        # intitialize vairable to capture to repeat length and repeat substring
        temp_ga_str = ""
        temp_ga_len = 0

        #determine orientation and standardize
        if regex.search(left_query, str(read.seq)):
            temp_seq = str(read.seq)
        else:
            temp_seq = str(read.seq.reverse_complement())

        # cut seqence at left anchor
        temp_seq_trim = temp_seq[regex.search(left_query, temp_seq).span()[1]:]

        # cut seqence at right anchor
        temp_seq_trim = temp_seq_trim[:regex.search(right_query, temp_seq_trim).span()[0]]

        # find GA seed and loop through sequence to count GAs
        temp_match_seed = regex.search(args.repeat_seed, temp_seq_trim)
        if temp_match_seed:
            temp_ga_start = temp_match_seed.span()[0]
            temp_seq_trim_ga = temp_seq_trim[temp_ga_start+6:]
            ga_found_count += 1
            temp_ga_len += 6
            ga_end = False
            while not ga_end:
                # note that regex.match() requires the match to be at the start of the string, whereas regex.search() will find a match anywhere in the string
                temp_match_ag = regex.match('AG', temp_seq_trim_ga)
                if not temp_match_ag:
                    ga_end = True
                else:
                    temp_ga_len += 2
                    temp_seq_trim_ga = temp_seq_trim_ga[2:]
            if(regex.match('A', temp_seq_trim_ga)):
                temp_ga_len += 1
            temp_ga_str = temp_seq_trim[temp_ga_start:temp_ga_start+temp_ga_len]
        else:
            no_ga_count += 1
            temp_ga_len = "NA"
            temp_ga_str = "NA"

        # write info to tsv
        tsv_out.writerow([temp_ga_len, temp_ga_str, temp_seq, read.name])

handle.close()
out_file.close()

logging.info("### Total seqs examined:\t" + str(total_count))
logging.info("### GA repeat seed found:\t" + str(ga_found_count))
logging.info("### GA repeat seed not found:\t" + str(no_ga_count))
