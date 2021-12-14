#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Usage: chomper.v0.1.py -H <STR> -b <STR> -s <STR> -q <FLT> -o <STR> [-h]

  [Options]
    -H, --haplotypes <STR>                      Haplotypes file from HAPCUT2
    -b, --bam <STR>                             Sorted and indexed bam file, containing primary alignments only
    -s, --sequence <STR>                        Sequence/chromosome to analyse
    -q, --quality <FLT>                         Minimum hapcut2 mismatch score for phased variant
    -o, --outprefix <STR>                       Output prefix
    -h, --help                                  Show this message

"""

import sys
from docopt import docopt
import pysam
import collections

# get the base of each haplotype at a SNP
def get_base_of_haplotype(haplotype, ref, alt):
	if haplotype == 0:
		return ref
	else:
		return alt.split(",")[haplotype-1] # this allow tri-allelic calls

def read_hapcut(haplotypes_f):
	# collect the phased variants
	print("[+] Reading in phased variants...")
	phased_variants = []
	with open(haplotypes_f, "r") as haplotypes:
		for line in haplotypes:
			fields = line.rstrip().split()
			quality = float(fields[10])
			# only interested in SNPs with quality above some value
			if quality >= float(args['--quality']):
				# record information about this phased SNP
				pos = int(fields[4]) # position
				ref = str(fields[5]) # ref allele
				alt = str(fields[6]) # alt allele
				haplotype_1 = int(fields[1]) # haplotype_1_allele, ref or alt
				haplotype_2 = int(fields[2]) # haplotype_2_allele, ref or alt
				haplotype_1_base = get_base_of_haplotype(haplotype_1, ref, alt) # base of haplotype 1
				haplotype_2_base = get_base_of_haplotype(haplotype_2, ref, alt) # base of haplotype 2
				phased_variants.append([pos, haplotype_1_base, haplotype_2_base])
	return phased_variants

def bam_parser(phased_variants, bam_f, sequence):
	# open the bam
	print("[+] Reading in bam...")
	with pysam.AlignmentFile(bam_f, "rb") as samfile:
		print("[+] Checking reads against phased variants...")
		# store phase information in a set
		read_dict = collections.defaultdict(set) 
		for position, haplotype_1_base, haplotype_2_base in phased_variants:
			# get reads that overlap the phased SNP
			for column in samfile.pileup(sequence, position-1, position, ignore_overlaps=False, truncate=True):
				#print(position)
				for pileupread in column.pileups:
					# don't process read if del at the phased variant
					if pileupread.query_position:
						# need to ask if mate is on --sequence too, if not we are not interested
						if pileupread.alignment.next_reference_name == sequence:
							read_name = pileupread.alignment.query_name
							read_base = pileupread.alignment.query_sequence[pileupread.query_position]
							# record phase information of the read
							if read_base == haplotype_1_base:
								read_dict[read_name].add(1)
							elif read_base == haplotype_2_base:
								read_dict[read_name].add(2)
	return read_dict


def parition_reads(read_names, bam_f):
	haplotype_1_read_1 = [] # lists to store reads
	haplotype_1_read_2 = []
	haplotype_2_read_1 = []
	haplotype_2_read_2 = []
	haplotype_1_count = 0 # counters
	haplotype_2_count = 0
	conflict_count = 0
	print("[+] Indexing bam by read name...")
	bam_file = pysam.AlignmentFile(bam_f, "rb")
	indexed_bam_file = pysam.IndexedReads(bam_file)
	indexed_bam_file.build() # build index by read name
	print("[+] Finding reads by name...")
	for read_name in read_names: # go through reads with phased variant info
		for read in indexed_bam_file.find(read_name):
			nucleotides = read.get_forward_sequence() # sequence
			qualities = "".join([chr(qual + 33) for qual in read.get_forward_qualities()]) # qualities
			if read.is_read1:
				if read_names[read_name] == {1}:
					haplotype_1_read_1.append("@" + read_name)
					haplotype_1_read_1.append(nucleotides)
					haplotype_1_read_1.append("+")
					haplotype_1_read_1.append(qualities)
					haplotype_1_count += 1
				elif read_names[read_name] == {2}:
					haplotype_2_read_1.append("@" + read_name)
					haplotype_2_read_1.append(nucleotides)
					haplotype_2_read_1.append("+")
					haplotype_2_read_1.append(qualities)
					haplotype_2_count += 1
				else:
					conflict_count += 1
			elif read.is_read2:
				if read_names[read_name] == {1}:
					haplotype_1_read_2.append("@" + read_name)
					haplotype_1_read_2.append(nucleotides)
					haplotype_1_read_2.append("+")
					haplotype_1_read_2.append(qualities)
				elif read_names[read_name] == {2}:
					haplotype_2_read_2.append("@" + read_name)
					haplotype_2_read_2.append(nucleotides)
					haplotype_2_read_2.append("+")
					haplotype_2_read_2.append(qualities)
	print("[=] haplotype_1 count:", haplotype_1_count)
	print("[=] haplotype_2 count:", haplotype_2_count)
	print("[=] conflict count:", conflict_count)
	return [haplotype_1_read_1, haplotype_1_read_2, haplotype_2_read_1, haplotype_2_read_2]

def write_fastq(partioned_reads, outprefix):
	print("[+] Writing fastq files...")
	with open(str(outprefix) + ".haplotype_1.R1.fastq", "w") as h1r1:
		h1r1.write("\n".join(partioned_reads[0]))
	with open(str(outprefix) + ".haplotype_1.R2.fastq", "w") as h1r2:
		h1r2.write("\n".join(partioned_reads[1]))
	with open(str(outprefix) + ".haplotype_2.R1.fastq", "w") as h2r1:
		h2r1.write("\n".join(partioned_reads[2]))
	with open(str(outprefix) + ".haplotype_2.R2.fastq", "w") as h2r2:
		h2r2.write("\n".join(partioned_reads[3]))


if __name__ == '__main__':
	__version__ = '0.1'
	args = docopt(__doc__)
	print(args)
	phased_variants = read_hapcut(args['--haplotypes'])
	read_names = bam_parser(phased_variants, args['--bam'], args['--sequence'])
	partioned_reads = parition_reads(read_names, args['--bam'])
	write_fastq(partioned_reads, args['--outprefix'])
	print("[=] Done chomping")