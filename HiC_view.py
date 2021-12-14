#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Usage: HiC_view.py -g <STR> -b <INT> -m <STR> -s <INT> [-t <BOL> -f <INT> -h]

  [Options]
    -g, --genomefile <STR>                     Genomefile for the sequences you want to view. Tab delimited with sequence, length, orientation (+/-).
    -b, --bins <INT>                           Number of bins in the heatmap (determines resolution) [default: 1000]
    -m, --mergenodups <STR>                    mergenodups.txt file from Juicer
    -s, --saturation <INT>                     Number of reads in a bin to reach colour saturation (determines saturation) [default: 50]
    -t, --stagger <BOL>                        Stagger sequence names over two columns/rows. Use this when there are many sequences. Choose from True or False [default: False]
    -f, --fsize <INT>                          Front size for labels [default: 8]
    -h, --help                                 Show this message

"""

import sys
from docopt import docopt
import numpy as np
import math
from matplotlib import pyplot as plt

# read in genome file and record cumulative lengths and sign in a dictionary
def generate_genomefile_dict(genomefile):
	genomefile_dict = {}
	sign_dict = {}
	cumulative_genome = 0
	with open(genomefile, "r") as fin:
		for line in fin:
			line = line.rstrip()
			chromosome = line.split("\t")[0]
			chromosome_length = int(line.split("\t")[1])
			sign = line.split("\t")[2]
			genomefile_dict[chromosome] = [cumulative_genome, cumulative_genome + chromosome_length]
			sign_dict[chromosome] = sign
			cumulative_genome += chromosome_length
	return genomefile_dict, cumulative_genome, sign_dict

# convert coordinates to bin numbers for plotting
def get_bin_number(chromosome, coord, genomefile_dict, bin_size):
	transformed_coord = coord + genomefile_dict[chromosome][0]
	bin_number = int(math.floor(transformed_coord / bin_size))
	return bin_number


args = docopt(__doc__)

# just in case to numpy array gets too big
if int(args['--bins']) > 40000:
	sys.exit("[X] Bin number requires {} Gb of memory".format((int(args['--bins'])**2 * 8) / 1000000000))

genomefile_dict, cumulative_genome, sign_dict = generate_genomefile_dict(args['--genomefile'])
bin_size = int(math.ceil(cumulative_genome / int(args['--bins'])))
print("[=] Bin size set to: {}".format(bin_size))

# numpy array for recording contact counts
HiC_matrix = np.zeros((int(args['--bins']), int(args['--bins'])), dtype=np.int64)

# loop through mnd file record HiC contacts in different bins
with open(args['--mergenodups'], "r") as mnd:
	for line in mnd:
		line = line.rstrip()
		split_line = line.split(" ")
		chromosome_A = split_line[1]
		chromosome_B = split_line[5]
		coord_A = int(split_line[2])
		coord_B = int(split_line[6])
		if chromosome_A in genomefile_dict.keys() and chromosome_B in genomefile_dict.keys():
			if sign_dict[chromosome_A] == "-": # flip coords if sign is minus
				coord_A = (genomefile_dict[chromosome_A][1] - genomefile_dict[chromosome_A][0]) - coord_A + 1
			if sign_dict[chromosome_B] == "-":
				coord_B = (genomefile_dict[chromosome_B][1] - genomefile_dict[chromosome_B][0]) - coord_B + 1
			bin_A = get_bin_number(chromosome_A, coord_A, genomefile_dict, bin_size)
			bin_B = get_bin_number(chromosome_B, coord_B, genomefile_dict, bin_size)
			HiC_matrix[bin_A, bin_B] += 1
			HiC_matrix[bin_B, bin_A] += 1


# plot the numpy array counts using imshow
plt.imshow(HiC_matrix, cmap='Reds', vmin=0, vmax=args['--saturation'], interpolation='nearest')

# add sequence boundaries
for chromosome in genomefile_dict:
	if genomefile_dict[chromosome][0] == genomefile_dict[chromosome][1]:
		pass
	else:
		chrom_boundry = genomefile_dict[chromosome][0] / bin_size
		plt.axvline(x=chrom_boundry, color='black', lw=0.25, linestyle='--')
		plt.axhline(y=chrom_boundry, color='black', lw=0.25, linestyle='--')

plt.xlim(0, int(args['--bins']))
plt.ylim(int(args['--bins']), 0)

ax = plt.gca()
ax.axes.xaxis.set_visible(False)
ax.axes.yaxis.set_visible(False)

# add sequence names to the plot
for i, chromosome in enumerate(genomefile_dict):
	middle_of_chromosome = genomefile_dict[chromosome][0] + ((genomefile_dict[chromosome][1] - genomefile_dict[chromosome][0]) / 2)
	middle_of_chromosome = middle_of_chromosome / bin_size
	chromosome_label = chromosome.split("_")[-1]
	if i % 2 == 0 and args["--stagger"] == "True":
		coord = int(args['--bins']) * -0.01
	else:
		coord = int(args['--bins']) * -0.05
	plt.text(middle_of_chromosome, coord, chromosome_label, ha='center', wrap=True, fontsize=int(args["--fsize"]))
	plt.text(coord, middle_of_chromosome, chromosome_label, ha='right', va='center', wrap=True, fontsize=int(args["--fsize"]))


plt.savefig("HiC_view.png", format="png", dpi=500, bbox_inches='tight')