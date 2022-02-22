#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Usage: busco2synteny.py -A <STR> -a <STR> -B <STR> -b <STR> [-l <STR> -m <INT> -n <INT> -c <INT> -w <INT> -t <FLT> -h]

  [Options]
    -A, --genomefileA <STR>                     Genomefile for taxon A
    -B, --genomefileB <STR>                     Genomefile for taxon B
    -a, --tsvA <STR>                            BUSCO tsv for taxon A
    -b, --tsvB <STR>                            BUSCO tsv for taxon B
    -l, --labels <STR>                          Whether to plot labels, choose from True or False [default: True]
    -m, --gapA <INT>                            Gap between genome A chromosomes [default: 10_000_000]
    -n, --gapB <INT>                            Gap between genome B chromosomes [default: 10_000_000]
    -c, --chromosome_width <INT>                Chromosome width [default: 6]
    -w, --alignment_width <INT>                 Alignment width [default: 1]
    -t, --alpha <FLT>                           Alpha of alignments [default: 0.1]
    -h, --help                                  Show this message

"""

import sys
from docopt import docopt
import collections
from scipy.interpolate import make_interp_spline, BSpline
import numpy as np
import matplotlib.pyplot as plt

def generate_genomefile_dict(genomefile, offset):
	genomefile_dict = {}
	cumulative_genome = offset
	with open(genomefile, "r") as fin:
		# for each chromosome, record cumulative coordinates, orientation, and label
		for line in fin:
			line = line.rstrip()
			chromosome, chromosome_length, orientation, label = line.split("\t")
			chromosome_length = int(chromosome_length)
			genomefile_dict[chromosome] = [cumulative_genome, cumulative_genome + chromosome_length, orientation, label]
			cumulative_genome += chromosome_length
			cumulative_genome += offset
	return genomefile_dict

# plots chromosomes as lines and adds labels if arg is True
def plot_chromosomes(genomefile_dict, y_coord, labels, chromosome_width):
	for chromosome in genomefile_dict:
		plt.plot([genomefile_dict[chromosome][0], genomefile_dict[chromosome][1]], 
			[y_coord, y_coord], color='slategrey', alpha=1, linewidth=chromosome_width)
		middle_of_chromosome = genomefile_dict[chromosome][0] + \
		((genomefile_dict[chromosome][1] - genomefile_dict[chromosome][0]) / 2)
		if labels == "True":
			plt.text(middle_of_chromosome, y_coord*1.15, genomefile_dict[chromosome][3], ha='center', va='center', 
				wrap=True, fontsize=10)

def generate_BUSCO_dicts(tsv, genomefile_dict, y_coord, BUSCO_coords):
	with open(tsv, "r") as fin:
		# for each BUSCO, record coordinates in either genome
		# the y_coord is just is so we can double check later that the entries are from different genomes
		for line in fin:
			line = line.rstrip()
			BUSCO_ID = line.split("\t")[0]
			chromosome = line.split("\t")[2]
			mid_point = int(line.split("\t")[3]) + ((int(line.split("\t")[4]) - int(line.split("\t")[3])) / 2)
			# only interested in BUSCOs on sequences in the genomefile
			if chromosome in genomefile_dict.keys():
				# flip coords if orientation is -
				if genomefile_dict[chromosome][2] == '-':
					mid_point = genomefile_dict[chromosome][1] - mid_point
				if genomefile_dict[chromosome][2] == '+':
					mid_point += genomefile_dict[chromosome][0]
				BUSCO_coords[BUSCO_ID].append(mid_point)
				BUSCO_coords[BUSCO_ID].append(y_coord)
	return BUSCO_coords


args = docopt(__doc__)

# generate dicts for each genome with cumulative coordinates
genomefile_A_dict = generate_genomefile_dict(args['--genomefileA'], int(args['--gapA']))
genomefile_B_dict = generate_genomefile_dict(args['--genomefileB'], int(args['--gapB']))

# set up plot
fig = plt.figure(figsize=(28,4), frameon=False)
ax = fig.add_subplot(111)
ax.axis('off')
plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
plt.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)

# each BUSCO has a list of coordinates
BUSCO_coords = collections.defaultdict(list)
# genome A is always plotted with y coordinate = 1
# genome B is always plotted with y coordinate = -1
BUSCO_coords = generate_BUSCO_dicts(args['--tsvA'], genomefile_A_dict, 1, BUSCO_coords)
BUSCO_coords = generate_BUSCO_dicts(args['--tsvB'], genomefile_B_dict, -1, BUSCO_coords)

# for each BUSCO, if the two entries are from the two genomes, plot a curve between them
for BUSCO in BUSCO_coords:
	if len(BUSCO_coords[BUSCO]) == 4 and BUSCO_coords[BUSCO][1] != BUSCO_coords[BUSCO][3]:
		# code is adapted from https://gist.github.com/andrewgiessel/5684769
		lower_x = BUSCO_coords[BUSCO][2]
		upper_x = BUSCO_coords[BUSCO][0]
		Y = np.array([-1, -0.568, -0.32, -0.16, -0.056, 0, 0.056, 0.16, 0.32, 0.568, 1])
		X = np.array([lower_x, lower_x + (0.1*(upper_x - lower_x)), 
		lower_x + (0.2*(upper_x - lower_x)), lower_x + (0.3*(upper_x - lower_x)), 
		lower_x + (0.4*(upper_x - lower_x)), lower_x + (0.5*(upper_x - lower_x)), 
		lower_x + (0.6*(upper_x - lower_x)), lower_x + (0.7*(upper_x - lower_x)), 
		lower_x + (0.8*(upper_x - lower_x)), lower_x + (0.9*(upper_x - lower_x)), upper_x])
		# requires sorted arrays, so flip if in the wrong order
		if lower_x > upper_x:
			X = np.flip(X)
			Y = np.flip(Y)
		xnew = np.linspace(X.min(), X.max(), 300) 
		spl = make_interp_spline(X, Y, k=3)
		power_smooth = spl(xnew)
		plt.plot(xnew, power_smooth, color='slategrey', alpha=float(args['--alpha']), linewidth=int(args['--alignment_width']))

# plot the chromosomes
plot_chromosomes(genomefile_A_dict, 1, args["--labels"], int(args['--chromosome_width']))
plot_chromosomes(genomefile_B_dict, -1, args["--labels"], int(args['--chromosome_width']))

plt.savefig("busco2synteny.pdf", format="pdf", bbox_inches='tight')

"""
# example command
~/software/busco2synteny.py -A brenthis_ino.SP_BI_364.v2_0.sequences.chromosomes.genomefile -B melitaea_cinxia.genomefile -a brenthis_ino.SP_BI_364.v2_0.sequences.BUSCO.shared.tsv -b melitaea_cinxia.BUSCO.shared.tsv -l False -m 27_250_000 -t 0.1
"""