#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
from matplotlib.patches import Patch
import string
import sys
import math
import argparse
from Bio import SeqIO
import json
import gzip

def plothist(his, bin, tit, col):
	plot, ax1 = plt.subplots()
	bars1 = plt.bar(bin[:-1], his, width = 0.2, color = col)
	ax1.set_title(f'{tit}')
	ax1.set_xlabel('Reactivity')
	ax1.set_ylabel('Count')
	plt.grid()
	plt.tight_layout()
	plt.savefig(f'{tit}.png')
	ax1.set_ylim(0,50)
	plt.savefig(f'{tit}_zoom.png')
	plt.close()

def plotfreqs(re, pr, se, bi, dump = False):
	rcounts = {}
	pcounts = {}
	rnanless = []
	pnanless = []
	rc = 0
	pc = 0


	for r, p, s in zip(re, pr, se):
		if s not in rcounts:
			rcounts[s] = []
		if s not in pcounts:
			pcounts[s] = []
		if math.isnan(r) != True:
			rcounts[s].append(r)
			rnanless.append(r)
		if math.isnan(p) != True:
			pcounts[s].append(p)
			pnanless.append(p)

	phist, pbins = np.histogram(pnanless, bins = bi, range = [0, 3])
	print("Histogram for file 1:", sys.argv[1])
	print("Values:", phist)
	print("Edges:", pbins)
	rhist, rbins = np.histogram(rnanless, bins = bi, range = [0, 3])
	print("Histogram for file 2:", sys.argv[2])
	print("Values:", rhist)
	print("Edges:", rbins)



	plot, ax1 = plt.subplots()
	bars1 = plt.bar(rbins[:-1], rhist, width = 0.2, color = 'red')
	ax1.set_title('U Reactivity Distribution')
	ax1.set_xlabel('Reactivity')
	ax1.set_ylabel('Count')
	plt.grid()
	plt.tight_layout()
	plt.savefig('Udistro_bigbin.png')
	ax1.set_ylim(0, 150)
	plt.savefig('Udistrozoom_bigbin.png')
	plt.close()

	plot, ax2 = plt.subplots()
	bars2 = plt.bar(pbins[:-1], phist, width = 0.2, color = 'blue')
	ax2.set_title('pseudo-U Reactivity Distribution')
	ax2.set_xlabel('Reactivity')
	ax2.set_ylabel('Count')
	plt.grid()
	plt.tight_layout()
	plt.savefig('pseudoUdistro_bigbin.png')
	ax2.set_ylim(0, 150)
	plt.savefig('pseudoUdistrozoom_bigbin.png')
	plt.close()



	for s in pcounts:
		shist, sbins, = np.histogram(pcounts[s], bins = 15 ,range = [0,3])
		print(f'Histogram for {s} in {sys.argv[1]}')
		print("Values:", shist)
		print("Edges:", sbins)
		plothist(shist, sbins, f'Distribution_of_Reactivities_for_{s}_in_pseudo-U_data_bigbin', 'blue')

	for s in rcounts:
		shist, sbins, = np.histogram(rcounts[s], bins = 15 ,range = [0,3])
		print(f'Histogram for {s} in {sys.argv[2]}')
		print("Values:", shist)
		print("Edges:", sbins)
		plothist(shist, sbins, f'Distribution_of_Reactivities_for_{s}_in_U_data_bigbin', 'red')

def plot_seqs(re, pr, se):
	plot, axs = plt.subplots()
	leng = [i for i in range(len(se) - 1)]
	pl1 = axs.plot(leng, re, 'r-', label = "U")
	pl2 = axs.plot(leng, pr, 'b-', label = "pseudo-U")
	patch1 = Patch(facecolor ='red', label = "Missing U Data")
	patch2 = Patch(facecolor = 'blue', label = "Missing pseudo-U Data")


	axs.set_title(f'Raw Reactivities across Sequence')

	plt.xticks(range(1, len(se) + 1), se)
	plt.grid()

	plt.gca().margins(x=0)
	plt.gcf().canvas.draw()
	tl = plt.gca().get_xticklabels()
	maxsize = max([t.get_window_extent().width for t in tl])
	marg = 0.2 #inch margin
	s = maxsize/plt.gcf().dpi*len(se)+2*marg
	margin = marg/plt.gcf().get_size_inches()[0]

	plt.gcf().subplots_adjust(left = margin, right = 1.0 - margin)
	plt.gcf().set_size_inches(s, plt.gcf().get_size_inches()[1])
	axs.set_ylabel("P(Paired)")
	axs.set_xlabel(f'Reactivity')

	axs.legend(loc = 'upper center')

	for i in range(1,len(re)+1):
		if math.isnan(re[i - 1]):
			plt.axvspan(i, i + 1, color = 'red', alpha = 0.5)

	for i in range(1, len(pr)+1):
		if math.isnan(pr[i - 1]):
			plt.axvspan(i, i + 1, color = 'blue', alpha = 0.5)

	plt.savefig('reactsoverlaid.png')
	plt.show()



def get_filepointer(filename):
	"""
	Returns a filepointer to a file based on file name (.gz or - for stdin).
	"""

	fp = None
	if   filenadumpme.endswith('.gz'): fp = gzip.open(filename, 'rt')
	elif filename == '-':              fp = sys.stdin
	else:                              fp = open(filename)
	return fp

def read_fasta(filename):
	"""
	Simple fasta reader that returns name, seq for a filename.

	Parameters
	----------
	+ filename
	"""

	name = None
	seqs = []

	fp = get_filepointer(filename)

	while True:
		line = fp.readline()
		if line == '': break
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield(name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield(name, ''.join(seqs))
	fp.close()

pseudoreacts = sys.argv[1]
regreacts = sys.argv[2]
seqfile = sys.argv[3]
bins = sys.argv[4]

pfp = open(pseudoreacts)
rfp = open(regreacts)
sfp = open(seqfile)

preacts = []
reacts = []
seqs = []

for line1 in pfp:
	if line1.startswith('>'):
		continue
	preacts.extend(line1.strip().split(' '))
for line2 in rfp:
	if line2.startswith('>'):
		continue
	reacts.extend(line2.strip().split(' '))
for line3 in sfp:
	if line3.startswith('>'):
		continue
	seqs.extend(list(line3))

for i in range(len(preacts)):
	preacts[i] = float(preacts[i])

for j in range(len(reacts)):
	reacts[j] = float(reacts[j])



if sys.argv[5].find('r') != -1:
	plot_seqs(reacts, preacts, seqs)

if sys.argv[5].find('f') != -1:
	plotfreqs(reacts, preacts, seqs, int(bins))

if sys.argv[5].find('d') != -1:
	plotfreqs(reacts, preacts, seqs, int(sys.argv[6]), dump = True)
