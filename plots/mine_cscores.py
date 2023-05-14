import gzip
import sys

def readscores(fn, colnames = True):
    sc = {}
    f = get_filepointer(fn)


    for line in f:

        if colnames == True:
            names = line.strip().split(' ')
            colnames = False
            for name in names:
                sc[name] = []
            continue

        elif len(sc) == 0:
            data = line.strip().split(' ')
            for i in range(len(data)):
                key = f'col_{i}'
                sc[key] = []

        data = line.strip().split(' ')
        for val, col in zip(data, sc):
            sc[col].append(val)


    return sc


def get_filepointer(filename):
	"""
	Returns a filepointer to a file based on file name (.gz or - for stdin).
	"""

	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	elif filename == '-':          fp = sys.stdin
	else:                          fp = open(filename)
	return fp

scorefile = sys.argv[1]

scores = readscores(scorefile)

print('> M-seq')

cscs = []

for csc in scores['c-score']:
    cscs.append(csc)

print(' '.join(cscs))
