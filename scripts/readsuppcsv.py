import sys


csvfile = sys.argv[1]

pos = []
seq = []
shape = []
shape_sd = []
DMS = []
DMS_sd = []
CMCT = []
CMCT_sd = []

cols = [
    pos,
    seq,
    shape,
    shape_sd,
    DMS,
    DMS_sd,
    CMCT,
    CMCT_sd
]

nameline = 'M-seq'

fp = open(csvfile, 'r')

names = fp.readline().strip().split(',')

for line in fp:
    line = line.strip().split(",")
    print(line)

    for i in range(len(line)):
        num = False
        entry = None
        line[i] = line[i].strip()
        entry = line[i]

        if line[i] == '-':
            entry = "nan"
        elif entry == '-999':
            entry = "nan"
        elif entry.startswith('-'):
            entry = '0'

        cols[i].append(entry)

print(cols)
for i in range(1, len(cols)):
    pp = open(f'{names[i]}_{i}.fa', 'wt')
    pp.write(f'> {nameline}\n')







    for j in range(0, len(cols[i]), 80):
        linewrite = ' '.join(cols[i][j:j+80])

        pp.write(f'{linewrite}\n')

    pp.close()
