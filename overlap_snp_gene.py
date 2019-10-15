import sys
from bx.intervals.intersection import Intersecter, Interval
import numpy as np
import operator

input_file = sys.argv[1]
check_file = sys.argv[2]
mm = open(input_file, 'r')
ref = open(check_file, 'r')
intv = int(sys.argv[3])

intersect_dict = {}

good_chrom = set([])
for x in range(1, 11):
    y = 'Chr' + str(x).zfill(2)
    intersect_dict[y] = Intersecter()
    good_chrom.add(y)

# this is the snp file in vcf format
ref.readline()
for y in ref:
    if y.startswith('#'):
        continue
    new = y.strip().split('\t')
    mychr = new[0]
    pos = new[1]
    st = int(new[1])
    sp = int(new[1])
    intersect_dict[mychr].add_interval(Interval(st, sp, value=str(st)))
# this is the output file
out = open(sys.argv[4], 'w')
# this is the annotation file
num = 0
xdict = {}
for line in mm:
    if line.startswith('##'):
        continue
    new = line.strip().split('\t')
    if new[2] != 'gene':
        continue
    x = new[-1].split(';')[0].split('.')
    gene = x[0].replace('ID=', '') + '.' + x[1]
    mychr = new[0]
    if mychr not in good_chrom:
        continue
    mystart = int(new[3]) - intv
    mystop = int(new[4]) + intv
    if mystart < 0:
        mystart = 0
    else:
        mystart = mystart
    a = intersect_dict[mychr].find(mystart, mystop)
    nlist = []
    if len(a) > 0:
        for b in a:
            nlist.append(b.value)
    num += len(nlist)
    if gene not in xdict:
        xdict[gene] = 0
    xdict[gene] = len(nlist)
    if len(nlist) == 0:
        out.write(gene + '\t' + new[0] + '\t' + new[3] + '\t' + new[4] + '\t' + 'NA' + '\n')
    else:
        out.write(gene + '\t' + new[0] + '\t' + new[3] + '\t' + new[4] + '\t' + ','.join(nlist) + '\n')

xlist = []
for m in sorted(xdict):
    xlist.append(xdict[m])
print num,np.median(xlist)
ref.close()
mm.close()
out.close()
