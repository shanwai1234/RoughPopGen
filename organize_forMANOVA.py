import sys

fh = open(sys.argv[1], 'r')
sh = open(sys.argv[2], 'r')
out = open(sys.argv[3], 'w')

mdict = {}
m = 0
for line in fh:
    new = line.strip().split('\t')
    if new[-1] == 'NA':
        continue
    nlist = new[-1].split(',')
    for i in nlist:
        nn = new[1] + '_' + i
        m += 1
        if nn not in mdict:
            mdict[nn] = []
        mdict[nn].append(new[0])
fh.close()

print m, len(mdict)

head = sh.readline().strip().split('\t')
out.write('SAP_Genes' + '\t' + head[0] + '\t' + '\t'.join(head[1:]) + '\n')
x = 0
for line in sh:
    new = line.strip().split('\t')
    snp = new[0] + '_' + new[1]
    if snp not in mdict:
        continue
    if len(mdict[snp]) == 1:
        out.write(mdict[snp][0] + '\t' + new[0] + '\t' + '\t'.join(new[1:]) + '\n')
    else:
        x += 1
print x
sh.close()
out.close()
