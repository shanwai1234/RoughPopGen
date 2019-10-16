import sys

fh = open(sys.argv[1],'r')
mp = open(sys.argv[2],'w')
geno = open(sys.argv[3],'w')

head = fh.readline()
mp.write('SNP\tChr\tPos\n')
mdict = {}
sdict = {}
plist = []
for line in fh:
	new = line.strip().split('\t')
	snp = new[1].split('_')
	chrom = snp[0]
	pos = snp[1]
	if new[1] not in mdict:
		mdict[new[1]] = new[2:]
	if chrom not in sdict:
		sdict[chrom] = []
	sdict[chrom].append(int(pos))

for i in range(1,11):
	chrom = 'Chr'+str(i).zfill(2)
	for j in sorted(sdict[chrom]):
		snp = chrom+'_'+str(j)
		mp.write(snp+'\t'+str(i)+'\t'+str(j)+'\n')
		nlist = []
		for k in mdict[snp]:
			nlist.append(str(int(float(k))))
		geno.write('\t'.join(nlist)+'\n')
fh.close()
mp.close()
geno.close()
