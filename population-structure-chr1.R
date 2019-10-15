fh <- read.table('SAP.10chrs.zlorganized.maf005.vcf',sep='\t',head=T,check.names=FALSE)
sh <- fh[fh$CHROM!='Chr01',]
th <- sh[,-c(1:2)]
m <- t(th)
fit <- prcomp(m)
pc <- fit$x[,c(1:20)]
write.table(file='SAP-top20-chr1.txt',pc)
