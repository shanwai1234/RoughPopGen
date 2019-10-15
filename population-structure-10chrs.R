fh <- read.table('SAP.10chrs.zlorganized.maf005.vcf',sep='\t',head=T,check.names=FALSE)
sh <- fh
th <- sh[,-c(1:2)]
m <- t(th)
fit <- prcomp(m)
print (summary(fit))
pc <- fit$x[,c(1:20)]
print (summary(fit))
write.table(file='SAP-top20-10chrs.txt',pc)
