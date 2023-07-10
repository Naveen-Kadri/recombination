##grr  <- read.table ('/cluster/work/pausch/naveen/RECOMBINATION/LINKPHASE_RUN2/bv/cleaned_GRR.txt',head=T)

grr  <- snakemake@input[["infile"]]
outfile  <- snakemake@output[["outfile"]]
##remove outliers within sexes
out  <- data.frame ()
grr <- read.table (grr, head=T)
for (sex in 1:2) {
    mysex  <- grr [grr$sex ==sex,]
    n1 <- nrow (mysex)
    mysex  <- mysex [mysex$gp ==1 | mysex$noff >=6,]
    n2 <- nrow (mysex)
    
    mymean <- mean (mysex$GRR)
    mysd  <- sd(mysex$GRR)
    std <- abs ((mysex$GRR - mymean)/mysd)
    mysex  <- mysex [std <=4,]
    n3 <- nrow (mysex)
    cat ( sex, n1,n2,n3,'\n')
    out <- rbind (out, mysex[,c('gam', 'sex')])
}
write.table (out, outfile, col.names=T,row.names=F,quote=F,sep="\t")
