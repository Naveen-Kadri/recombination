##rm (list = ls (all=T))

#coposfiles  <- paste0('/cluster/work/pausch/naveen/RECOMBINATION/LINKPHASE_RUN1/fv/CHR',1:29,'/cleaned_copos.txt')
#origrrfile  <- '~/RECOMBINATION/LINKPHASE/ori_grr.txt'
#outfile  <- 'xx'



coposfiles  <- snakemake@input[["coposfiles"]]
origrrfile  <-  snakemake@input[["origrr"]]
outfile  <- snakemake@output[["outfile"]]


for (chr in 1:length (coposfiles) ) {
    inf  <- read.table (coposfiles[chr], head=T)
    mytab  <- aggregate (inf$gam, by=list (inf$gam), length)
    colnames (mytab)  <- c('gam', paste0('chr',chr))

    if (chr ==1) {
        genome <- mytab
    }else {
        genome  <- merge (genome, mytab, by='gam', all=T)
        genome[is.na(genome)]  <- 0
    }
    
    cat (chr, nrow (genome), "\n")
}


genome$GRR  <- apply (genome [,-c(1)],1,sum, na.rm=T)
genome  <- genome [,c("gam", "GRR")]


origrr  <- read.table (origrrfile, head=T)
colnames (origrr) [ colnames(origrr) == 'grr' ]  <- 'origrr'
origrr$gam  <- paste(origrr$id, origrr$par,sep=":")
out  <- merge (origrr, genome, by='gam')


##add the std dev to the file
mymean <- mean (out$GRR)
mysd <- sd (out$GRR)

out$zscore <- abs (out$GRR - mymean) /mysd

write.table (out, outfile, col.names=T,row.names=F,quote=F,sep="\t")





