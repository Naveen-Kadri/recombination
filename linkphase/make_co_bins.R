copos_files  <- snakemake@input[["copos"]]
good_gams  <- snakemake@input[["good_gams"]]
chrlengs_file  <- snakemake@input[["chrlengs"]]
sex <- snakemake@wildcards[["sex"]]
nbins <- snakemake@params[["nbins"]]
binfile<- snakemake@output[["binfile"]]
numfile<- snakemake@output[["numfile"]]



if (sex == 'male') {
    sex <- 1
}else{
    sex <- 2
}

#copos_files  <- paste0('/cluster/work/pausch/naveen/RECOMBINATION/REFALT/LINKPHASE_RUN2/bv/CHR',1:29,'/cleaned_copos.txt')
#chrlengs_file  <- '/cluster/work/pausch/naveen/RECOMBINATION/REFALT/LINKPHASE_RUN2/bv/chr_lengths.txt'
#good_gams  <- '/cluster/work/pausch/naveen/RECOMBINATION/LINKPHASE_RUN2/bv/good_gams.txt'


chrlengths <- scan (chrlengs_file)
goodgams  <- read.table (good_gams, head=T)


####----------------

RES <-  list ()
ns  <- c()
gams  <- c()
tokeep  <- goodgams [goodgams$sex == sex,1]
result <- data.frame (matrix (0,nrow=29, ncol=nbins)) 
for (chr in 1:29) {
    mychr  <- read.table (copos_files [chr], head=T)
    ##remove outliers
    mychr  <- mychr [mychr$gam %in% tokeep,]
    ##fiter for sex
    ##mychr  <- mychr [mychr$sex == sex,]
    gams  <- unique (c(gams, mychr$gam))
    myseq <- seq (0, chrlengths[chr], length.out=nbins+1)
    bin.size <- chrlengths [chr] /nbins
    for (i in 1:nbins) {
        mybin <- mychr[mychr$start <= myseq [i+1] & mychr$end > myseq [i] ,]
        if ( nrow (mybin) >0 ) {
            starts <- pmax (myseq [i], mybin$start)
            ends <- pmin (myseq [i+1], mybin$end)
            result [chr,i] <-  sum ((ends - starts +1)/mybin$leng )   
        }
    }
    err <- nrow (mychr) - sum (result [chr,])
    cat (sex,chr , err , "\n")
}

write.table (result, binfile, col.names=F, row.names=F,quote=F,sep="\t")
n<- length (unique (gams))
write.table  (n,numfile, col.names=F,row.names=F,quote=F,sep="\t")

