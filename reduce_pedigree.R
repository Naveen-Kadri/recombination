pedfile  <- snakemake@input[["pedigree"]]
genophenofile <- snakemake@input [["genotyped"]]
unk  <- snakemake@params[["missing_code"]]
outfile  <- snakemake@output [["outfile"]]


cat ("Reading the pedigree file\n")
ped <- read.table (pedfile,  colClasses='character')
cat ('number of lines in the pedfile ', nrow (ped), "\n")



cat ("Readind the list of animals to be kept in the reduced pedigree\n")
genopheno <- scan (genophenofile, what='character', sep="\t")

cat ("The unknowns are coded as ",  unk,  "\n")

add.parents <- function (x) {
  parents <-c(ped [ped[, 1]%in%x,  2],  ped [ped[, 1]%in%x, 3]) #sires and dams
  parents <- parents [parents != unk] #only known
  withanc <- unique(c(x,  parents)) #given a vector gives the same vector with parents
  return (withanc)
}

base <- genopheno #assumes that all genopheno are in the pedigree which is nto the case atleast in the CRV pedigree !
cat ("Pegigree needs to be built for ",  length (genopheno),  "animals\n")
withanc <- add.parents (base)

i <- 1
while ( length(withanc) != length (base)  ) {
    cat ("---------------------------------------------------\n")
    cat ("iteration:",  i, "\n")
    base <- withanc
    withanc <- add.parents (base)
    cat (length (withanc),  "||", 100* length (withanc) /nrow (ped),   " % %have to be kept in the pedigree\n")
    cat ("---------------------------------------------------\n")
    i <- i+1
}

myped <-ped[ped [, 1] %in% base , ] #assumes that all pedigreed animals are on the first column !
cat ("Pedigree reduced from :",   nrow(ped),  "to :",  nrow (myped),  "\n\n\n")
nf <- withanc [!withanc %in% ped [, 1]]

cat ("Pedigree not found for ::", length (nf), "\n")

write.table (myped,  outfile, col.names=F, row.names=F, quote=F, sep="\t")

cat ("Reduced pedigree is written to Reduced.ped", "\n\n\n")
