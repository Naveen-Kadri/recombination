pedfile  <- snakemake@input[["pedigree"]]
unk  <- snakemake@params[["missing_code"]]
output_key  <- snakemake@output [["key"]]
output_ped  <- snakemake@output [["pedigree"]]


##pedfile  <-  '/cluster/work/pausch/naveen/RECOMBINATION/RECODE/bv/reduced.ped'
##unk  <-  '49c79968ae1ce80754d1c95a1bf30b1c'



ped <- read.table (pedfile, colClasses="character")
sum (ped [,2] %in% ped [,1])
missingsires <- ped [!ped [,2] %in% ped [,1] ,2  ]
missingdams <- ped [!ped [,3] %in% ped [,1] ,3  ]
missingsires  <- missingsires [missingsires != unk]
missingdams  <- missingdams [missingdams != unk]
missing <- c (missingsires,missingdams)
cat ("completing the pedigree\n")
if (length (missing) >0 ) {
    toadd <- data.frame(missing, 0,0)
    colnames (toadd) <- colnames (ped)
    ped <- rbind (toadd, ped)
    cat (nrow (toadd), "entries are added to complete the pedigree\n")
}else {
    cat ("Pedigree is complete\n")
}

colnames (ped) <-  c ("id", "sire", "dam", 'year')
ped$tmp <- 0  
ped$gen <- 0


ngen=-1
check=1


while (check==1){
      check=0
      ngen <- ngen+1
      cat ("Generation :", ngen, "\n")
      
      ##parents of the current geration are marked for gen++
      ped[  ped$id  %in%  (ped [ped$gen==ngen & ped$sire !=0, "sire"]), "tmp"] <- 1  #the id [first column] is a prent ? is he/she in second or third column? if he is then he blongs to the previous generation !
      ped[  ped$id  %in%  (ped [ped$gen==ngen & ped$dam  !=0, "dam" ]), "tmp"] <- 1
      
      if (any (ped$tmp==1))  {
          ped[ped$tmp==1, "gen"] <- ngen+1
          check=1
          ped$tmp=0;
      }
}

ped <- ped [order(ped[,  "gen"],  decreasing=T),  ]
#end
##plot (ped$year, ped$newcode)
ped$newcode <- 1:nrow (ped)
write.table (ped, output_ped, col.names=F,row.names=F,quote=F,sep="\t")
key <- ped [,c("id","newcode")]
##add the missing code
key  <- rbind (key, data.frame (id=unk, newcode=0))
write.table (key, output_key, col.names=F,row.names=F,quote=F,sep="\t")


