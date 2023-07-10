'''
Check for Mendelian errors
Errors are set to missing -- for all members of a trio for the SNP showing error.
Trios and Duos are tested.
MISSING GENOTYPES -- If genotype is missing for one parent in a trio it is considered as a duo and opposite homozygosity is tested 
'''

import gzip

##vcf = '/cluster/work/pausch/naveen/RECOMBINATION/REFALT/bv/CHR25/recoded.vcf.gz'
##pedigree_file = '/cluster/work/pausch/naveen/RECOMBINATION/REFALT/bv/recoded.pedigree'
##log_file ='test.log'
##out_file ='test.vcf'

vcf = snakemake.input.vcf
pedigree_file = snakemake.input.pedigree
log_file = snakemake.output.log
out_file = snakemake.output.vcf



log =open (log_file, "w")
header =  "\t".join ( ['rel', 'position', 'off', 'sire', 'dam', 'geno_off', 'geno_sire', 'geno_dam'] )
log.write (f'{header}\n')
out = open (out_file, "w")


additive = {    "0/0" :0,    "0|0" :0,    "1/1" :2,    "1|1" :2,    "0/1" :1,    "0|1" :1,    "1/0" :1,    "1|0" :1 ,    "./." : -9,    ".|." : -9 }
vcfgeno = {    0  : '0/0',    1  : '0/1',    2  : '1/1',   -9  :  './.'}
recoded = {    (0,0) : 0,    (0,1) : 1,    (1,0) : 1,    (1,1) : 2
}
table = {}
genotypes = [(0,0), (0,1), (1,1) ]
for f in genotypes:
    for m in genotypes:
        mend =[]
        for fal in f:
            for mal in m:
                mend.append (   recoded.get (  (fal, mal) ))
        table [ (recoded.get (f), recoded.get (m))       ] = set(mend)
                
for mykey, myvalue in table.items ():
    print (mykey, myvalue)
print ("scanning the genotype file")

trios ={}
ntrio=0
nvar =0
merr=0
with gzip.open (vcf, "rt") as inf:
    for line in inf:
        if line[0:2]!="##":
            spl = line.rstrip().split()
            if line[0:6] == "#CHROM":
                ids = spl [9:]
                genotyped = {el:1 for el in ids}
                out.write (f'{line}')
                print ('reading the pedigree file')
                with open (pedigree_file) as ped:
                    for line in ped:
                        i,s,d =line.rstrip().split() [0:3]
                        if i in genotyped:
                            if s in genotyped or d in genotyped: ##WHEN AT LEAST ONE PARENT IS GENOTYPED.
                                trios [ (i,s,d) ] = 1
                                ntrio+=1
                                #if ntrio>10:break
                print (f'Number of relationships = {ntrio}')
            else:
                nvar +=1
                pos = spl [1]
                gtindex = spl [8].split (":").index ('GT')
                gts = spl [9:]
                gts = [additive.get ( gt.split (":") [gtindex], gt.split(":") [gtindex] ) for gt in gts]
                genotypes = dict (zip (ids, gts))
                for mytrio in trios:
                    myid, mysire, mydam = mytrio
                    mygenotypes=list (map(     lambda el : genotypes.get (el, -9), mytrio))
                    if mygenotypes [0] != -9: ## OFFSPRING'S GENOTYPE SHOULD ALWAYS BE KNOWN
                        if -9 not in mygenotypes: ##ALL SHOULD BE NON-MISSING! > TRIO
                            exp = table [ ( mygenotypes [1], mygenotypes [2]  ) ]
                            if mygenotypes [0] not in exp:
                                merr+=1
                                log.write (f'trio\t{pos}\t{mytrio[0]}\t{mytrio[1]}\t{mytrio[2]}\t{mygenotypes[0]}\t{mygenotypes[1]}\t{mygenotypes[2]}\n')
                                for sample in mytrio:
                                    genotypes [sample] =-9
                        elif 0 in mygenotypes and 2 in mygenotypes: #OPPOSITE HOMOZYGOTES
                            merr+=1
                            log.write (f'duo\t{pos}\t{mytrio[0]}\t{mytrio[1]}\t{mytrio[2]}\t{mygenotypes[0]}\t{mygenotypes[1]}\t{mygenotypes[2]}\n')
                            for sample in mytrio:
                                genotypes [sample] =-9
                genotypes = [vcfgeno.get (genotypes[myid])  for myid in ids] ##IMPORTANT TO KEEP THE ID ORDER
                tw="\t".join( spl[:9] + genotypes )
                out.write (f'{tw}\n')
                print (nvar, merr)
        else:
            out.write (f'{line}')

out.close()
log.close()
