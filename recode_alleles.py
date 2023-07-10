import gzip

# vcf_file='/cluster/work/pausch/naveen/RECOMBINATION/HAPBASED/CHR25/cleaned.vcf.gz'
# allele_file='ref_top_alleles.txt'
# out_file='todel.vcf'

vcf_file =snakemake.input.infile
allele_file =snakemake.input.allele_file
out_file=snakemake.output.outfile

vcf=gzip.open(vcf_file, "rt")
table=open (allele_file, "rt")

alleles=dict ()

for line in table:
    snp, ref, top, alok, ambi  =line.rstrip().split()
    if alok == "OK" and ambi =="OK":
        alleles[snp] = "".join ([ref[0],ref[1],top[0],top[1]])
    


comp = {
    "A" :"T",
    "T" :"A",
    "C" :"G",
    "G": "C"
}


flip = {
    "0/0" : "1/1",
    "1/1" : "0/0",
    "0/1" : "1/0",
    "1/0" : "0/1",
}




nprobs=0
out = open (out_file, "w")

seen=dict ()

for line in vcf:
    if line[0] == "#":
        out.write (line)
        #I Dont worry about the ids here
    else:
        spl=line.rstrip().split ("\t")
        snp=spl[2]
        ref=spl[3]
        alt=spl[4]
        myid=spl[0] + ":" + spl [1]

        if myid not in seen:
            seen [myid] =1
            if snp in alleles:    #ambigous variants and those wher ethe ref and top don't match are removed!
                alleleinfo=alleles.get (snp)
                ##change it in all cases !
                spl [3] = alleleinfo [0] #ref
                spl [4] = alleleinfo [1] #alt
                if ref == alleleinfo [0] or  ref == comp.get (alleleinfo [0]) : ##if the ref in vcf matches the real ref  or the compliment of the real ref then there is no need to flip the genotypes 
                    toprint = "\t".join (spl)                


                    out.write (f"{toprint}\n")
                elif ref == alleleinfo [1] or ref == comp.get (alleleinfo [1]):  #if the vcf ref matches the real alt or the compliment of the real alt then switch
                    toprint="\t".join(  [flip.get (el, el) for el in spl ])


                    out.write (f"{toprint}\n")
                else: #if there is not match in the alleles ; this should not happen, I have taken care of this while preparing the allele tables
                    nprobs+=1
                    print (f"Something's wrong at snp : {snp}\n")

            

out.close ()

print (f"Number of problems : {nprobs}")
