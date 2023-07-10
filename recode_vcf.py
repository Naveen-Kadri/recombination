import gzip
key_file = snakemake.input.key
vcf = snakemake.input.vcf

out_vcf = open (snakemake.output.vcf, 'w')



key ={}
print ('reading the key file')
with open (key_file) as inf:
    for line in inf:
        oid,nid=line.rstrip().split()
        key [oid] = nid
        
print ('recoding the vcf')
with gzip.open (vcf, "rt") as inf:
    for line in inf:
        if line[0:6]=="#CHROM":
            spl=line.rstrip().split()
            for i, el in enumerate (spl):
                if i>8:
                    try:
                        spl[i]=key [ spl[i] ]
                    except KeyError:
                        print(f'key not found for {spl[i]}')
            tw = "\t".join (spl)
            out_vcf.write (f'{tw}\n')
                
        else:
            out_vcf.write (f'{line}')
        
out_vcf.close ()

