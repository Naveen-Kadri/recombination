
# xchr = 1
# mcs_file = f'/cluster/work/pausch/naveen/RECOMBINATION/LINKPHASE1/bv/CHR{xchr}/MCS.txt'
# map_file = f'/cluster/work/pausch/naveen/RECOMBINATION/LINKPHASE1/bv/CHR{xchr}/linkphase.map'
# typ_file=f'/cluster/work/pausch/naveen/RECOMBINATION/LINKPHASE1/bv/CHR{xchr}/linkphase.typ'
# out_map_file ='todel.map'
# out_typ_file ='todel.typ'
# thresh = 0.99

mcs_file = snakemake.input.mcs
map_file = snakemake.input.map
typ_file = snakemake.input.typ

out_map_file = snakemake.output.map
out_typ_file = snakemake.output.typ

thresh = snakemake.params.thresh


out_map = open (out_map_file, "w")
out_typ = open (out_typ_file, "w")


tormv = {}
with open (mcs_file) as inf:
    for i, line in enumerate (inf):
        if i>0:
            spl = line.rstrip().split()
            if float(spl [8] ) < thresh:
                ##when mcs is low.. the snp and the next should be removed.
                mysnpnum = int (spl [0])
                tormv [mysnpnum] =1
                tormv [mysnpnum+1] =1
#print (f'Number of markers to remove : {len (tormv) }')
#print ('max in tormv')
#print (max (list(tormv.keys()) ))
nkept = 0

nsnp=0
print (f'Removing bad markers from the map file')
with open (map_file) as inf:
    for line in inf:
        spl=line.rstrip().split()
        nsnp +=1
        if nsnp not in tormv:
            nkept+=1
            spl [0] = str(nkept)
            tw = "\t".join (spl)
            out_map.write  (f'{tw}\n')
print (f'Number of markers in the map, kept : {nsnp}, {nkept}')


print ('Reading the typ file')
rmvd = []
with open (typ_file, "rt") as inf:
    for i,line in enumerate(inf):
        ###--------------
        ##if i>0:break
        ###--------------
        genotypes=line.rstrip().split()
        myid = genotypes.pop(0)
        out_typ.write (f'{myid}')
        for snpnum, mycol in enumerate (range (0, nsnp*2, 2)):
            if snpnum+1 not in tormv:
                out_typ.write (f' {genotypes [mycol] } {genotypes[mycol+1]}') #zero based
            else:
                rmvd.append(snpnum+1)
        out_typ.write ('\n')
print (f'Number of markers : {nsnp}')


out_map.close()                
out_typ.close()
