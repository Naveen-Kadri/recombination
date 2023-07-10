pedfile = snakemake.input.ped
mapfile = snakemake.input.map


out_map = open(snakemake.output.map, 'w')
out_typ = open (snakemake.output.typ, 'w')

with open (pedfile, "rt") as inf:
    for line in inf:
        spl=line.rstrip().split()
        tw = " ".join ( [ spl [0] ] + spl [6:] )
        out_typ.write (f'{tw}\n')
out_typ.close()

print ('formatting the mapfile')
nm=0
with open (mapfile, "rt") as inf:
    for line in inf:
        nm+=1
        xchr, name, link, pos=line.rstrip().split()
        pos = float(pos)/1_000_000
        #tw = " ".join (  [nm, name, pos, xchr])
        #out_map.write (f'{tw}\n')
        out_map.write (f'{nm}\t{name}\t{pos}\t{xchr}\n')
        
out_map.close()
        
        
