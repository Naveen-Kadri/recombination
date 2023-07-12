''' Assumes the recombinations_hmm is sorted for gametes and then for CO positions -- this is the case so it works
'''

# grr_file = 'nrec_hmm.txt'
# coposfile = 'recombinations_hmm'
# mapfile = 'linkphase.map'
# outfile = 'todel.txt'


grr_file =  snakemake.input.grr
coposfile = snakemake.input.copos
mapfile = snakemake.input.map
outfile_cleaned = snakemake.output.cleaned
outfile_dcos = snakemake.output.dcos
thresh = snakemake.params.thresh

out_clean = open (outfile_cleaned, 'w')
out_dcos = open (outfile_dcos, 'w')
pos = {}
infos = {}

with open (grr_file, "rt") as inf:
    for line in inf:
        off, par, sex, noff, gp,mate,nhet,nhomo,ninfo, grr=line.rstrip ().split()
        gam = off + ":" + par
        info = [sex, gp, noff, mate, grr]
        infos  [gam] = info


with open (mapfile) as inf:
    for line in inf:
        xnum, xsnp, xpos, xchr=line.rstrip().split()
        pos [xnum] = float(xpos) * 1_000_000


prev_gam = 'NA'
prev_end = thresh *-2
nall = 0
nrmvd = 0
tormv = {}
with open (coposfile) as inf:
    for coindex, line in enumerate(inf):
        nall+=1
        myid, mypar, start, end = line.rstrip().split()
        gam = myid + ":" +  mypar
        start = pos [start]
        end = pos [end]
        dist = start - prev_end
        if prev_gam == gam and dist <= thresh:
            nrmvd +=1
            tormv[coindex-1] =1
            tormv[coindex] =1 
            #print (f'Removed a co for gamete {gam}...gap between the consecutive co {prev_end} = {start} = {dist} ')
        prev_gam = gam
        prev_end = end

        
print (f'Number of total COs, to be removed : {nall},{nrmvd}')

header = "\t".join (['gam', 'id', 'par', 'st', 'en', 'start', 'end', 'chr', 'leng', 'sex', 'gp', 'noff', 'mate', 'origrr'])
out_clean.write (f'{header}\n')
out_dcos.write (f'{header}\n')
with open (coposfile, "rt") as inf:
    for coindex,line in enumerate(inf):

        myid, mypar, start, end = line.rstrip().split()
        gam = myid + ":" + mypar
        leng = str(pos [end] - pos [start] )
        tw = [gam, myid, mypar, start, end, str(pos [start]), str(pos[end]), xchr, leng]
        tw = "\t".join (  tw + infos [gam] )
        if coindex not in tormv:
            out_clean.write (f'{tw}\n')
        else:
            out_dcos.write (f'{tw}\n')

out_clean.close()
out_dcos.close()
            
            

        
    
