from collections import defaultdict

##infiles = [f'/cluster/work/pausch/naveen/RECOMBINATION/LINKPHASE_RUN1/fv/CHR{xchr}/nrec_hmm.txt' for xchr in range (1,30) ]
##outfile = 'ori_grr.txt'


infiles = snakemake.input.infiles
outfile = snakemake.output.outfile

out = open (outfile, 'w')
allinfo = defaultdict (dict)

infos = ['sex', 'mate', 'gp', 'noff', 'nhet', 'nhomo', 'ninfo', 'grr']
for infile in infiles:
    print (infile)
    with open (infile) as inf:
        for line in inf:
            off, par, sex, noff, gp,mate,nhet,nhomo,ninfo, grr=line.rstrip ().split()
            myid = (off, par)
            for i, info in enumerate(infos):
                if i <= 3:
                    allinfo [myid] [info] = globals().get (info)
                else:
                    allinfo [myid] [info] = allinfo [myid].get (info,0) + int ( globals().get (info))
            

header = "\t".join ( ['id', 'par'] + infos )
out.write (f'{header}\n')

for myid in allinfo:
    out.write (f'{myid[0]}\t{myid[1]}')
    for info in infos:
        out.write (f'\t{allinfo [myid] .get(info) }')
    out.write ('\n')

out.close()
