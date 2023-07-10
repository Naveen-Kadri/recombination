#infiles = [f'/cluster/work/pausch/naveen/RECOMBINATION/LINKPHASE_RUN1/fv/CHR{xchr}/linkphase.map' for xchr in range (1,30)]
#outfile = 'todel.txt'
infiles = snakemake.input.infiles
outfile = snakemake.output.outfile

out = open (outfile, 'w')

for infile in infiles:
    with open (infile, "rt") as inf:
        for line in inf:
            spl=line.rstrip().split()
        maxi=int (float(spl [2]) *1_000_000)
    out.write (f'{maxi}\n')

out.close()
