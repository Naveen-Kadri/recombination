'''
Run LINKPHASE to get the GRR counts

'''
import os
localrules: get_chrlengs, good_gams
LINKPHASE = "/cluster/home/nkadri/RECOMBINATION/LINKPHASE/LINKPHASE3/LINKPHASE3 "
linkinfile = "/cluster/home/nkadri/RECOMBINATION/LINKPHASE/linkin.txt"

OUT_DIR = "/cluster/work/pausch/naveen/RECOMBINATION"
plink = "/cluster/work/pausch/group_bin/plink "
breeds = ['bv', 'fv']
chromosomes = range(1, 30)
sexes = ['male', 'female']
runs = [1, 2]

wildcard_constraints:
    num = "\d"


# rule all:
#     input:
#         expand(
#             OUT_DIR + "/REFALT/LINKPHASE_RUN{run}/{breed}/cleaned_GRR.txt", run=runs, breed=breeds)

# rule all:
#     input:
#         expand(
#             OUT_DIR + '/REFALT/LINKPHASE_RUN{run}/{breed}/co_bins_{sex}.txt', run=runs, breed=breeds, sex=sexes)


rule linkphase1:
    input:
        typ = rules.format_for_linkphase.output.typ,
        map = rules.format_for_linkphase.output.map,
        pedigree = rules.recode_pedigree.output.pedigree
    output:
        outfile_nrec = OUT_DIR + \
            "/REFALT/LINKPHASE_RUN1/{breed}/CHR{chr}/nrec_hmm.txt",
        outfile_pos = OUT_DIR + \
            "/REFALT/LINKPHASE_RUN1/{breed}/CHR{chr}/recombinations_hmm",
        outfile_mcs = OUT_DIR + \
            "/REFALT/LINKPHASE_RUN1/{breed}/CHR{chr}/MCS.txt"
    params:
        lambda wildcards, input: os.path.dirname(input.typ)
    resources:
        mem_mb = 16000,
        walltime = "04:00"
    shell:
        " cd {params[0]} ; "
        " ln -fs {input.pedigree} linkphase.pedigree ;"
        " cp {LINKPHASE} . ;"
        " cp {linkinfile} . ;"
        " ./LINKPHASE3 "


rule clean_map_errors:
    input:
        mcs = rules.linkphase1.output.outfile_mcs,
        map = OUT_DIR + \
            "/REFALT/LINKPHASE_RUN1/{breed}/CHR{chr}/linkphase.map",
        typ = OUT_DIR + "/REFALT/LINKPHASE_RUN1/{breed}/CHR{chr}/linkphase.typ"
    output:
        map = OUT_DIR + \
            '/REFALT/LINKPHASE_RUN2/{breed}/CHR{chr}/linkphase.map',
        typ = OUT_DIR + '/REFALT/LINKPHASE_RUN2/{breed}/CHR{chr}/linkphase.typ'
    params:
        thresh = 0.99
    resources:
        mem_mb = 2000,
        walltime = "04:00"
    script:
        'clean_map_errors.py'

rule linkphase2:
    input:
        typ = rules.clean_map_errors.output.typ,
        map = rules.clean_map_errors.output.map,
        #pedigree = OUT_DIR + "/RECODE/{breed}/recoded.pedigree"
        pedigree = rules.recode_pedigree.output.pedigree
    output:
        outfile_grr = OUT_DIR + \
            "/REFALT/LINKPHASE_RUN2/{breed}/CHR{chr}/nrec_hmm.txt",
        outfile_pos = OUT_DIR + \
            "/REFALT/LINKPHASE_RUN2/{breed}/CHR{chr}/recombinations_hmm"
    params:
        lambda wildcards, input: os.path.dirname(input.typ)
    resources:
        mem_mb = 16000,
        walltime = "04:00"
    shell:
        " cd {params[0]} ; "
        " ln -fs {input.pedigree} linkphase.pedigree ;"
        " cp {LINKPHASE} . ;"
        " cp {linkinfile} . ;"
        " ./LINKPHASE3 "

rule count_GRR:
    input:
        infiles = expand(
            OUT_DIR + "/REFALT/LINKPHASE_RUN{{run}}/{{breed}}/CHR{chr}/nrec_hmm.txt", chr=chromosomes)
    output:
        outfile = OUT_DIR + "/REFALT/LINKPHASE_RUN{run}/{breed}/GRR.txt"
    script:
        "count_GRR.py"

rule clean_dco:
    ''' Remove recombinations that are possibly erros -- recombinations places within the thresh base pairs (see params) -- done chrwise to save time'''
    input:
        grr = OUT_DIR + "/REFALT/LINKPHASE_RUN{run}/{breed}/GRR.txt",
        copos = OUT_DIR + \
            "/REFALT/LINKPHASE_RUN{run}/{breed}/CHR{chr}/recombinations_hmm",
        map = OUT_DIR + \
            "/REFALT/LINKPHASE_RUN{run}/{breed}/CHR{chr}/linkphase.map"
    output:
        cleaned = OUT_DIR + \
            "/REFALT/LINKPHASE_RUN{run}/{breed}/CHR{chr}/cleaned_copos.txt",
        dcos = OUT_DIR + \
            "/REFALT/LINKPHASE_RUN{run}/{breed}/CHR{chr}/dco_positions.txt"
    params:
        thresh = 1_000_000
    resources:
        mem_mb = 2000,
        walltime = "00:20"
    script:
        "clean_dco.py"


rule cleaned_GRR:
    input:
        coposfiles = expand(
            OUT_DIR + "/REFALT/LINKPHASE_RUN{{run}}/{{breed}}/CHR{chr}/cleaned_copos.txt", chr=chromosomes),
        origrr = rules.count_GRR.output.outfile
    output:
        outfile = OUT_DIR + \
            "/REFALT/LINKPHASE_RUN{run}/{breed}/cleaned_GRR.txt"
    resources:
        mem_mb = 16000,
        walltime = "02:00"
    script:
        "cleaned_GRR.R"

rule get_chrlengs:
    input:
        infiles = expand(
            OUT_DIR + '/REFALT/LINKPHASE_RUN{{run}}/{{breed}}/CHR{chr}/linkphase.map', chr=chromosomes)
    output:
        outfile = OUT_DIR + \
            '/REFALT/LINKPHASE_RUN{run}/{breed}/chr_lengths.txt'
    script:
        'get_chrlengs.py'

rule good_gams:
    ''' Filter GRR data for design.. gp ==1 and noff >=6 and also remove outliers std >=4 '''
    input:
        infile = rules.cleaned_GRR.output.outfile
    output:
        outfile = OUT_DIR + '/REFALT/LINKPHASE_RUN{run}/{breed}/good_gams.txt'
    script:
        'good_gams.R'

rule make_co_bins:
    input:
        copos = expand(
            OUT_DIR + "/REFALT/LINKPHASE_RUN{{run}}/{{breed}}/CHR{chr}/cleaned_copos.txt", chr=chromosomes),
        chrlengs = rules.get_chrlengs.output.outfile,
        good_gams = rules.good_gams.output.outfile
    output:
        binfile = OUT_DIR + \
            '/REFALT/LINKPHASE_RUN{run}/{breed}/co_bins_{sex}.txt',
        numfile = OUT_DIR + \
            '/REFALT/LINKPHASE_RUN{run}/{breed}/co_bins_ngams_{sex}.txt'
    params:
        nbins = 100
    resources:
        mem_mb = 8000,
        walltime = '02:00'
    script:
        'make_co_bins.R'

rule plot_bins:
    input:
        copos = expand(
            OUT_DIR + "/REFALT/LINKPHASE_RUN{{run}}/{{breed}}/CHR{chr}/cleaned_copos.txt", chr=chromosomes),
        chrlengs = rules.get_chrlengs.output.outfile,
        good_gams = rules.good_gams.output.outfile
    output:
        plot_file = OUT_DIR + \
            "/REFALT/LINKPHASE_RUN{run}/{breed}/co_bins_{breed}.pdf"
    resources:
        mem_mb = 16000,
        walltime = "02:00"
    script:
        "plot_bins.R"
