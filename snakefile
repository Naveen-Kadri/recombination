# WILDCARDS
breeds = ['fv', 'bv']
chromosomes = range(1, 30)
runs = [1, 2]
sexes = ['male', 'female']
plink_binaries = [
    f'/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/asr_transfer/{{breed}}/geno2.{ext}' for ext in ['bed', 'bim', 'fam']]
pedigree_file = '/nfs/nas12.ethz.ch/fs1201/green_groups_tg_public/data/BTA/asr_transfer/{breed}/geno.ped'
allele_table = 'ref_top_alleles.txt'


OUT_DIR = '/cluster/work/pausch/naveen/RECOMBINATION'
# HOW MISSING INFO IS CODED IN THE FILE
missing_codes = "-9,49c79968ae1ce80754d1c95a1bf30b1c"

# PROGRAMS
plink = "/cluster/work/pausch/group_bin/plink "

# PARAMETERS TO CLEAN [PLINK]
mind = 0.1
geno = 0.1
hwe = 1e-08


rule all:
    input:
        expand(
            OUT_DIR + '/REFALT/LINKPHASE_RUN{run}/{breed}/co_bins_{breed}.pdf', run=runs, breed=breeds),
        expand(
            OUT_DIR + '/REFALT/LINKPHASE_RUN{run}/{breed}/cleaned_GRR.txt', run=runs, breed=breeds),
        expand(
            OUT_DIR + '/REFALT/LINKPHASE_RUN{run}/{breed}/co_bins_{sex}.txt', run=runs, breed=breeds, sex=sexes)


rule get_vcf:
    ''' Get vcf .. simultaneosly also clean for individual-wise missingness and HWE [missingness is cleaned SNP-WISE first, I dont't want that .. so cleaning is doen in two steps]'''
    input:
        infiles = plink_binaries,
    output:
        outfile = OUT_DIR + '/REFALT/{breed}/genotypes.vcf.gz'
    params:
        lambda wildcards, input: input.infiles[0][:-4],
        lambda wildcards, output: output.outfile[:-7],
    resources:
        mem_mb = 16000,
        walltime = "04:00"
    threads:
        1
    shell:
        " {plink} --cow "
        " --chr 1-29 "
        " --missing_code {missing_codes} "
        " --bfile {params[0]} "
        " --geno {geno}"
        " --hwe {hwe}"
        " --keep-allele-order "
        " --recode vcf-iid bgz "
        " --out {params[1]} "
        " --threads {threads}"

rule recode_alleles:
    input:
        infile = rules.get_vcf.output.outfile,
        allele_file = allele_table
    output:
        outfile = OUT_DIR + '/REFALT/{breed}/ref_alt_genotypes.vcf'
    script:
        "recode_alleles.py"

rule zip:
    input:
        infile = rules.recode_alleles.output.outfile
    output:
        outfile = OUT_DIR + '/REFALT/{breed}/ref_alt_genotypes.vcf.gz'
    shell:
        'module load htslib ; bgzip {input.infile}'


rule clean:
    ''' get the ped format and simultaneosly clean for Individual wise missingness'''
    input:
        infile = rules.zip.output.outfile
    output:
        outfile = OUT_DIR + '/REFALT/{breed}/cleaned.vcf.gz'
    params:
        lambda wildcards, output: output.outfile[:-7],
        mind = mind
    resources:
        mem_mb = 16000,
        walltime = "04:00"
    threads:
        1
    shell:
        " {plink} --cow "
        " --vcf {input.infile}"
        " --mind {params.mind}"
        " --keep-allele-order "
        " --recode vcf-iid bgz "
        " --out {params[0]} "
        " --threads {threads}"

rule split:
    '''
    split into chromosomewise files
    '''
    input:
        infile = rules.clean.output.outfile
    output:
        outfile = OUT_DIR + \
            "/REFALT/{breed}/CHR{chr}/cleaned_ref_alt_genotypes.vcf.gz"
    params:
        lambda wildcards, output:  output.outfile[:-7],
    resources:
        mem_mb = 4000,
        walltime = '00:30'
    threads:
        1
    shell:
        " {plink} --cow "
        " --chr {wildcards.chr}"
        " --missing_code {missing_codes} "
        " --vcf {input.infile}"
        " --recode vcf-iid bgz "
        " --out {params[0]} "
        " --threads {threads}"
        " --keep-allele-order"

rule list_genotyped:
    input:
        infiles = plink_binaries
    output:
        outfile = OUT_DIR + "/REFALT/{breed}/genotyped.txt"
    params:
        lambda wildcards, input: [
            myfile for myfile in input.infiles if ".fam" in myfile][0]
    shell:
        "cut -d' ' -f2 {params[0]} > {output.outfile}"


rule reduce_pedigree:
    '''
    Linkphase uses sample ids to get the number of samples and has to be numeric!, so to reduce linkphase memory only useful info in pedigree is retained  
    '''
    input:
        pedigree = pedigree_file,
        genotyped = rules.list_genotyped.output.outfile
    output:
        outfile = OUT_DIR + "/REFALT/{breed}/reduced.ped"
    params:
        missing_code = "49c79968ae1ce80754d1c95a1bf30b1c"
    script:
        "reduce_pedigree.R"

rule sort_pedigree:
    input:
        pedigree = rules.reduce_pedigree.output.outfile
    output:
        pedigree = OUT_DIR + "/REFALT/{breed}/sorted.ped",
        key = OUT_DIR + "/REFALT/{breed}/id_key.txt"
    params:
        missing_code = "49c79968ae1ce80754d1c95a1bf30b1c"
    script:
        "sort_pedigree.R"

rule recode_vcf:
    input:
        key = rules.sort_pedigree.output.key,
        vcf = rules.split.output.outfile,
    output:
        vcf = OUT_DIR + '/REFALT/{breed}/CHR{chr}/recoded.vcf',
    resources:
        mem_mb = 3000,
        walltime = "00:15"
    script:
        'recode_vcf.py'


rule zip2:
    input:
        infile = rules.recode_vcf.output.vcf
    output:
        vcf = OUT_DIR + '/REFALT/{breed}/CHR{chr}/recoded.vcf.gz'
    shell:
        'module load htslib ; bgzip {input.infile}'

rule recode_pedigree:
    input:
        key = rules.sort_pedigree.output.key,
        pedigree = rules.sort_pedigree.output.pedigree,
    output:
        pedigree = OUT_DIR + '/REFALT/{breed}/recoded.pedigree',
    resources:
        mem_mb = 1000,
        walltime = "00:30"
    script:
        'recode_pedigree.py'


rule mendelian:
    input:
        vcf = rules.zip2.output.vcf,
        pedigree = rules.recode_pedigree.output.pedigree,
    output:
        log = OUT_DIR + '/REFALT/{breed}/CHR{chr}/mendelian.log',
        vcf = OUT_DIR + '/REFALT/{breed}/CHR{chr}/clean.vcf'
    resources:
        mem_mb = 4000,
        walltime = "02:00"
    script:
        'mendelian.py'


rule zip3:
    input:
        infile = rules.mendelian.output.vcf
    output:
        vcf = OUT_DIR + '/REFALT/{breed}/CHR{chr}/clean.vcf.gz'
    shell:
        'module load htslib ; bgzip {input.infile}'

rule recode_12:
    input:
        vcf = rules.zip3.output.vcf
    output:
        outfiles = expand(
            OUT_DIR + '/REFALT/{{breed}}/CHR{{chr}}/clean.{ext}', ext=['ped', 'map'])
    params:
        lambda wildcards, output: output.outfiles[0][:-4]
    resources:
        mem_mb = 16000,
        walltime = '00:20'
    shell:
        " {plink} --cow "
        " --vcf {input.vcf} "
        " --keep-allele-order "
        " --recode 12 "
        " --out {params[0]} "
        " --threads {threads}"

rule format_for_linkphase:
    input:
        ped = rules.recode_12.output.outfiles[0],
        map = rules.recode_12.output.outfiles[1]
    output:
        typ = OUT_DIR + \
            '/REFALT/LINKPHASE_RUN1/{breed}/CHR{chr}/linkphase.typ',
        map = OUT_DIR + '/REFALT/LINKPHASE_RUN1/{breed}/CHR{chr}/linkphase.map'
    resources:
        mem_mb = 4000,
        walltime = '00:30'
    script:
        "format_for_linkphase.py"


include: 'linkphase/linkphase.smk'
