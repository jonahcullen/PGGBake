
localrules: mc_seqfile
rule mc_seqfile:
    input:
        fas = expand(
            '{bucket}/public/refgen/linear/clean/{u.breed}/{u.name}/{u.name}.clean.fa',
            u=genseqs.itertuples(),
            bucket=config['bucket']
        )
    output:
        seqfile = '{bucket}/public/mc/seqfile.tsv',
        sf_done = '{bucket}/public/mc/.seqfile.done',
    run:
        with open(output.seqfile, 'w') as f_out:
            for i in input.fas:
                basename = Path(i).name.split('.clean.fa')[0]
                if basename in config['ref_assem']:
                    basename = f'Ref{basename}'
                print(basename, i, sep='\t', file=f_out)

        # create a blank file to be used by cactus_minigraph and
        # allowing seqfile to be a param
        Path(output.sf_done).touch()


           #set +u
           #source {params.mc_venv}
           #set -u
# MODIFY WORKDIR TO INCLUDE A VERSION DIRECTORY FROM THE CONFIG
# AND USE THE VERSION CONFIG THROUGHOUT IN NEXT ITERATION
rule cactus_minigraph:
    input:
       #seqfile = '{bucket}/public/mc/seqfile.tsv'
        rules.mc_seqfile.output.sf_done
    output:
        gfa = '{bucket}/public/mc/cactus_minigraph/equine.sv.gfa',
    params:
        seqfile    = rules.mc_seqfile.output.seqfile,
        jobstore   = lambda wildcards, output: Path(output.gfa).parent / 'js',
        mc_workdir = lambda wildcards, output: Path(output.gfa).parent / 'gfa_workdir',
        ref_assem  = ' '.join(['Ref' + i for i in config['ref_assem']]),
    singularity: config['mc']['image'] 
    benchmark:
        '{bucket}/public/mc/cactus_minigraph/cactus_minigraph.benchmark.txt'
    threads: 24
    resources:
        time   = 1440,
        mem_mb = 80000
    shell:
        '''
            mkdir -p {params.mc_workdir}
            
            cactus-minigraph \
                {params.jobstore} \
                {params.seqfile} \
                {output.gfa}  \
                --workDir {params.mc_workdir} \
                --reference {params.ref_assem} \
                --maxCores {threads} \
                --maxMemory {resources.mem_mb}M \
                --binariesMode local \
                --clean always
        '''
               #python3 -c "import os; os.system('toil clean {params.jobstore}')"

           ## if pythonpath not set, set to empty string to prevent 'unbound
           ## variable' error
           #if [ -z ${{PYTHONPATH+x}} ]; then export PYTHONPATH=""; fi
           #source {params.mc_venv}
           #mkdir -p {params.mc_workdir}

           #if [ -d {params.jobstore} ]; then
           #    toil clean {params.jobstore}
           #fi
rule cactus_graphmap:
    input:
        gfa = '{bucket}/public/mc/cactus_minigraph/equine.sv.gfa',
    output:
        paf     = '{bucket}/public/mc/cactus_graphmap/equine.paf',
        fa      = '{bucket}/public/mc/cactus_graphmap/equine.sv.gfa.fa',
        paf_log = '{bucket}/public/mc/cactus_graphmap/equine.paf.log',
    params:
        seqfile    = rules.mc_seqfile.output.seqfile,
        jobstore   = lambda wildcards, output: Path(output.paf).parent / 'js',
        mc_workdir = lambda wildcards, output: Path(output.paf).parent / 'paf_workdir',
        del_filter = config['mc']['del_filter'],
        ref_assem  = ' '.join(['Ref' + i for i in config['ref_assem']]),
        map_cores  = 16
    singularity: config['mc']['image'] 
    benchmark:
        '{bucket}/public/mc/cactus_graphmap/cactus_graphmap.benchmark.txt'
    threads: 24
    resources:
        time   = 1440,
        mem_mb = 120000
    shell:
        '''
            mkdir -p {params.mc_workdir}

            cactus-graphmap \
                {params.jobstore} \
                {params.seqfile} \
                {input.gfa} \
                {output.paf} \
                --outputFasta {output.fa} \
                --delFilter {params.del_filter} \
                --workDir {params.mc_workdir} \
                --reference {params.ref_assem} \
                --maxCores {threads} \
                --mapCores {params.map_cores} \
                --maxMemory {resources.mem_mb}M \
                --binariesMode local \
                --logFile {output.paf_log} \
                --clean always
        '''

rule cactus_graphmap_split:
    input:
       #seqfile = '{bucket}/public/mc/seqfile.tsv',
        gfa = '{bucket}/public/mc/cactus_minigraph/equine.sv.gfa',
        paf = '{bucket}/public/mc/cactus_graphmap/equine.paf',
    output:
        chrom_file = '{bucket}/public/mc/cactus_graphmap_split/contigs/chromfile.txt',
        split_log = '{bucket}/public/mc/cactus_graphmap_split/equine.paf.log',
    params:
        seqfile    = rules.mc_seqfile.output.seqfile,
        jobstore   = lambda wildcards, output: Path(output.split_log).parent / 'js',
        mc_workdir = lambda wildcards, output: Path(output.split_log).parent / 'split_workdir',
        contig_dir = lambda wildcards, output: Path(output.chrom_file).parent,
        ref_assem  = ' '.join(['Ref' + i for i in config['ref_assem']]),
        chr_other  = config['mc']['other_contig']
    singularity: config['mc']['image'] 
    benchmark:
        '{bucket}/public/mc/cactus_graphmap_split/cactus_graphmap_split.benchmark.txt'
    threads: 24
    resources:
        time   = 360,
        mem_mb = 60000
    shell:
        '''
            mkdir -p {params.mc_workdir}

            cactus-graphmap-split \
                {params.jobstore} \
                {params.seqfile} \
                {input.gfa} \
                {input.paf} \
                --outDir {params.contig_dir} \
                --refContigs $(for i in $(seq 31); do printf "chr$i "; done ; echo "chrX") \
                --otherContig {params.chr_other} \
                --reference {params.ref_assem} \
                --workDir {params.mc_workdir} \
                --maxCores {threads} \
                --maxMemory {resources.mem_mb}M \
                --binariesMode local \
                --logFile {output.split_log} \
                --clean always
        '''

# NOTE cactus-align has a --gpu option
# NOTE change the output to be something "better" than just a dumb log...
checkpoint cactus_align:
    input:
        chrom_file = '{bucket}/public/mc/cactus_graphmap_split/contigs/chromfile.txt',
    output:
        align_dir = directory('{bucket}/public/mc/cactus_align/vg_hal/'),
        align_log = '{bucket}/public/mc/cactus_align/equine.align.log',
    params:
        jobstore   = lambda wildcards, output: Path(output.align_log).parent / 'js',
        mc_workdir = lambda wildcards, output: Path(output.align_log).parent / 'align_workdir',
        ref_assem  = ' '.join(['Ref' + i for i in config['ref_assem']]),
        max_len    = config['mc']['max_len']
    singularity: config['mc']['image'] 
    benchmark:
        '{bucket}/public/mc/cactus_align/cactus_align.benchmark.txt'
    threads: 24
    resources:
        time   = 1440,
        mem_mb = 60000
    shell:
        '''
            mkdir -p {params.mc_workdir}

            cactus-align \
                {params.jobstore} \
                {input.chrom_file} \
                {output.align_dir} \
                --batch \
                --pangenome \
                --reference {params.ref_assem} \
                --workDir {params.mc_workdir} \
                --outVG \
                --maxLen {params.max_len} \
                --consCores {threads} \
                --maxCores {threads} \
                --maxMemory {resources.mem_mb}M \
                --binariesMode local \
                --logFile {output.align_log} \
                --clean always
        '''
           ## capture exit status and clean up if failed
           #exit_status=$?
           #echo "$exit_status"
           #echo ABOVE

           #if [ $exit_status -ne 0 ]; then
           #    toil clean {params.jobstore}
           #fi

def get_chrom_vgs(wildcards):
    outdir = checkpoints.cactus_align.get(**wildcards).output.align_dir
    vgs = str(Path(outdir) / '{chrom}.vg')
    CHROMS, = glob_wildcards(vgs)
    return sorted(expand(
        '{bucket}/public/mc/cactus_align/vg_hal/{chrom}.vg',
        bucket=config['bucket'],
        chrom=CHROMS,
    ))

def get_chrom_hals(wildcards):
    outdir = checkpoints.cactus_align.get(**wildcards).output.align_dir
    hals = str(Path(outdir) / '{chrom}.hal')
    CHROMS, = glob_wildcards(hals)
    return sorted(expand(
        '{bucket}/public/mc/cactus_align/vg_hal/{chrom}.hal',
        bucket=config['bucket'],
        chrom=CHROMS,
    ))

# NOTE - DO WE WANT TO BOTHER WITH ANY VIZ (DRAW) HERE AS AT THE CHROM LEVEL
# SHOULD WE NOT PLAN TO SORT/OPTIMIZE PRIOR TO ANY 2-D REPS
# NOTE - NEED TO INCLUDE THE CLIP GFA AND VCFWAVE VCF AS OUTPUTS
rule cactus_graphmap_join:
    input:
        vgs  = get_chrom_vgs,
        hals = get_chrom_hals
    output:
       #'{bucket}/public/mc/cactus_graphmap_join/epic-v0.2-mc.gfa.gz',
       #'{bucket}/public/mc/cactus_graphmap_join/epic-v0.2-mc.wave.vcf.gz',
        join_log = '{bucket}/public/mc/cactus_graphmap_join/equine.join.log',
    params:
        jobstore    = lambda wildcards, output: Path(output.join_log).parent / 'js',
        mc_workdir  = lambda wildcards, output: Path(output.join_log).parent / 'join_workdir',
        join_dir    = lambda wildcards, output: Path(output.join_log).parent,
        out_name    = config['mc']['out_name'],
        ref_assem   = ' '.join(['Ref' + i for i in config['ref_assem']]),
        vcf_ref     = ' '.join(config['mc']['vcf_ref']),
        clip        = config['mc']['clip'],
        freq_filter = config['mc']['filter']
    singularity: config['mc']['image'] 
    benchmark:
        '{bucket}/public/mc/cactus_graphmap_join/cactus_graphmap_join.benchmark.txt'
    threads: 32
    resources:
        time   = 2880,
        mem_mb = 160000
    shell:
        '''
            mkdir -p {params.mc_workdir}

            cactus-graphmap-join \
                {params.jobstore} \
                --vg {input.vgs} \
                --hal {input.hals} \
                --outDir {params.join_dir} \
                --outName {params.out_name} \
                --reference {params.ref_assem} \
                --clip {params.clip} \
                --filter {params.freq_filter} \
                --chop \
                --gfa {{full,clip}} \
                --odgi \
                --viz \
                --chrom-vg {{clip,filter}} \
                --chrom-og \
                --vcf \
                --vcfReference {params.vcf_ref} \
                --vcfwave \
                --gbz {{full,clip,filter}} \
                --haplo clip \
                --workDir {params.mc_workdir} \
                --indexCores $(( {threads} -1 )) \
                --indexMemory {resources.mem_mb}M \
                --maxCores {threads} \
                --maxMemory {resources.mem_mb}M \
                --binariesMode local \
                --logFile {output.join_log} \
                --stats
        '''
#### UPDATE TO V2.9.0 AND DROP --giraffe
#### DO WE NEED ODGI FULL (BIG FILE 72G - ANY IMMEDIATE VALUE?)
#### DO WE NEED TO GENERATE THE CHROM-VG AND CHROM-OG AS PART OF 
#### THIS RULE OR NEW RULE TO GEN SEPARATELY USING THE SAME COMMAND?
#### SIMILAR QUESTIONS WITH BASICALLY ALL OUTPUTS AND OUTPUT TYPES
#### remove --draw (its currently experimental and takes forever)
#### and use --chrom-og and do drawing by hand in separate rule

#### AGAIN WITH ALL REFS TO BE REF-SENSE ??

rule normalize_mc_vcf:
    input:
        vcf = '{bucket}/public/mc/cactus_graphmap_join/ecab_pg.vcf.gz',
        tbi = '{bucket}/public/mc/cactus_graphmap_join/ecab_pg.vcf.gz.tbi',
        ref = expand(
            '{bucket}/public/refgen/linear/clean/{name}/{name}.clean.fa',
            bucket=config['bucket'],
            name=config['ref_assem'][0]
        )
    output:
        '{bucket}/public/mc/norm/ecab_pg.norm.vcf.gz.tbi',
        vcf = '{bucket}/public/mc/norm/ecab_pg.norm.vcf.gz',
    threads: 12
    resources:
        time   = 1440,
        mem_mb = 30000
    shell:
        '''
            # NOTE the change from -d none to -d both
            bcftools norm --threads {threads} -m -any -f {input.ref} {input.vcf} \
                | bcftools norm --threads 12 -d both \
                | bcftools sort -Oz -o {output.vcf}
            
            tabix -p vcf {output.vcf}
        '''

rule variant_length_annotation:
    input:
        '{bucket}/public/mc/norm/ecab_pg.norm.vcf.gz.tbi',
        vcf = '{bucket}/public/mc/norm/ecab_pg.norm.vcf.gz',
    output:
        ilens_tab = '{bucket}/public/mc/norm/ilens.tab',
        ilens_gz  = '{bucket}/public/mc/norm/ilens.tab.gz',
        ilens_tbi = '{bucket}/public/mc/norm/ilens.tab.gz.tbi',
        header    = '{bucket}/public/mc/norm/ilens_info.txt'
    threads: 4
    resources:
        time   = 1440,
        mem_mb = 30000
    shell:
        '''
            bcftools query -f'%CHROM\t%POS\t%REF\t%ALT\n' {input.vcf} \
                | awk -v OFS='\t' '{{print $1, $2, $3, $4, length($4)-length($3)}}' \
                > {output.ilens_tab}

            bgzip -c {output.ilens_tab} > {output.ilens_gz}
            tabix -s1 -b2 -e2 {output.ilens_gz}

            echo '##INFO=<ID=ILEN,Number=1,Type=Integer,Description="Length difference between REF and ALT alleles">' >> {output.header}
        '''

rule length_annotate_mc_vcf:
    input:
        '{bucket}/public/mc/norm/ecab_pg.norm.vcf.gz.tbi',
        vcf       = '{bucket}/public/mc/norm/ecab_pg.norm.vcf.gz',
        ilens_gz  = '{bucket}/public/mc/norm/ilens.tab.gz',
        ilens_tbi = '{bucket}/public/mc/norm/ilens.tab.gz.tbi',
        header    = '{bucket}/public/mc/norm/ilens_info.txt'
    output:
        '{bucket}/public/mc/norm/ecab_pg.ilens.norm.vcf.gz.tbi',
        vcf = '{bucket}/public/mc/norm/ecab_pg.ilens.norm.vcf.gz',
    threads: 12
    resources:
        time   = 1440,
        mem_mb = 30000
    shell:
        '''
            bcftools annotate \
                -h {input.header} \
                -a {input.ilens_gz} \
                -c CHROM,POS,REF,ALT,ILEN \
                -Oz -o {output.vcf} \
                {input.vcf} \

            tabix -p vcf {output.vcf}
        '''

