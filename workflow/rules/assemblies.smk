
rule sketch_assembly:
    input:
        '{bucket}/public/refgen/linear/clean/{breed}/{name}/{name}.clean.fa'
    output:
        '{bucket}/public/refgen/linear/clean/{breed}/{name}/{name}.clean.msh'
    params:
        outname = lambda wildcards, input: os.path.splitext(str(input))[0]
    threads: 10
    resources:
        time   = 10,
        mem_mb = 4000
    shell:
        '''
            mash sketch -p {threads} -o {params.outname} {input}
        '''

rule paste_assemblies:
    input:
        expand(
            '{bucket}/public/refgen/linear/clean/{u.breed}/{u.name}/{u.name}.clean.msh',
            u=genseqs.itertuples(),
            bucket=config['bucket']
        )
    output:
        '{bucket}/public/refgen/linear/clean/paste/all.clean.msh'
    threads: 4
    resources:
        time   = 10,
        mem_mb = 4000
    shell:
        '''
            mash paste {output} {input}
        '''

rule quast:
    input:
        expand(
            '{bucket}/public/refgen/linear/clean/{u.breed}/{u.name}/{u.name}.clean.fa',
            u=genseqs.itertuples(),
            bucket=config['bucket']
        )
    output:
        quast = directory('{bucket}/public/refgen/linear/clean/quast_report/report.html')
    params:
        conda_env = config['conda_envs']['assem'],
    threads: 10
    resources:
        time   = 180,
        mem_mb = 30000
    shell:
        '''
            source activate {params.conda_env}
            quast.py {input} -o {output.quast}
        '''

rule assemblies_dist:
    input:
        '{bucket}/public/refgen/linear/clean/paste/all.clean.msh'
    output:
        '{bucket}/public/refgen/linear/clean/paste/all.clean.dist'
    threads: 4
    resources:
        time   = 10,
        mem_mb = 4000
    shell:
        '''
            mash dist {input} {input} > {output}
        '''

rule dist_tree:
    input:
        '{bucket}/public/refgen/linear/clean/paste/all.clean.dist'
    output:
        '{bucket}/public/refgen/linear/clean/paste/pairwise_dist.tsv',
        '{bucket}/public/refgen/linear/clean/paste/pairwise_dist.tiff',
    params:
        ref = config['ref_assem'][0]
    script:
        '../src/phylo_tree_assembly.R'

#'Rscript /opt/panimal/bin/phylo_tree_assembly.R'

# NOTE the compleasm environment
rule download_lineage:
    output:
        '{bucket}/public/lineage/{lineage}/{lineage}_odb10.done'
    params:
        conda_env = config['conda_envs']['assem'],
        outdir    = lambda wildcards, output: Path(output[0]).parent
    threads: 1
    resources:
        time   = 30,
        mem_mb = 4000
    shell:
        '''
            source activate {params.conda_env}

            compleasm download \
                {wildcards.lineage} \
                -L {params.outdir}
        '''

rule compleasm_run:
    input:
        unpack(get_fasta),
        rules.download_lineage.output
    output:
        '{bucket}/public/refgen/linear/raw/{breed}/{name}/busco/{lineage}/summary.txt',
    singularity: config['image']['panimal']
    params:
        conda_env = config['conda_envs']['assem'],
        outdir    = lambda wildcards, output: Path(output[0]).parent,
        lindir    = lambda wildcards, input: Path(input[1]).parent
    threads: 24
    resources:
        time   = 360,
        mem_mb = 60000
    shell:
        '''
            set +eu
            source activate {params.conda_env}
            set -e

            compleasm run \
                -a {input.assembly} \
                -o {params.outdir} \
                -l {wildcards.lineage} \
                -L {params.lindir} \
                -t {threads}
        '''

localrules: summarize_compleasm
rule summarize_compleasm:
    input:
        sums = expand(
            '{bucket}/public/refgen/linear/raw/{u.breed}/{u.name}/busco/{lineage}/summary.txt',
            bucket=config['bucket'],
            u=genseqs.itertuples(),
            lineage=config['lineages']
        ),
    output:
        sums_out = '{bucket}/public/refgen/linear/raw/combine/all_compleasm_summaries.csv'
    run:
        l = []
        for f in input.sums:
            with open(f, 'r') as f_in:
                lines = f_in.readlines()

            lineage = f.split('/')[-2]
            d = {}
            for line in lines[1:]:
                if ':' in line and ',' in line:
                    k,v = line.strip().split(',')
                    k = k.split(':')[0]
                    v = int(v.strip())
                    d[k] = v

            name = f.split('/')[-4]
            single_copy_complete = d.get('S', 0)
            duplicated_complete = d.get('D', 0)
            fragmented = d.get('F', 0) + d.get('I', 0)
            missing = d.get('M', 0)

            l.append([name, lineage, single_copy_complete, duplicated_complete, fragmented, missing])

        df = pd.DataFrame(
            l, 
            columns=[
                'name', 'lineage', 'single_copy_complete', 
                'duplicated_complete', 'fragmented', 'missing'
            ]
        )

        df.to_csv(output.sums_out, sep=',', index=False)

rule create_blob:
    input:
        '{bucket}/public/refgen/linear/clean/{breed}/{name}/{name}.clean.fa'
    output: 
        blob = directory('{bucket}/public/refgen/linear/clean/{breed}/{name}/blob/')
    params:
        conda_env = config['conda_envs']['assem'],
    threads: 1
    resources:
        time   = 10,
        mem_mb = 4000
    shell:
        ''' 
            source activate {params.conda_env}

            blobtools create --fasta {input} {output.blob}
        '''

rule snail_plot:
    input:
        rules.create_blob.output.blob,
    output: 
        snail      = '{bucket}/public/refgen/linear/clean/{breed}/{name}/plots/{name}.clean.snail.png',
        cumulative = '{bucket}/public/refgen/linear/clean/{breed}/{name}/plots/{name}.clean.cumulative.png'
    params:
        conda_env = config['conda_envs']['assem'],
    threads: 1
    resources:
        time   = 10,
        mem_mb = 4000
    shell:
        '''
            source activate {params.conda_env}
            
            blobtk plot -v snail -d {input} -o {output.snail}
            blobtk plot -v cumulative -d {input} -o {output.cumulative}
        '''

rule bbstats:
    input:
       #'{bucket}/public/refgen/linear/clean/{name}/{name}.clean.fa'
        '{bucket}/public/refgen/linear/clean/{breed}/{name}/{name}.clean.fa'
    output:
        stats  = '{bucket}/public/refgen/linear/clean/{breed}/{name}/bbstats/{name}.stats.tsv',
        gchist = '{bucket}/public/refgen/linear/clean/{breed}/{name}/bbstats/{name}.gchist.tsv',
        shist  = '{bucket}/public/refgen/linear/clean/{breed}/{name}/bbstats/{name}.shist.tsv',
    params:
        conda_env = config['conda_envs']['assem'],
    threads: 1
    resources:
        time   = 10,
        mem_mb = 4000
    shell:
        '''
            source activate {params.conda_env}
            
            stats.sh \
                in={input} \
                out={output.stats} \
                gchist={output.gchist} \
                shist={output.shist} \
                format=3
        '''

localrules: assembly_table
rule assembly_table:
    input:
        stats = expand(
            '{bucket}/public/refgen/linear/clean/{u.breed}/{u.name}/bbstats/{u.name}.stats.tsv',
            bucket=config['bucket'],
            u=genseqs.itertuples(),
        )
    output:
        table = '{bucket}/public/refgen/linear/clean/paste/assembly_stats.csv',
    run:
        # bbmap's stats.sh has the N and L swapped so changed here
        cols = [
            'n_scaffolds', 'n_contigs', 'scaf_bp', 'contig_bp', 
            'gap_pct', 'scaf_L50', 'scaf_N50', 'ctg_L50', 'ctg_N50', 
            'scaf_L90', 'scaf_N90', 'ctg_L90', 'ctg_N90', 'scaf_max', 
            'ctg_max', 'scaf_n_gt50K', 'scaf_pct_gt50K', 'gc_avg', 'gc_std'
        ]

        dfs = []
        for f in input.stats:
            df = pd.read_csv(f, sep='\t', header=None, skiprows=1, names=cols)
            df['assembly'] = Path(f).name.split('.stats.tsv')[0]
            dfs.append(df)

        # concat to df
        df = pd.concat(dfs, ignore_index=True)
        df.to_csv(output.table, sep=',', index=False)

rule nx_tables:
    input:
       #clean = '{bucket}/public/refgen/linear/clean/{name}/{name}.clean.fa'
        clean = '{bucket}/public/refgen/linear/clean/{breed}/{name}/{name}.clean.fa'
    output:
        lengths  = '{bucket}/public/refgen/linear/clean/{breed}/{name}/{name}.lens.csv',
    threads: 1
    resources:
        time   = 20,
        mem_mb = 4000
    run:
        assembly = Path(input.clean).name.split('.clean.fa')[0]
        seq_lens = calc_lengths(input.clean)
        c_lens = np.cumsum(seq_lens)
        t_len = c_lens[-1]
        # write to output
        with open(output.lengths, 'w', newline='') as f_out:
            writer = csv.writer(f_out)
            writer.writerow(['assembly', 'n_x', 'scaf_len'])
            for perc, length in zip(c_lens/t_len*100, seq_lens):
                writer.writerow([assembly, perc, length])

localrules: combine_nx_tables
rule combine_nx_tables:
    input:
        expand(
            '{bucket}/public/refgen/linear/clean/{u.breed}/{u.name}/{u.name}.lens.csv',
            bucket=config['bucket'],
            u=genseqs.itertuples(),
        )
    output:
        '{bucket}/public/refgen/linear/clean/paste/nx_tables.csv'
    shell:
        '''
            awk 'NR == 1 || FNR > 1' {input} > {output}
        '''

localrules: all_done_assemblies
rule all_done_assemblies:
    input:
        expand(
            '{bucket}/public/refgen/linear/clean/paste/pairwise_dist.tiff',
            bucket=config['bucket'],
        ),
        expand(
            '{bucket}/public/refgen/linear/clean/{u.breed}/{u.name}/plots/{u.name}.clean.snail.png',
            bucket=config['bucket'],
            u=genseqs.itertuples(),
        ),
        expand(
            '{bucket}/public/refgen/linear/raw/combine/all_compleasm_summaries.csv',
            bucket=config['bucket'],
        ),
        expand(
            '{bucket}/public/refgen/linear/clean/paste/assembly_stats.csv',
            bucket=config['bucket'],
        ),
        expand(
            '{bucket}/public/refgen/linear/clean/paste/nx_tables.csv',
            bucket=config['bucket'],
        ),
    output:
        "you_did_it.done"
    shell:
        '''
            touch {output}
        '''
