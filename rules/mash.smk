
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
        phylo_tree = Path(workflow.basedir) / 'src' / 'phylo_tree_assembly.R'
    script:
        '{params.phylo_tree}'

