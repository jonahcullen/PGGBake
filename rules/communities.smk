
localrules: combine_seqs
rule combine_seqs:
    input:
        expand(
            '{bucket}/public/{u.name}/{u.name}.clean.fa',
            bucket=config['bucket'],
            u=genseqs.itertuples()
        )
    output:
        fa = '{bucket}/public/combine/all.fa',
        gz = '{bucket}/public/combine/all.fa.gz',
    singularity: config['pggb']['image']
    threads: 4
    resources:
        time   = 60,
        mem_mb = 4000
    shell:
        '''
            cat {input} > {output.fa}
            bgzip -@ {threads} -c {output.fa} > {output.gz}
            samtools faidx {output.gz}
        '''

rule sequence_dist:
    input:
        rules.combine_seqs.output.gz,
    output:
        '{bucket}/public/combine/all.dist.tsv'
    singularity: config['image']['mash']
    threads: 16
    resources:
        time   = 60,
        mem_mb = 8000
    shell:
        '''
            # calculate mash distance
            mash dist {input[0]} {input[0]} -p {threads} -s 10000 -i > {output}
        '''

localrules: mash_to_network
rule mash_to_network:
    input:
        rules.sequence_dist.output[0]
    output:
        verts     = '{bucket}/public/combine/all.dist.tsv.vertices.id2name.txt',
        edge_wts  = '{bucket}/public/combine/all.dist.tsv.edges.weights.txt',
        edge_list = '{bucket}/public/combine/all.dist.tsv.edges.list.txt'
    params:
        mash2net = Path(workflow.basedir) / 'src' / 'mash2net.py'
    threads: 1
    resources:
        time   = 20,
        mem_mb = 4000
    shell:
        '''
            python3 {params.mash2net} -m {input}
        '''

localrules: identify_communities
rule identify_communities:
    input:
        edge_list = rules.mash_to_network.output.edge_list,
        edge_wts  = rules.mash_to_network.output.edge_wts,
        verts     = rules.mash_to_network.output.verts,
    output:
        '{bucket}/public/combine/all.dist.tsv.edges.weights.txt.communities.pdf',
    params:
        net2comm = Path(workflow.basedir) / 'src' / 'net2communities.py'
    threads: 1
    resources:
        time   = 20,
        mem_mb = 4000
    shell:
        '''
            python3 {params.net2comm} \
                -e {input.edge_list} \
                -w {input.edge_wts} \
                -n {input.verts} \
                --plot
        '''

