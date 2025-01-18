
rule combine_seqs:
    input:
        expand(
            '{bucket}/public/refgen/linear/clean/{u.breed}/{u.name}/{u.name}.clean.fa'
            u=genseqs.itertuples(),
            bucket=config['bucket']
        )
    output:
        '{bucket}/public/refgen/linear/clean/sequence/all_seqs.fa',
    threads: 4
    resources:
        time   = 60,
        mem_mb = 4000
    shell:
        '''
            cat {input} > {output}
        '''

rule sequence_dist:
    input:
        fltr_fa = '{bucket}/public/refgen/linear/clean/sequence/filter_seqs.fa',
    output:
        '{bucket}/public/refgen/linear/clean/sequence/filter_seqs.fa.gz.fai',
        '{bucket}/public/refgen/linear/clean/sequence/filter_seqs.fa.gz.gzi',
        fa   = '{bucket}/public/refgen/linear/clean/sequence/filter_seqs.fa.gz',
        dist = '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv'
    threads: 4
    resources:
        time   = 360,
        mem_mb = 8000
    shell:
        '''
            # combine all clean fastas and index
            bgzip -@ {threads} -c {input.fltr_fa} > {output.fa}
            samtools faidx {output.fa}

            # calculate mash distance
            mash dist {output.fa} {output.fa} -s 10000 -i > {output.dist}
        '''

rule mash_to_network:
    input:
        '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv'
    output:
        '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.vertices.id2name.txt',
        '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.edges.weights.txt',
        '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.edges.list.txt'
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

rule identify_communities:
    input:
        verts = '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.vertices.id2name.txt',
        wts   = '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.edges.weights.txt',
        edges ='{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.edges.list.txt'
    output:
        '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.edges.weights.txt.communities.pdf',
    params:
        net2comm = Path(workflow.basedir) / 'src' / 'net2communities.py'
    threads: 1
    resources:
        time   = 20,
        mem_mb = 4000
    shell:
        '''
            python3 {params.net2comm} \
                -e {input.edges} \
                -w {input.wts} \
                -n {input.verts} \
                --plot
        '''

