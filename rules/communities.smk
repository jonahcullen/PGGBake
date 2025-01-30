
localrules: combine_seqs
rule combine_seqs:
    input:
        expand(
            '{bucket}/public/{u.name}/{u.name}.clean.fa',
            bucket=config['bucket'],
            u=genseqs.itertuples()
        )
    output:
        '{bucket}/public/combine/all.fa.gz.gzi',
        fa  = '{bucket}/public/combine/all.fa',
        gz  = '{bucket}/public/combine/all.fa.gz',
        fai = '{bucket}/public/combine/all.fa.gz.fai',
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

#rule combine_assemblies:
#    input:
#        expand(
#            '{bucket}/public/refgen/linear/clean/{u.name}/{u.name}.clean.fa',
#            u=genseqs.itertuples(),
#            bucket=config['bucket']
#        )
#    output:
#        '{bucket}/public/refgen/linear/clean/comms/all_contigs.fa',
#        '{bucket}/public/refgen/linear/clean/comms/all_contigs.fa.gz',
#        '{bucket}/public/refgen/linear/clean/comms/all_contigs.fa.gz.fai',
#        '{bucket}/public/refgen/linear/clean/comms/all_contigs.fa.gz.gzi',
#    threads: 4
#    resources:
#        time   = 480,
#        mem_mb = 4000
#    shell:
#        '''
#            # combine all clean fastas and index
#            cat {input} > {output[0]}
#            bgzip -@ {threads} -c {output[0]} > {output[1]}
#            samtools faidx {output[1]}
#
#        '''
#
#rule sequence_dist:
#    input:
#        '{bucket}/public/refgen/linear/clean/comms/all_contigs.fa.gz.fai',
#        '{bucket}/public/refgen/linear/clean/comms/all_contigs.fa.gz.gzi',
#        fa = '{bucket}/public/refgen/linear/clean/comms/all_contigs.fa.gz',
#    output:
#        '{bucket}/public/refgen/linear/clean/comms/all_contigs.dist.tsv'
#    threads: 4
#    resources:
#        time   = 2880,
#        mem_mb = 24000
#    shell:
#        '''
#            # calculate mash distance
#            mash dist {input.fa} {input.fa} -s 10000 -i > {output}
#        '''
#
#rule mash_to_network:
#    input:
#        '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv'
#    output:
#        '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.vertices.id2name.txt',
#        '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.edges.weights.txt',
#        '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.edges.list.txt'
#    threads: 1
#    resources:
#        time   = 20,
#        mem_mb = 4000
#    shell:
#        '''
#            python3 ./src/mash2net.py -m {input}
#        '''
#
#rule identify_communities:
#    input:
#        verts = '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.vertices.id2name.txt',
#        wts   = '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.edges.weights.txt',
#        edges ='{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.edges.list.txt'
#    output:
#        '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.edges.weights.txt.communities.pdf',
#    threads: 1
#    resources:
#        time   = 20,
#        mem_mb = 4000
#    shell:
#        '''
#            python3 ./src/net2communities.py \
#                -e {input.edges} \
#                -w {input.wts} \
#                -n {input.verts} \
#                --plot
#        '''


# sequence divergence to find ideal value of of p (default 90)
#localrules: chrom_split
#checkpoint chrom_split:
#    input:
#        fa  = rules.combine_seqs.output.fa,
#        fai = rules.combine_seqs.output.fai
#    output:
#        directory('{bucket}/public/combine/chroms/')
#    threads: 1
#    resources:
#        time   = 600,
#        mem_mb = 4000
#    shell:
#        '''
#            mkdir -p {output}
#
#            cut -f 1 {input.fai} | cut -f 3 -d '#' | sort | uniq | while read CHROM; do
#                CHR_FASTA={output}/$CHROM.fa.gz
#                samtools faidx {input.fa} $(grep -P $"$CHROM\t" {input.fai} | cut -f 1) | bgzip -@ {threads} > $CHR_FASTA
#                echo "Generated $CHR_FASTA"
#            done
#        '''
