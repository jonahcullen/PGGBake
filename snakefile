import pandas as pd

workdir: config['workdir']
#singularity: config['sif']
#configfile: 'config.yaml'

include: 'src/utils.py'

# load assemblies into dataframe 
genseqs = pd.read_csv(
        Path(workflow.basedir) / config['assems'],
        sep='\t', 
        header=None, 
        names=['name', 'fasta']
)

rule all:
    input:
       ## leaving this here as an example
       #expand(
       #    '{bucket}/public/{u.name}/{u.name}.clean.fa',
       #    bucket=config['bucket'],
       #    u=genseqs.itertuples()
       #),
        expand(
            '{bucket}/public/combine/all.dist.tsv.edges.weights.txt.communities.pdf',
            bucket=config['bucket']
        ),

include: 'rules/prep_seqs.smk'
include: 'rules/communities.smk'
#include: 'rules/pggb.smk'
