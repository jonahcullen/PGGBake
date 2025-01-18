import re
import os
import csv
import glob
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.Seq import Seq
import seaborn as sns

#configfile: 'config.yaml'
workdir: config['workdir']
#singularity: config['sif']

include: 'src/utils.py'

# load assemblies into dataframe with cols sample name, breed, and fasta path
genseqs = pd.read_csv(
        Path(workflow.basedir) / 'assemblies.tsv', 
        sep='\t', 
        header=None, 
        names=['name', 'breed', 'fasta']
)

rule all:
    input:
       ## ASSEMBLIES
       #'you_did_it.done',
       ## PGGB
       #expand(
       #   #'{bucket}/public/refgen/linear/clean/comms/all_contigs.dist.tsv',
       #    u=genseqs.itertuples(),
       #    bucket=config['bucket'],
       #),

include: 'rules/prep_seqs.smk'
include: 'rules/assemblies.smk'
include: 'rules/pggb.smk'
