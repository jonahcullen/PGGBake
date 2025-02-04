
import re
import bz2
import gzip
import lzma

from Bio import SeqIO
from Bio.Seq import Seq

class RawFile(object):
    def __init__(self,filename):
        self.filename = filename
        if filename.endswith('.gz'):
            if self.is_gzipped(filename):
                self.handle = gzip.open(filename, 'rt')
            else:
                self.handle = open(filename,'r')
        elif filename.endswith('bz2'):
            self.handle = bz2.open(filename,'rt')
        elif filename.endswith('xz'):
            self.handle = lzma.open(filenaem,'rt')
        else:
            self.handle = open(filename,'r')

    def is_gzipped(self, filename):
        try:
            with gzip.open(filename, 'rt') as f_in:
                f_in.read(1)
            return True
        except gzip.BadGzipFile:
            return False

    def __enter__(self):
        return self.handle

    def __exit__(self,dtype,value,traceback):
        self.handle.close()

def get_fasta(wildcards):
    '''
    get fasta file for each assemblies
    '''
    seq = genseqs.loc[genseqs['name'] == (wildcards.name), ['fasta']].dropna()

    return {'assembly': seq.fasta}

def clean_sequence(seq):
    '''
    clean sequence by removing invalid characters (only A, C, G, T, U, N valid)
    '''
    return Seq(''.join([c for c in seq if c in 'ACGTUNacgtun']))

def calc_lengths(f):
    '''
    calculate sequence sequences from fasta input
    '''
    l = []
    for rec in SeqIO.parse(f, 'fasta'):
        l.append(len(rec.seq))
    l.sort(reverse=True)
    return l
