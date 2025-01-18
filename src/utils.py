
import re
import bz2
import gzip
import lzma

from Bio import SeqIO

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
    Get fasta file for each assemblies.
    '''
    seq = genseqs.loc[genseqs['name'] == (wildcards.name), ['fasta']].dropna()

    return {'assembly': seq.fasta}

def clean_sequence(seq):
    '''
    Clean sequence by removing invalid characters (only A, C, G, T, U, N valid).
    '''
    return Seq(''.join([c for c in seq if c in 'ACGTUNacgtun']))

def calc_lengths(f):
    '''
    Calculate sequence sequences from fasta input.
    '''
    l = []
    for rec in SeqIO.parse(f, 'fasta'):
        l.append(len(rec.seq))
    l.sort(reverse=True)
    return l

def eval_sample(wildcards):
    '''
    get fastq files for each sample
    '''
    df = map_evals.loc[(wildcards.sample), ['r1', 'r2']].dropna()

    return {'r1': df.r1, 'r2': df.r2}

def extract_rg(f, sample_name):
    with gzip.open(f, 'rb') as f:
        header = f.readline().decode('utf-8').strip()

    try:
        # attempt to match illumina header
       #ill_match = re.match(r'^@(\S+):(\d+):(\S+):(\d+):\d+:\d+:\d+ \S+:\S+:\S+:(\S+)', header)
        ill_match_full = re.match(
            r'^@(\S+):(\d+):(\S+):(\d+):\d+:\d+:\d+ \S+:\S+:\S+:(\S+)$', 
            header
        )
        if ill_match_full:
            instrument_id = ill_match_full.group(1)
            run_number = ill_match_full.group(2)
            flowcell_id = ill_match_full.group(3)
            lane = ill_match_full.group(4)
            sample_barcode = ill_match_full.group(5)
        
            return {
                'ID': f'{sample_name}.{flowcell_id}.{lane}',
                'PL': 'ILLUMINA',
                'LB': f'{sample_name}_lib',
                'SM': sample_name,
                'PU': f'{flowcell_id}.{lane}.{sample_barcode}'
            }
        # if no clean match, try simplified version
        ill_match_simp = re.match(
            r'^@(\S+):(\d+):(\S+):(\d+):\d+:\d+:\d+/\d+$',
            header
        )
        if ill_match_simp:
            instrument_id = ill_match_simp.group(1)
            run_number = ill_match_simp.group(2)
            flowcell_id = ill_match_simp.group(3)
            lane = ill_match_simp.group(4)

            return {
                'ID': f'{sample_name}.{flowcell_id}.{lane}',
                'PL': 'ILLUMINA',
                'LB': f'{sample_name}_lib',
                'SM': sample_name,
                'PU': f'{flowcell_id}.{lane}.NA'
            }

        # handle generic public header
        pub_match = re.match(r'^@(\S+)$', header)
        if pub_match:
            extract_name = pub_match.group(1)

            return {
                'ID': f'{extract_name}.NA.NA',
                'PL': 'UNKNOWN',
                'LB': f'{extract_name}_lib',
                'SM': sample_name,
                'PU': 'NA.NA.NA'
            }

        # other generic headers
        other_match = re.match(r'^@(\S+)\s+\S+\s+length=\d+', header)
        if other_match:
            extract_name = other_match.group(1)
            
            return {
                'ID': f'{extract_name}.NA.NA',
                'PL': 'UNKNOWN',
                'LB': f'{extract_name}_lib',
                'SM': sample_name,
                'PU': f'NA.NA.NA'
            }

        raise ValueError(f'header format not recognized: {header}')

    except Exception as e:
        raise ValueError(f'unexpected error while parsing header: {e}. header: {header}')

