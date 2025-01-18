#!/usr/bin/env python3

import pandas as pd
import argparse

def merge_kmer_files(shared_file, bwauniq_file, giruniq_file, output_file, filtered_output_file, min_total=10):
    # read kmer counts into dfs
    shared = pd.read_csv(shared_file, sep='\t', names=['kmer', 'shared'], dtype={'kmer': str, 'shared': int}) \
        .groupby('kmer')['shared'].sum().reset_index()
    bwauniq = pd.read_csv(bwauniq_file, sep='\t', names=['kmer', 'bwauniq'], dtype={'kmer': str, 'bwauniq': int}) \
        .groupby('kmer')['bwauniq'].sum().reset_index()
    giruniq = pd.read_csv(giruniq_file, sep='\t', names=['kmer', 'giruniq'], dtype={'kmer': str, 'giruniq': int}) \
        .groupby('kmer')['giruniq'].sum().reset_index()

    # merge all three
    merged = shared.merge(bwauniq, on='kmer', how='outer').merge(giruniq, on='kmer', how='outer')
    merged.fillna(0, inplace=True)  # Replace NaN with 0

    # counts to integers
    merged['shared'] = merged['shared'].astype(int)
    merged['bwauniq'] = merged['bwauniq'].astype(int)
    merged['giruniq'] = merged['giruniq'].astype(int)

    # add a totals column
    merged['total'] = merged['shared'] + merged['bwauniq'] + merged['giruniq']

    # sort
    merged = merged.sort_values('total', ascending=False)

    # write to file
    merged.to_csv(output_file, sep='\t', index=False)
    print(f'full merged kmer file saved to: {output_file}')

    # filter kmers with total less than min_total
    filtered = merged[merged['total'] >= min_total]
    filtered.to_csv(filtered_output_file, sep='\t', index=False)
    print(f'filtered kmer file saved to: {filtered_output_file} (k-mers with total >= {min_total})')

def main():
    parser = argparse.ArgumentParser(description='merge and filter k-mer datasets.')
    parser.add_argument('-s', '--shared', required=True, help='path to shared k-mers file.')
    parser.add_argument('-b', '--bwauniq', required=True, help='path to BWA-unique k-mers file.')
    parser.add_argument('-g', '--giruniq', required=True, help='path to Giraffe-unique k-mers file.')
    parser.add_argument('-o', '--output', default='merged_kmers.txt', help='output file for merged k-mers (default: merged_kmers.txt).')
    parser.add_argument('-f', '--filtered', default='filtered_kmers.txt', help='output file for filtered k-mers (default: filtered_kmers.txt).')
    parser.add_argument('-t', '--min-total', type=int, default=10, help='minimum total count across categories to retain k-mers (default: 10).')

    args = parser.parse_args()
    merge_kmer_files(args.shared, args.bwauniq, args.giruniq, args.output, args.filtered, args.min_total)

if __name__ == '__main__':
    main()

