
def extract_haplotype(s):
    match = re.search(r'\.(\d+)$', s)
    if match:
        return match.group(1)
    return '0'

# NOTE: need to force users to ensure sample names do not contain
#       any special characters (match mc and pggb instructions) and
#       haplotypes specified with SAMPLE.1, SAMPLE.2, or SAMPLE
rule clean_assembly:
    input:
        unpack(get_fasta)
    output:
        clean_fa = '{bucket}/public/refgen/linear/clean/{breed}/{name}/{name}.clean.fa'
    threads: 1
    resources:
        time   = 10,
        mem_mb = 4000
    run:
        with open(output.clean_fa, 'w') as f_out:
            with RawFile(str(input.assembly[0])) as f_in:
                for rec in SeqIO.parse(f_in, 'fasta'):
                    # get haploid number if
                    hap_num = extract_haplotype(wildcards.name)
                    # clean and write to output
                    rec.id = f'{wildcards.name.rsplit(".", 1)[0]}#{hap_num}#{rec.id.strip()}'
                    rec.description = ''
                    rec.seq = clean_sequence(rec.seq)
                    SeqIO.write(rec, f_out, 'fasta')

rule index_clean_assembly:
    input:
        '{bucket}/public/refgen/linear/clean/{breed}/{name}/{name}.clean.fa'
    output:
        '{bucket}/public/refgen/linear/clean/{breed}/{name}/{name}.clean.fa.fai',
        '{bucket}/public/refgen/linear/clean/{breed}/{name}/{name}.clean.contigs.bed'
    threads: 1
    resources:
        time   = 10,
        mem_mb = 4000
    shell:
        '''
            samtools faidx {input}
            awk 'BEGIN {{OFS="\\t"}} {{print $1, 0, $2}}' {output[0]} > {output[1]}
        '''

localrules: collect_contigs
rule collect_contigs:
    input:
        expand(
            '{bucket}/public/refgen/linear/clean/{u.breed}/{u.name}/{u.name}.clean.contigs.bed',
            u=genseqs.itertuples(),
            bucket=config['bucket']
        ),
    output:
        '{bucket}/public/refgen/linear/clean/paste/all.clean.contigs.bed',
    threads: 1
    shell:
        '''
            cat {input} > {output}
        '''
