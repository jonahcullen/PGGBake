
rule combine_assemblies:
    input:
        expand(
            '{bucket}/public/refgen/linear/clean/{u.breed}/{u.name}/{u.name}.clean.fa',
            u=genseqs.itertuples(),
            bucket=config['bucket']
        )
    output:
        '{bucket}/public/refgen/linear/clean/comms/all_contigs.fa',
        '{bucket}/public/refgen/linear/clean/comms/all_contigs.fa.gz',
        '{bucket}/public/refgen/linear/clean/comms/all_contigs.fa.gz.fai',
        '{bucket}/public/refgen/linear/clean/comms/all_contigs.fa.gz.gzi',
    threads: 4
    resources:
        time   = 480,
        mem_mb = 4000
    shell:
        '''
            # combine all clean fastas and index
            cat {input} > {output[0]}
            bgzip -@ {threads} -c {output[0]} > {output[1]}
            samtools faidx {output[1]}

        '''

rule sequence_dist:
    input:
        '{bucket}/public/refgen/linear/clean/comms/all_contigs.fa.gz.fai',
        '{bucket}/public/refgen/linear/clean/comms/all_contigs.fa.gz.gzi',
        fa = '{bucket}/public/refgen/linear/clean/comms/all_contigs.fa.gz',
    output:
        '{bucket}/public/refgen/linear/clean/comms/all_contigs.dist.tsv'
    threads: 4
    resources:
        time   = 2880,
        mem_mb = 24000
    shell:
        '''
            # calculate mash distance
            mash dist {input.fa} {input.fa} -s 10000 -i > {output}
        '''

rule mash_to_network:
    input:
        '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv'
    output:
        '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.vertices.id2name.txt',
        '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.edges.weights.txt',
        '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.edges.list.txt'
    threads: 1
    resources:
        time   = 20,
        mem_mb = 4000
    shell:
        '''
            python3 ./src/mash2net.py -m {input}
        '''

rule identify_communities:
    input:
        verts = '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.vertices.id2name.txt',
        wts   = '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.edges.weights.txt',
        edges ='{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.edges.list.txt'
    output:
        '{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.edges.weights.txt.communities.pdf',
       #'{bucket}/public/refgen/linear/clean/sequence/seqs.dist.tsv.edges.weights.txt.communities_filtered.pdf',
    threads: 1
    resources:
        time   = 20,
        mem_mb = 4000
    shell:
        '''
            python3 ./src/net2communities.py \
                -e {input.edges} \
                -w {input.wts} \
                -n {input.verts} \
                --plot
        '''

# sequence divergence to find ideal value of of p (default 90)
localrules: contig_split
checkpoint contig_split:
    input:
        fa  = '{bucket}/public/refgen/linear/clean/comms/all_contigs.fa.gz',
        fai = '{bucket}/public/refgen/linear/clean/comms/all_contigs.fa.gz.fai',
       #fa = '{bucket}/public/refgen/linear/clean/sequence/filter_seqs.fa.gz',
       #fai = '{bucket}/public/refgen/linear/clean/sequence/filter_seqs.fa.gz.fai',
    output:
        directory('{bucket}/public/refgen/linear/clean/comms/split')
    threads: 1
    resources:
        time   = 600,
        mem_mb = 4000
    shell:
        '''
            mkdir -p {output}

            cut -f 1 {input.fai} | cut -f 3 -d '#' | sort | uniq | while read CHROM; do
                CHR_FASTA={output}/$CHROM.fa.gz
                samtools faidx {input.fa} $(grep -P $"$CHROM\t" {input.fai} | cut -f 1) | bgzip -@ {threads} > $CHR_FASTA
                echo "Generated $CHR_FASTA"
            done
        '''

def get_chrom_fastas(wildcards):
    fas_dir = checkpoints.contig_split.get(**wildcards).output[0]
    fastas = str(Path(fas_dir) / '{chrom}.fa.gz')
    CHROMS, = glob_wildcards(fastas)
    return sorted(expand(
        '{bucket}/public/refgen/linear/clean/chroms/{chrom}.fa.gz',
        bucket=config['bucket'],
        chrom=CHROMS,
    ))

localrules: max_divergence
rule max_divergence:
    input:
        get_chrom_fastas
    output:
        '{bucket}/public/refgen/linear/clean/divergence/filter_seqs.divergence.txt',
    threads: 4
    resources:
        time   = 60,
        mem_mb = 8000
    shell:
        '''
            ls {input} | while read CHR_FASTA; do
                CHROM=$(basename $CHR_FASTA | cut -f 1 -d '.')
                MAX_DIVERGENCE=$(mash triangle -p 4 $CHR_FASTA | sed 1,1d | tr '\t' '\n' | grep chr -v | LC_ALL=C  sort -g -k 1nr | uniq | head -n 1)

                echo -e "$CHROM\t$MAX_DIVERGENCE" >> {output}
            done
        '''

# for pggb chrom (at least for now) use the value calculated from the max
# divergence
rule pggb_chrom_builds:
    input:
        '{bucket}/public/refgen/linear/clean/chroms/{chrom}.fa.gz',
    output:
       #'{bucket}/public/pggb/chroms/builds/{chrom}/{wfmash}.{seqwish}.{smoothxg}/{chrom}.fa.gz.{wfmash}.{seqwish}.{smoothxg}.smooth.final.gfa'
        gfa          = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.gfa',
        og           = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og',
        lay          = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.lay',
        lay_tsv      = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.lay.tsv',
        lay_draw     = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.lay.draw.png',
        mqc_lay_draw = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.lay.draw_multiqc.png',
        og_stats     = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.stats.yaml',
        viz_o        = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.viz_O_multiqc.png',
        viz_depth    = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.viz_depth_multiqc.png',
        viz_inv      = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.viz_inv_multiqc.png',
        viz_mqc      = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.viz_multiqc.png',
        viz_pos      = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.viz_pos_multiqc.png',
        viz_uncall   = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.viz_uncalled_multiqc.png',
        vcf          = f'{{bucket}}/public/pggb/chroms/builds/{{chrom}}/s{{segm_len}}_p{{perc_ident}}_k{{min_mat_len}}/G{{poa_len}}/{{chrom}}.s{{segm_len}}_p{{perc_ident}}_k{{min_mat_len}}.{config["pggb"]["ref_call"]}.smooth.final.vcf',
        vcf_stats    = f'{{bucket}}/public/pggb/chroms/builds/{{chrom}}/s{{segm_len}}_p{{perc_ident}}_k{{min_mat_len}}/G{{poa_len}}/{{chrom}}.s{{segm_len}}_p{{perc_ident}}_k{{min_mat_len}}.{config["pggb"]["ref_call"]}.smooth.final.vcf.stats',
        affixes      = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.fix.affixes.tsv.gz',
        seqw_stats   = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.seqwish.og.stats.yaml',
        paf          = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.alignments.wfmash.paf',
        mqc_html     = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc_report.html',
        mqc_yaml     = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc_config.yaml',
        mqc_log      = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/multiqc_data/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc.log',
        mqc_bcft     = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/multiqc_data/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc_bcftools_stats.txt',
        mqc_cits     = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/multiqc_data/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc_citations.txt',
        mqc_json     = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/multiqc_data/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc_data.json',
        mqc_stats    = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/multiqc_data/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc_general_stats.txt',
        mqc_odgi     = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/multiqc_data/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc_odgi_stats.txt',
        mqc_vers     = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/multiqc_data/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc_software_versions.txt',
        mqc_src      = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/multiqc_data/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc_sources.txt',
        params_yaml  = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.params.yaml',
        build_log    = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.build.log',
    params:
        fa_name      = lambda wildcards, input: str(Path(input[0]).name),
        outdir       = lambda wildcards, output: str(Path(output.gfa).parent),
        ref_call     = config['pggb']['ref_call'],
        n_haps       = int(config['pggb']['haps']),
        n_haps_m1    = int(config['pggb']['haps'])-1,
        poa_len      = lambda wildcards, input: str(wildcards.poa_len.replace('.', ',')),
       #block_id_min = lambda wildcards, input: "{:.4f}".format(float(wildcards.perc_ident)/100.0)
       #segm_len    = config['pggb']['segm_len'],
       #perc_ident  = config['pggb']['perc_ident'],
       #min_mat_len = config['pggb']['min_mat_len'],
    threads: 24
    resources:
        time   = 1080,
        mem_mb = 180000
    shell:
        '''
            # remove the snakemake generated multiqc directory to avoid multiqc
            # creating a multiqc_data_1 directory
            rm -rf {params.outdir}/multiqc_data/

            block_id_min=$(echo "scale=4; {wildcards.perc_ident} / 100.0" | bc)

            wfmash=$(echo W-s{wildcards.segm_len}-l25000-p{wildcards.perc_ident}-n{params.n_haps_m1}-K19-F0.001-xfalse-X | sha256sum | head -c 7)
            seqwish=$(echo k{wildcards.min_mat_len}-f0-B10000000 | sha256sum | head -c 7)
            echo GOODLUCK
            echo h{params.n_haps}-G{params.poa_len}-j0-e0-d100-I"$block_id_min"-R0-p1,19,39,3,81,1-O0.001
            smoothxg=$(echo h{params.n_haps}-G{params.poa_len}-j0-e0-d100-I"$block_id_min"-R0-p1,19,39,3,81,1-O0.001 | sha256sum | head -c 7)

            gfa_src={params.fa_name}."$wfmash"."$seqwish"."$smoothxg"
            
            ./pggb-0.5.4/pggb \
                -i {input} \
                -o {params.outdir} \
                -s {wildcards.segm_len} -p {wildcards.perc_ident} -l 25000 -n {params.n_haps} -K 19 -F 0.001 \
                -k {wildcards.min_mat_len} -f 0 -B 10000000 \
                -H {params.n_haps} -j 0 -e 0 -G {params.poa_len} -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ \
                -V {params.ref_call}:# \
                --multiqc -S \
                --threads {threads} --poa-threads {threads}

            mv {params.outdir}/"$gfa_src".smooth.final.gfa {output.gfa}
            
            mv {params.outdir}/"$gfa_src".smooth.final.og {output.og}
            mv {params.outdir}/"$gfa_src".smooth.final.og.lay {output.lay}
            mv {params.outdir}/"$gfa_src".smooth.final.og.lay.tsv {output.lay_tsv}
            mv {params.outdir}/"$gfa_src".smooth.final.og.lay.draw.png {output.lay_draw}
            mv {params.outdir}/"$gfa_src".smooth.final.og.lay.draw_multiqc.png {output.mqc_lay_draw}
            mv {params.outdir}/"$gfa_src".smooth.final.og.stats.yaml {output.og_stats}
            mv {params.outdir}/"$gfa_src".smooth.final.og.viz_O_multiqc.png {output.viz_o}
            mv {params.outdir}/"$gfa_src".smooth.final.og.viz_depth_multiqc.png {output.viz_depth}
            mv {params.outdir}/"$gfa_src".smooth.final.og.viz_inv_multiqc.png {output.viz_inv}
            mv {params.outdir}/"$gfa_src".smooth.final.og.viz_multiqc.png {output.viz_mqc}
            mv {params.outdir}/"$gfa_src".smooth.final.og.viz_pos_multiqc.png {output.viz_pos}
            mv {params.outdir}/"$gfa_src".smooth.final.og.viz_uncalled_multiqc.png {output.viz_uncall}
            
            mv {params.outdir}/"$gfa_src".smooth.final.{params.ref_call}.vcf {output.vcf}
            mv {params.outdir}/"$gfa_src".smooth.final.{params.ref_call}.vcf.stats {output.vcf_stats}
            mv {params.outdir}/"$gfa_src".smooth.fix.affixes.tsv.gz {output.affixes}
            
            mv {params.outdir}/{params.fa_name}."$wfmash"."$seqwish".seqwish.og.stats.yaml {output.seqw_stats}
            mv {params.outdir}/{params.fa_name}."$wfmash".alignments.wfmash.paf {output.paf}
            
            mv {params.outdir}/multiqc_report.html {output.mqc_html}
            mv {params.outdir}/multiqc_config.yaml {output.mqc_yaml}
            mv {params.outdir}/multiqc_data/multiqc.log {output.mqc_log}
            mv {params.outdir}/multiqc_data/multiqc_bcftools_stats.txt {output.mqc_bcft}
            mv {params.outdir}/multiqc_data/multiqc_citations.txt {output.mqc_cits}
            mv {params.outdir}/multiqc_data/multiqc_data.json {output.mqc_json}
            mv {params.outdir}/multiqc_data/multiqc_general_stats.txt {output.mqc_stats}
            mv {params.outdir}/multiqc_data/multiqc_odgi_stats.txt {output.mqc_odgi}
            mv {params.outdir}/multiqc_data/multiqc_software_versions.txt {output.mqc_vers}
            mv {params.outdir}/multiqc_data/multiqc_sources.txt {output.mqc_src}

            cat {params.outdir}/"$gfa_src".smooth.*.params.yml > {output.params_yaml}
            cat {params.outdir}/"$gfa_src".smooth.*.log > {output.build_log}
        '''
# should we adjust -G (poa-length-target) - the chicken paper used 3079,3559

def get_chrom_gfas(wildcards):
    chroms_dir = checkpoints.chrom_split.get(**wildcards).output[0]
    fas = str(Path(chroms_dir) / "{chrom}.fa.gz")
    CHROMS, = glob_wildcards(fas)
   #CHROMS = ['chr30']
   
   # return list of community gfas
    return sorted(expand(
        '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.gfa',
        bucket=config['bucket'],
        chrom=CHROMS,
        segm_len=config['pggb']['segm_len'],
        perc_ident=config['pggb']['perc_ident'],
        min_mat_len=config['pggb']['min_mat_len'],
        poa_len=config['pggb']['poa_len'],
        ref_call=config['pggb']['ref_call']
    ))

localrules: CHROMFAKE
rule CHROMFAKE:
    input:
        get_chrom_gfas
       #get_renamed_fastas
       #'{bucket}/public/pggb/communities/structure/community_structure.tsv'
    output:
        '{bucket}/okay.txt'
    shell:
        '''
            touch {output}
        '''

