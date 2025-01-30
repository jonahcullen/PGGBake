
checkpoint par_before_pggb:
    input:
        rules.combine_seqs.output.gz
    output:
        directory('{bucket}/public/combine/communities/')
    params:
        n_haps     = config['pggb']['haps'],
        par_before = Path(workflow.basedir) / 'src' / 'partition-before-pggb'
    singularity: config['pggb']['image']
    threads: 16
    resources:
        time   = 360,
        mem_mb = 120000
    shell:
        '''
            which() {{
                command -v "$1" 2>/dev/null
            }}
            export -f which

            bash {params.par_before} \
                -i {input} \
                -o {output} \
                -n {params.n_haps}
        '''

def get_community_fastas(wildcards):
    comm_dir = checkpoints.par_before_pggb.get(**wildcards).output[0]
    fastas = str(Path(comm_dir) / 'all.fa.{hash}.community.{comm}.fa')
    VHASH,COMMS, = glob_wildcards(fastas)
    return sorted(expand(
        '{bucket}/public/combine/communities/all.fa.{vhash}.community.{comm}.fa',
        bucket=config['bucket'],
        vhash=VHASH[0],
        chrom=CHROMS,
    ))

localrules: max_divergence
rule max_divergence:
    input:
        get_community_fastas
    output:
        '{bucket}/public/combine/chrom.max_divergence.txt',
    threads: 1
    resources:
        time   = 10,
        mem_mb = 100
    shell:
        '''
            touch {output}
        '''

# for pggb chrom (at least for now) use the value calculated from the max
# divergence
#rule pggb_chrom_builds:
#    input:
#        '{bucket}/public/refgen/linear/clean/chroms/{chrom}.fa.gz',
#    output:
#       #'{bucket}/public/pggb/chroms/builds/{chrom}/{wfmash}.{seqwish}.{smoothxg}/{chrom}.fa.gz.{wfmash}.{seqwish}.{smoothxg}.smooth.final.gfa'
#        gfa          = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.gfa',
#        og           = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og',
#        lay          = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.lay',
#        lay_tsv      = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.lay.tsv',
#        lay_draw     = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.lay.draw.png',
#        mqc_lay_draw = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.lay.draw_multiqc.png',
#        og_stats     = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.stats.yaml',
#        viz_o        = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.viz_O_multiqc.png',
#        viz_depth    = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.viz_depth_multiqc.png',
#        viz_inv      = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.viz_inv_multiqc.png',
#        viz_mqc      = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.viz_multiqc.png',
#        viz_pos      = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.viz_pos_multiqc.png',
#        viz_uncall   = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.og.viz_uncalled_multiqc.png',
#        vcf          = f'{{bucket}}/public/pggb/chroms/builds/{{chrom}}/s{{segm_len}}_p{{perc_ident}}_k{{min_mat_len}}/G{{poa_len}}/{{chrom}}.s{{segm_len}}_p{{perc_ident}}_k{{min_mat_len}}.{config["pggb"]["ref_call"]}.smooth.final.vcf',
#        vcf_stats    = f'{{bucket}}/public/pggb/chroms/builds/{{chrom}}/s{{segm_len}}_p{{perc_ident}}_k{{min_mat_len}}/G{{poa_len}}/{{chrom}}.s{{segm_len}}_p{{perc_ident}}_k{{min_mat_len}}.{config["pggb"]["ref_call"]}.smooth.final.vcf.stats',
#        affixes      = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.fix.affixes.tsv.gz',
#        seqw_stats   = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.seqwish.og.stats.yaml',
#        paf          = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.alignments.wfmash.paf',
#        mqc_html     = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc_report.html',
#        mqc_yaml     = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc_config.yaml',
#        mqc_log      = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/multiqc_data/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc.log',
#        mqc_bcft     = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/multiqc_data/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc_bcftools_stats.txt',
#        mqc_cits     = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/multiqc_data/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc_citations.txt',
#        mqc_json     = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/multiqc_data/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc_data.json',
#        mqc_stats    = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/multiqc_data/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc_general_stats.txt',
#        mqc_odgi     = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/multiqc_data/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc_odgi_stats.txt',
#        mqc_vers     = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/multiqc_data/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc_software_versions.txt',
#        mqc_src      = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/multiqc_data/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.multiqc_sources.txt',
#        params_yaml  = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.params.yaml',
#        build_log    = '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.build.log',
#    params:
#        fa_name      = lambda wildcards, input: str(Path(input[0]).name),
#        outdir       = lambda wildcards, output: str(Path(output.gfa).parent),
#        ref_call     = config['pggb']['ref_call'],
#        n_haps       = int(config['pggb']['haps']),
#        n_haps_m1    = int(config['pggb']['haps'])-1,
#        poa_len      = lambda wildcards, input: str(wildcards.poa_len.replace('.', ',')),
#       #block_id_min = lambda wildcards, input: "{:.4f}".format(float(wildcards.perc_ident)/100.0)
#       #segm_len    = config['pggb']['segm_len'],
#       #perc_ident  = config['pggb']['perc_ident'],
#       #min_mat_len = config['pggb']['min_mat_len'],
#    threads: 24
#    resources:
#        time   = 1080,
#        mem_mb = 180000
#    shell:
#        '''
#            # remove the snakemake generated multiqc directory to avoid multiqc
#            # creating a multiqc_data_1 directory
#            rm -rf {params.outdir}/multiqc_data/
#
#            block_id_min=$(echo "scale=4; {wildcards.perc_ident} / 100.0" | bc)
#
#            wfmash=$(echo W-s{wildcards.segm_len}-l25000-p{wildcards.perc_ident}-n{params.n_haps_m1}-K19-F0.001-xfalse-X | sha256sum | head -c 7)
#            seqwish=$(echo k{wildcards.min_mat_len}-f0-B10000000 | sha256sum | head -c 7)
#            echo GOODLUCK
#            echo h{params.n_haps}-G{params.poa_len}-j0-e0-d100-I"$block_id_min"-R0-p1,19,39,3,81,1-O0.001
#            smoothxg=$(echo h{params.n_haps}-G{params.poa_len}-j0-e0-d100-I"$block_id_min"-R0-p1,19,39,3,81,1-O0.001 | sha256sum | head -c 7)
#
#            gfa_src={params.fa_name}."$wfmash"."$seqwish"."$smoothxg"
#            
#            ./pggb-0.5.4/pggb \
#                -i {input} \
#                -o {params.outdir} \
#                -s {wildcards.segm_len} -p {wildcards.perc_ident} -l 25000 -n {params.n_haps} -K 19 -F 0.001 \
#                -k {wildcards.min_mat_len} -f 0 -B 10000000 \
#                -H {params.n_haps} -j 0 -e 0 -G {params.poa_len} -P 1,19,39,3,81,1 -O 0.001 -d 100 -Q Consensus_ \
#                -V {params.ref_call}:# \
#                --multiqc -S \
#                --threads {threads} --poa-threads {threads}
#
#            mv {params.outdir}/"$gfa_src".smooth.final.gfa {output.gfa}
#            
#            mv {params.outdir}/"$gfa_src".smooth.final.og {output.og}
#            mv {params.outdir}/"$gfa_src".smooth.final.og.lay {output.lay}
#            mv {params.outdir}/"$gfa_src".smooth.final.og.lay.tsv {output.lay_tsv}
#            mv {params.outdir}/"$gfa_src".smooth.final.og.lay.draw.png {output.lay_draw}
#            mv {params.outdir}/"$gfa_src".smooth.final.og.lay.draw_multiqc.png {output.mqc_lay_draw}
#            mv {params.outdir}/"$gfa_src".smooth.final.og.stats.yaml {output.og_stats}
#            mv {params.outdir}/"$gfa_src".smooth.final.og.viz_O_multiqc.png {output.viz_o}
#            mv {params.outdir}/"$gfa_src".smooth.final.og.viz_depth_multiqc.png {output.viz_depth}
#            mv {params.outdir}/"$gfa_src".smooth.final.og.viz_inv_multiqc.png {output.viz_inv}
#            mv {params.outdir}/"$gfa_src".smooth.final.og.viz_multiqc.png {output.viz_mqc}
#            mv {params.outdir}/"$gfa_src".smooth.final.og.viz_pos_multiqc.png {output.viz_pos}
#            mv {params.outdir}/"$gfa_src".smooth.final.og.viz_uncalled_multiqc.png {output.viz_uncall}
#            
#            mv {params.outdir}/"$gfa_src".smooth.final.{params.ref_call}.vcf {output.vcf}
#            mv {params.outdir}/"$gfa_src".smooth.final.{params.ref_call}.vcf.stats {output.vcf_stats}
#            mv {params.outdir}/"$gfa_src".smooth.fix.affixes.tsv.gz {output.affixes}
#            
#            mv {params.outdir}/{params.fa_name}."$wfmash"."$seqwish".seqwish.og.stats.yaml {output.seqw_stats}
#            mv {params.outdir}/{params.fa_name}."$wfmash".alignments.wfmash.paf {output.paf}
#            
#            mv {params.outdir}/multiqc_report.html {output.mqc_html}
#            mv {params.outdir}/multiqc_config.yaml {output.mqc_yaml}
#            mv {params.outdir}/multiqc_data/multiqc.log {output.mqc_log}
#            mv {params.outdir}/multiqc_data/multiqc_bcftools_stats.txt {output.mqc_bcft}
#            mv {params.outdir}/multiqc_data/multiqc_citations.txt {output.mqc_cits}
#            mv {params.outdir}/multiqc_data/multiqc_data.json {output.mqc_json}
#            mv {params.outdir}/multiqc_data/multiqc_general_stats.txt {output.mqc_stats}
#            mv {params.outdir}/multiqc_data/multiqc_odgi_stats.txt {output.mqc_odgi}
#            mv {params.outdir}/multiqc_data/multiqc_software_versions.txt {output.mqc_vers}
#            mv {params.outdir}/multiqc_data/multiqc_sources.txt {output.mqc_src}
#
#            cat {params.outdir}/"$gfa_src".smooth.*.params.yml > {output.params_yaml}
#            cat {params.outdir}/"$gfa_src".smooth.*.log > {output.build_log}
#        '''
## should we adjust -G (poa-length-target) - the chicken paper used 3079,3559
#
#def get_chrom_gfas(wildcards):
#    chroms_dir = checkpoints.chrom_split.get(**wildcards).output[0]
#    fas = str(Path(chroms_dir) / "{chrom}.fa.gz")
#    CHROMS, = glob_wildcards(fas)
#   #CHROMS = ['chr30']
#   
#   # return list of community gfas
#    return sorted(expand(
#        '{bucket}/public/pggb/chroms/builds/{chrom}/s{segm_len}_p{perc_ident}_k{min_mat_len}/G{poa_len}/{chrom}.s{segm_len}_p{perc_ident}_k{min_mat_len}.smooth.final.gfa',
#        bucket=config['bucket'],
#        chrom=CHROMS,
#        segm_len=config['pggb']['segm_len'],
#        perc_ident=config['pggb']['perc_ident'],
#        min_mat_len=config['pggb']['min_mat_len'],
#        poa_len=config['pggb']['poa_len'],
#        ref_call=config['pggb']['ref_call']
#    ))
#
#localrules: CHROMFAKE
#rule CHROMFAKE:
#    input:
#        get_chrom_gfas
#       #get_renamed_fastas
#       #'{bucket}/public/pggb/communities/structure/community_structure.tsv'
#    output:
#        '{bucket}/okay.txt'
#    shell:
#        '''
#            touch {output}
#        '''
#
