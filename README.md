# PGGBake

To get started, you will need:
- `snakemake ~=7.0` (recommended `mamba` install)
- `python-igraph`
- `pycairo`

```
mamba create -n snake -c bioconda -c conda-forge snakemake=7.19 python-igraph pycairo
```

## Assembly data
Prepare _chm13.fa_ and then two haplotypes following MemPanG24 Day2 instructions
```
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG01978/assemblies/year1_f1_assembly_v2_genbank/HG01978.paternal.f1_assembly_v2_genbank.fa.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC/HG01978/assemblies/year1_f1_assembly_v2_genbank/HG01978.maternal.f1_assembly_v2_genbank.fa.gz
```
(unarchive these assemblies for example reasons)

## Metadata
Create a tsv (no header) with sample name and absolute path to FASTAs
```
HG01978.1	/path/to/HG01978.paternal.f1_assembly_v2_genbank.fa
HG01978.2	/path/to/HG01978.maternal.f1_assembly_v2_genbank.fa
CHM13	/path/to/chm13.fa
```
and modify `assems` in `configs/config.hprc.yaml` to point to this tsv. While there, modify `workdir` to be your preferred working space for pipeline outputs. Please ensure directory exists. You can also modify `bucket` as this value will prefix most if not all input/outputs (`workdir \ config['bucket']`. This is not necessary but it has some benefits for different workflow structures.

## Running
To test if snakemake can build the final target(s) (specified in the _one rule to rule them all_), a "dry run" or test can be conducted as
```
snakemake \
    -s snakefile \
    --configfile configs/config.hprc.yaml \
    --rerun-triggers mtime \
    -np
```
Just want to point out `--rerun-triggers` here as it helps to be aware of [these options](https://snakemake.readthedocs.io/en/stable/executing/cli.html)
```
--rerun-triggers
Possible choices: code, input, mtime, params, software-env

    Define what triggers the rerunning of a job. By default, all triggers are used, which guarantees that results are consistent with the workflow code and configuration. If you rather prefer the traditional way of just considering file modification dates, use ‘–rerun-trigger mtime’
```
