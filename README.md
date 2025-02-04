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
Modify the `assemblies.tsv` (no header) with sample name and absolute path to FASTAs
```
HG01978.1	/path/to/HG01978.paternal.f1_assembly_v2_genbank.fa
HG01978.2	/path/to/HG01978.maternal.f1_assembly_v2_genbank.fa
CHM13	/path/to/chm13.fa
```
and if you changed the name (i.e. `assemblies.Jan25.tsv`), update `assems` in `config/config.hprc.yaml` to point to this tsv. While there, modify `workdir` to be your preferred working space for workflow outputs. Please ensure directory exists. You can also modify `bucket` as this value will prefix most if not all input/outputs (`workdir \ config['bucket']`. This is not necessary but it has some benefits for different workflow structures.

## Running
To test if snakemake can build the final target(s) (specified in the _one rule to rule them all_), a "dry run" or test can be conducted as
```
snakemake \
    --configfile config/config.hprc.yaml \
    --rerun-triggers mtime \
    -np
```
Just want to point out `--rerun-triggers` here as it helps to be aware of [these options](https://snakemake.readthedocs.io/en/stable/executing/cli.html)
```
--rerun-triggers
Possible choices: code, input, mtime, params, software-env

    Define what triggers the rerunning of a job. By default, all triggers are used, which guarantees that results are consistent with the workflow code and configuration. If you rather prefer the traditional way of just considering file modification dates, use ‘–rerun-trigger mtime’
```
If the dry run looks _good_, the workflow can be run locally like
```
ASSEM_DIR=/path/to/directory/of/assemblies/
PROFILE=/path/to/pggbake/profiles/local.pggbake

export APPTAINER_TMPDIR=/path/apptainer/tmpdir
export TMPDIR=/path/to/tmpdir

snakemake \
    --configfile config/config.hprc.yaml \
    --use-singularity \
    --singularity-args "-B $ASSEM_DIR" \
    --profile "$PROFILE" \
    --rerun-triggers mtime \
    --keep-going
```
or on an HPC with slurm (assumes the batch script within `pggbake` repo) like
```
#!/bin/bash -l
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2gb
#SBATCH --job-name submit_pggbake.slurm
#SBATCH -o %j.submit_pggbake.out
#SBATCH -e %j.submit_pggbake.err

set -e
cd $SLURM_SUBMIT_DIR
conda activate SNAKENV

ASSEM_DIR=/users/9/cull0084/projects/PGGBake/assemblies/
PROFILE="$SLURM_SUBMIT_DIR/profiles/slurm.pggbake"

export APPTAINER_TMPDIR=/scratch.global/mccuem_PGGBAKE
export TMPDIR=/scratch.global/mccuem_PGGBAKE

snakemake \
    --configfile configs/config.hprc.yaml \
    --use-singularity \
    --singularity-args "-B $PWD,$ASSEM_DIR" \
    --profile "$PROFILE" \
    --rerun-triggers mtime \
    --keep-going
```
