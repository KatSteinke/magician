# MAGICIAN
MAGICIAN is a tool for easily generating simulated metagenome-assembled genomes from a user-determined "community".
## Requirements
MAGICIAN is a Snakemake pipeline that uses conda to manage dependencies.
Thus, it primarily requires Snakemake and conda/mamba to be used.

It is also necessary to install [a fork of CAMISIM 1.2](https://github.com/KatSteinke/CAMISIM)
in order to use custom error profiles. This can be done by running
```commandline
git clone https://github.com/KatSteinke/CAMISIM
```
## Getting started
In order to get started with MAGICIAN, simply clone the repository:
```commandline
git clone https://github.com/KatSteinke/magician
```

You will also have to adapt the Snakefile. 
### Using conda (bringing your own CAMISIM)
When using your own copy of CAMISIM, set `CAMISIM_DIR` to the directory in which you installed CAMISIM.

## Running MAGICIAN
### Preparing the input 
MAGICIAN requires the following files to run:
* genome sequences of the organisms the simulated community should consist of, in genbank or fasta format
* a tab-separated file of sample distributions named `sample_distributions.tsv` in the directory where you wish to 
simulate your communities. The first column lists the paths to the genomes, the second lists sequence type 
(chromosome or plasmid) and all subsequent columns list community composition names and the relative abundance of the 
sequences in these communities:
```
| genomes          | seq_type   | community1 | community2 | ... |
|------------------|------------|------------|------------|-----|
| /path/to/genome1 | chromosome | 1          | 1.5        |     |
| /path/to/genome2 | chromosome | 1          | 0          |     |
| /path/to/plasmid | plasmid    | 1          | 1          |     |
  ```
### Starting MAGICIAN
MAGICIAN is started using `run_magician.py`:
```
run_magician.py [-h] [--profile_type {mbarc,hi,mi,hi150,own}]
                       [--profile_name PROFILE_NAME]
                       [--profile_readlength PROFILE_READLENGTH]
                       [--insert_size INSERT_SIZE] [--cluster CLUSTER]
                       target
                       --snake_flags "--cores [N_CORES] [SNAKE_FLAGS...]"

```
#### Required arguments
* `target`: the desired output file or rule. To run the entire workflow for all communities in 
`sample_distributions.tsv`, specify `all_bin_summaries` as the target here. To run the workflow for a single
community, give `summaries/bin_summary_[COMMUNITY].xlsx` here, replacing `[COMMUNITY]` with the name of the community
you wish to simulate. 
* `--snake_flags`: the flags to be passed on to Snakemake, enclosed in double quotes. As a minimum, this means `"-n "` 
for a dry run or `"--cores [N_CORES]"` (with `[N_CORES]` being the amount of cores Snakemake should use) for an actual run. \
To use conda or mamba, specify `--use-conda`
  (and `--conda-frontend conda` if required). For all else, refer to Snakemake's documentation.
#### Optional arguments
* `--profile_type`: the error profile CAMISIM should use for ART. This defaults to CAMISIM's default of `mbarc`; other choices
are `hi,mi,hi150,own` . The last allows users to specify their own profiles.
* `--profile_name`: required when specifying one's own profile. This is the base path to the forward/reverse reads' 
error profiles (e.g. `path/to/custom/profile_R` if  forward and reverse reads are located at 
`path/to/custom/profile_R1.txt` and `path/to/custom/profile_R2.txt` respectively)
* `--profile_readlength`: the read length used for the custom error profile; required when specifying one's own 
error profile.
* `--insert_size`: mean insert size for read simulation (defaults to 270 bp)
* `--cluster`: when using Snakemake's cluster mode, supply the command for submitting jobs as you would with Snakemake
