from argparse import ArgumentParser
from pathlib import Path
from textwrap import dedent
from typing import Optional

from Bio import SeqIO


# utility function for getting amount of genomes
def get_file_length(infile: Path) -> int:
    """Count the number of lines in a file.
    Arguments:
        infile: Path to the file to be counted
    Returns:
        The amount of lines in a file
    """
    count = 0
    with open(infile, "r") as file_read:
        for count, line in enumerate(file_read, 1):
            pass
    return count


def get_sample_size(file_record: Path, coverage: Optional[float]=1) -> float:
    """Parse a id_to_genome file, sum up sizes of all fasta files given and,
    if desired, multiply with an average coverage factor.
    Arguments:
        file_record:    Path to id_to_genome file
        coverage:       desired average coverage for the sample
    Returns:
        The sample size for obtaining a given average coverage level, given
        the total size of all genomes examined.
    """
    # sanity check: is the value for coverage above 0?
    if coverage <= 0:
        raise ValueError("Coverage must be above 0")
    # get all paths
    with open(file_record, "r") as idfile:
        fasta_paths = [Path(line.strip().split("\t")[1]) for line in idfile]
    # for each file, get genome size and add up
    total_size = sum([len(record) for fasta_file in fasta_paths
                      for record in SeqIO.parse(fasta_file, "fasta")])
    return round(coverage*(total_size/1000000000), 2)


def generate_config_file(camisim_dir: Path, meta_file: Path, id_file: Path, output_dir: str,
                         readsim: str, readsim_path: Path, sample: str,
                         amount_genomes: int, sample_size: float, error_profiles: Optional[Path] = "",
                         abundance_file: Optional[Path]="") -> str:
    """Generate a config file for CAMISIM and write it to a specified filename.
    Arguments:
        camisim_dir:    Path to the directory containing CAMISIM
        meta_file:      Path to the metadata file for CAMISIM
        id_file:        Path to the file linking genome IDs to fasta files
        output_dir:     output directory for CAMISIM
        readsim:        read simulator to use
        readsim_path:   Path to the read simulator to use
        sample:         type of sample to use
        amount_genomes: amount of genomes in the sample
        sample_size:    sample size in Gbp
        error_profiles: Path to error profiles, can be blank if using wgsim
        abundance_file: Path to file listing abundance for genomes, if given
    Returns:
        A CAMISIM config file with the chosen parameters.
    """
    # sanity checks
    # is the read simulator a valid choice?
    if not readsim in {"art", "wgsim", "nanosim", "pbsim"}:
        raise ValueError("{} is not a valid read simulator. Valid options are art, wgsim, nanosim, pbsim.".format(readsim))
    # is the read simulator not wgsim, but there's no error profile?
    if not error_profiles and not readsim == "wgsim":
        raise ValueError("Error profile can only be omitted with wgsim")
    # is the sample type a valid choice?
    if not sample in {'replicates', 'timeseries_lognormal', 'timeseries_normal', 'differential'}:
        raise ValueError(
            """{} is not a valid sample type. Valid options are 'replicates', 'timeseries_lognormal',\
             'timeseries_normal', 'differential'.""".format(sample))
    # TODO: introduce more defaults?
    config_string = '''\
    [Main]
    # maximum number of processes
    max_processors=8
    
    # 0: community design + read simulator,
    # 1: read simulator only
    phase=0
    
    # ouput directory, where the output will be stored (will be overwritten if set in from_profile)
    output_directory={outdir}
    
    # temporary directory
    temp_directory=/tmp
    
    # gold standard assembly
    gsa=True
    
    # gold standard for all samples combined
    pooled_gsa=True
    
    # anonymize sequences?
    anonymous=False
    
    # compress data (levels 0-9, recommended is 1 the gain of higher levels is not too high)
    compress=1
    
    # id of dataset, used in foldernames and is prefix in anonymous sequences
    dataset_id=RL
    
    # Read Simulation settings, relevant also for from_profile
    [ReadSimulator]
    # which readsimulator to use:
    #           Choice of 'art', 'wgsim', 'nanosim', 'pbsim'
    type={readsim}
    
    # Samtools (http://www.htslib.org/) takes care of sam/bam files. Version 1.0 or higher required!
    # file path to executable
    samtools={camidir}/tools/samtools-1.3/samtools
    
    # file path to read simulation executable
    readsim={simpath}
    
    #error profiles:
    #for ART:
    #HiSeq 150bp: hi150
    #MBARC-26 150bp: mbarc
    #for wgsim:
    #error rate as <float> (e.g. 0.05 for 5% error rate)
    #blank for nanosim and wgsim
    profile=mbarc
    
    # Directory containing error profiles (can be blank for wgsim)
    error_profiles={profile_path}
    
    #paired end read, insert size (not applicable for nanosim)
    fragments_size_mean=270
    fragment_size_standard_deviation=27
    
    # Only relevant if not from_profile is run:
    [CommunityDesign]
    # optional: give abundance of genomes
    distribution_file_paths={dist_file}
    # specify the samples size in Giga base pairs
    size={samplesize}
    
    # how many different samples?
    number_of_samples=1
    
    # how many communities
    num_communities=1
    
    # directory containing the taxdump of ncbi, version from 22.02.2017 is shipped
    # "nodes.dmp"
    # "merged.dmp"
    # "names.dmp"
    ncbi_taxdump={camidir}/tools/ncbi-taxonomy_20170222.tar.gz
    
    # the strain simulator for de novo strain creation
    strain_simulation_template={camidir}/scripts/StrainSimulationWrapper/sgEvolver/simulation_dir/
    
    # define communities: [community<integer>]
    [community0]
    # information about all included genomes:
    # can be used for multiple samples
    metadata={metafile}
    id_to_genome_file={genomefile}
    
    # how many genomes do you want to sample over all?
    genomes_total={amount}
    num_real_genomes={amount}
    
    # how many genomes per species taxon
    #   (species taxon will be replaced by OTU-cluster later on)
    max_strains_per_otu=1
    ratio=1
    
    # which kind of different samples do you need?
    #   replicates / timeseries_lognormal / timeseries_normal / differential
    mode={sampletype}
    
    # Part: community design
    # Set parameters of log-normal and normal distribution, number of samples
    # sigma > 0; influences shape (higher sigma -> smaller peak and longer tail),
    log_sigma=2
    
    # mu (real number) is a parameter for the log-scale
    log_mu=1
    
    # do you want to see a distribution before you decide to use it? yes/no
    view=no
    '''.format(
        camidir=camisim_dir,
        metafile=meta_file,
        genomefile=id_file,
        outdir=output_dir,
        readsim=readsim,
        simpath=readsim_path,
        sampletype=sample,
        profile_path=error_profiles,
        amount=amount_genomes,
        samplesize=sample_size,
        dist_file=abundance_file
    )
    config_string = dedent(config_string)
    return config_string


if __name__ == "__main__":
    parser = ArgumentParser(description="Generate a config file for CAMISIM from parameters")
    parser.add_argument("camisim_dir", help="Directory containing CAMISIM")
    parser.add_argument("metadata", help="Path to CAMISIM metadata file")
    parser.add_argument("genome_file", help="Path to file containing genome ID to fasta path data")
    parser.add_argument("-f", "--filename", help="Name for the config file (default: camisim_config.ini)",
                        default="camisim_config.ini")
    parser.add_argument('-o', '--out_dir', action="store", help="Output directory for CAMISIM (default: camisim_out)",
                        default="camisim_out")
    parser.add_argument('-a', '--abundance_file', action="store",
                        help="Optional: file giving relative abundance of genomes", default="")
    parser.add_argument('-c', '--coverage', action="store",
                        help="Desired average coverage for the sample (default: 1X)",
                        default=1)
    parser.add_argument('--read_sim', action="store", help="Read simulator to use",
                        choices=["art", "wgsim", "nanosim", "pbsim"],
                        default="art")
    parser.add_argument('--read_sim_path', action="store",
                        help="Path to read simulator executable (default: ART shipped with CAMISIM)",
                        default=False)
    parser.add_argument('--sample_type', action="store",
                        help="Type of different samples to be simulated (default: replicates)",
                        choices=["replicates", "timeseries_lognormal", "timeseries_normal", "differential"],
                        default="replicates")
    parser.add_argument("--error_profile", action="store",
                        help="Path to error profiles (overridden by --errorfree). Default: ART error profiles",
                        default=False)
    parser.add_argument('--errorfree', action="store_true", help="Don't use an error profile (only works with wgsim)")
    args = parser.parse_args()

    # establish location of CAMISIM dir
    camisim_dir = Path(args.camisim_dir).resolve()

    # check if no errors is specified without using wgsim - fail early if so
    if args.errorfree:
        if not args.read_sim == "wgsim":
            parser.error("Omitting an error profile is only possible with wgsim.")
        else:
            error_profile = ""
    else:
        # was the path left blank? Then use CAMISIM default
        if not args.error_profile:
            error_profile = camisim_dir / "tools" / "art_illumina-2.3.6" / "profiles"
        else:
            error_profile = Path(args.error_profile).resolve()

    # establish remaining paths
    metadata = Path(args.metadata).resolve()
    genome_file = Path(args.genome_file).resolve()
    if args.abundance_file:
        abundance_file = Path(args.abundance_file).resolve()
    else:
        abundance_file = args.abundance_file
    # if no read_sim_path, it's the absolute path to CAMISIM's ART
    if not args.read_sim_path:
        read_sim_path = camisim_dir / "tools" / "art_illumina-2.3.6" / "art_illumina"
    else:
        read_sim_path = Path(args.read_sim_path).resolve()

    # set remaining parameters
    filename = args.filename
    out_dir = args.out_dir
    read_sim = args.read_sim
    sample_type = args.sample_type
    coverage = float(args.coverage)

    #calculate amount of genomes and sample size
    genomes = get_file_length(genome_file)
    size_total = get_sample_size(genome_file, coverage)
    config_str = generate_config_file(camisim_dir, metadata, genome_file, out_dir, read_sim, read_sim_path, sample_type,
                                      genomes, size_total, error_profile, abundance_file)
    with open(filename, "w") as outfile:
        outfile.write(config_str)