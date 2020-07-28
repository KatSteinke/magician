import unittest

from pathlib import Path
from textwrap import dedent

import camisim_setup.generate_camisim_config as camiconf


class TestLineCount(unittest.TestCase):
    def test_file_with_lines(self):
        line_path = Path('test/data/genomefile_2line')
        assert camiconf.get_file_length(line_path) == 2

    def test_blank_file(self):
        blank_path = Path('test/data/blank_file')
        assert camiconf.get_file_length(blank_path) == 0


class TestGenerateSize(unittest.TestCase):
    def test_only_positive_coverage(self):
        fasta_path = Path('test/data/id_to_genome_file')
        with self.assertRaisesRegex(ValueError, "Coverage must be above 0"):
            camiconf.get_sample_size(fasta_path, 0)

    def test_correct_size(self):
        test_ecoli_size = 4641652
        test_bacillus_size = 4215606
        total_in_gbp_rounded = round((test_bacillus_size + test_ecoli_size)/1000000000, 2)
        fasta_path = Path('test/data/id_to_genome_file')
        assert camiconf.get_sample_size(fasta_path, 1) == total_in_gbp_rounded


class TestGenerateConfig(unittest.TestCase):
    def test_only_wgsim_errorfree(self):
        camisim_dir = Path("/home/people/katste/camisim/CAMISIM")
        metadata = Path("test/data/metadata")
        id_to_genome = Path("test/data/id_to_genome_file")
        output_dir = "camisim_out"
        readsim = "art"
        readsim_dir = camisim_dir / "tools" / "art_illumina-2.3.6" / "art_illumina"
        sample = "replicates"
        amount_genomes = 2
        samplesize = 0.1
        with self.assertRaisesRegex(ValueError, "Error profile can only be omitted with wgsim"):
            camiconf.generate_config_file(camisim_dir, metadata, id_to_genome, output_dir, readsim,
                                          readsim_dir, sample, amount_genomes, samplesize)

    def test_catch_invalid_readsim(self):
        camisim_dir = Path("/home/people/katste/camisim/CAMISIM")
        metadata = Path("test/data/metadata")
        id_to_genome = Path("test/data/id_to_genome_file")
        output_dir = "camisim_out"
        readsim = "blah"
        readsim_dir = camisim_dir / "tools" / "art_illumina-2.3.6" / "art_illumina"
        sample = "replicates"
        amount_genomes = 2
        samplesize = 0.1
        error_profiles = camisim_dir / "tools" / "art_illumina-2.3.6" / "profiles"
        with self.assertRaisesRegex(ValueError, "blah is not a valid read simulator"):
            camiconf.generate_config_file(camisim_dir, metadata, id_to_genome, output_dir, readsim,
                                          readsim_dir, sample, amount_genomes, samplesize,
                                          error_profiles=error_profiles)

    def test_catch_invalid_sample(self):
        camisim_dir = Path("/home/people/katste/camisim/CAMISIM")
        metadata = Path("test/data/metadata")
        id_to_genome = Path("test/data/id_to_genome_file")
        output_dir = "camisim_out"
        readsim = "art"
        readsim_dir = camisim_dir / "tools" / "art_illumina-2.3.6" / "art_illumina"
        sample = "blah"
        amount_genomes = 2
        samplesize = 0.1
        error_profiles = camisim_dir / "tools" / "art_illumina-2.3.6" / "profiles"
        with self.assertRaisesRegex(ValueError, "blah is not a valid sample type"):
            camiconf.generate_config_file(camisim_dir, metadata, id_to_genome, output_dir, readsim,
                                          readsim_dir, sample, amount_genomes, samplesize,
                                          error_profiles=error_profiles)

    def test_successful_config(self):
        camisim_dir = Path("/home/people/katste/camisim/CAMISIM")
        metadata = Path("test/data/metadata")
        id_to_genome = Path("test/data/id_to_genome_file")
        output_dir = "camisim_out"
        readsim = "art"
        readsim_dir = camisim_dir / "tools" / "art_illumina-2.3.6" / "art_illumina"
        sample = "replicates"
        amount_genomes = 2
        samplesize = 0.1
        error_profiles = camisim_dir / "tools" / "art_illumina-2.3.6" / "profiles"

        correct_config = '''\
        [Main]
        # maximum number of processes
        max_processors=8
        
        # 0: community design + read simulator,
        # 1: read simulator only
        phase=0
        
        # ouput directory, where the output will be stored (will be overwritten if set in from_profile)
        output_directory=camisim_out
        
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
        type=art
        
        # Samtools (http://www.htslib.org/) takes care of sam/bam files. Version 1.0 or higher required!
        # file path to executable
        samtools=/home/people/katste/camisim/CAMISIM/tools/samtools-1.3/samtools
        
        # file path to read simulation executable
        readsim=/home/people/katste/camisim/CAMISIM/tools/art_illumina-2.3.6/art_illumina
        
        #error profiles:
        #for ART:
        #HiSeq 150bp: hi150
        #MBARC-26 150bp: mbarc
        #for wgsim:
        #error rate as <float> (e.g. 0.05 for 5% error rate)
        #blank for nanosim and wgsim
        profile=mbarc
        
        # Directory containing error profiles (can be blank for wgsim)
        error_profiles=/home/people/katste/camisim/CAMISIM/tools/art_illumina-2.3.6/profiles
        
        #paired end read, insert size (not applicable for nanosim)
        fragments_size_mean=270
        fragment_size_standard_deviation=27
        
        # Only relevant if not from_profile is run:
        [CommunityDesign]
        # optional: give abundance of genomes
        distribution_file_paths=
        # specify the samples size in Giga base pairs
        size=0.1
        
        # how many different samples?
        number_of_samples=1
        
        # how many communities
        num_communities=1
        
        # directory containing the taxdump of ncbi, version from 22.02.2017 is shipped
        # "nodes.dmp"
        # "merged.dmp"
        # "names.dmp"
        ncbi_taxdump=/home/people/katste/camisim/CAMISIM/tools/ncbi-taxonomy_20170222.tar.gz
        
        # the strain simulator for de novo strain creation
        strain_simulation_template=/home/people/katste/camisim/CAMISIM/scripts/StrainSimulationWrapper/sgEvolver/simulation_dir/
        
        # define communities: [community<integer>]
        [community0]
        # information about all included genomes:
        # can be used for multiple samples
        metadata=test/data/metadata
        id_to_genome_file=test/data/id_to_genome_file
        
        # how many genomes do you want to sample over all?
        genomes_total=2
        num_real_genomes=2
        
        # how many genomes per species taxon
        #   (species taxon will be replaced by OTU-cluster later on)
        max_strains_per_otu=1
        ratio=1
        
        # which kind of different samples do you need?
        #   replicates / timeseries_lognormal / timeseries_normal / differential
        mode=replicates
        
        # Part: community design
        # Set parameters of log-normal and normal distribution, number of samples
        # sigma > 0; influences shape (higher sigma -> smaller peak and longer tail),
        log_sigma=2
        
        # mu (real number) is a parameter for the log-scale
        log_mu=1
        
        # do you want to see a distribution before you decide to use it? yes/no
        view=no
        '''
        correct_config = dedent(correct_config)
        generated_config = camiconf.generate_config_file(camisim_dir, metadata, id_to_genome, output_dir, readsim,
                                                         readsim_dir, sample, amount_genomes, samplesize,
                                                         error_profiles)

        assert generated_config == correct_config