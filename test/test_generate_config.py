import unittest

from pathlib import Path

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
        fasta_path = Path('test/data/NC_000913.3.fa')
        with self.assertRaisesRegex(ValueError, "Coverage must be above 0"):
            camiconf.get_sample_size(fasta_path, 0)
    # TODO implement success check


class TestGenerateConfig(unittest.TestCase):
    def test_only_wgsim_errorfree(self):
        camisim_dir = Path("/home/people/katste/camisim/CAMISIM")
        metadata = Path("test/data/metadata")
        id_to_genome = Path("test/data/id_to_genome_file")
        file_name = "camisim_config_TEST.ini"
        output_dir = "camisim_out"
        readsim = "art"
        readsim_dir = camisim_dir / "tools" / "art_illumina-2.3.6" / "art_illumina"
        sample = "replicates"
        amount_genomes = 2
        samplesize = 0.1
        with self.assertRaisesRegex(ValueError, "Error profile can only be omitted with wgsim"):
            camiconf.generate_config_file(camisim_dir, metadata, id_to_genome, file_name, output_dir, readsim,
                                          readsim_dir, sample, amount_genomes, samplesize)

    def test_catch_invalid_readsim(self):
        camisim_dir = Path("/home/people/katste/camisim/CAMISIM")
        metadata = Path("test/data/metadata")
        id_to_genome = Path("test/data/id_to_genome_file")
        file_name = "camisim_config_TEST.ini"
        output_dir = "camisim_out"
        readsim = "blah"
        readsim_dir = camisim_dir / "tools" / "art_illumina-2.3.6" / "art_illumina"
        sample = "replicates"
        amount_genomes = 2
        samplesize = 0.1
        error_profiles = camisim_dir / "tools" / "art_illumina-2.3.6" / "profiles"
        with self.assertRaisesRegex(ValueError, "blah is not a valid read simulator"):
            camiconf.generate_config_file(camisim_dir, metadata, id_to_genome, file_name, output_dir, readsim,
                                          readsim_dir, sample, amount_genomes, samplesize,
                                          error_profiles=error_profiles)

    def test_catch_invalid_sample(self):
        camisim_dir = Path("/home/people/katste/camisim/CAMISIM")
        metadata = Path("test/data/metadata")
        id_to_genome = Path("test/data/id_to_genome_file")
        file_name = "camisim_config_TEST.ini"
        output_dir = "camisim_out"
        readsim = "art"
        readsim_dir = camisim_dir / "tools" / "art_illumina-2.3.6" / "art_illumina"
        sample = "blah"
        amount_genomes = 2
        samplesize = 0.1
        error_profiles = camisim_dir / "tools" / "art_illumina-2.3.6" / "profiles"
        with self.assertRaisesRegex(ValueError, "blah is not a valid sample type"):
            camiconf.generate_config_file(camisim_dir, metadata, id_to_genome, file_name, output_dir, readsim,
                                          readsim_dir, sample, amount_genomes, samplesize,
                                          error_profiles=error_profiles)

    # TODO: make success testable?