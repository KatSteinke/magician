import unittest

import numpy as np
import pandas as pd

from pathlib import Path

import generate_summary.extract_stats as summary_stats


class TestBBStats(unittest.TestCase):
    def test_get_stats(self):
        true_stats = pd.DataFrame({"bin_name": ["bin_1", "bin_2"], "scaf_bp": [8667507, 8277261], "gc_avg": [0.6, 0.5],
                                "n_scaffolds": [4, 5], "n_contigs": [5, 6], "scaffold_L50": [2, 3],
                                "scaffold_N50": [1667507, 1277261]})
        stats_file = Path('test/data/fake_bb.csv')
        test_stats = summary_stats.get_bb_stats(stats_file)
        assert test_stats.equals(true_stats)


class TestCheckMStats(unittest.TestCase):
    def test_get_checkm(self):
        true_checkm = pd.DataFrame({"Bin Id": ["sample13.bin.1", "sample13.bin.2"],
                                "Marker lineage": ["f__Enterobacteriaceae (UID5124)", "f__Streptomycetaceae (UID2048)"],
                                "# genomes": [134, 60], "# markers": [1173, 460], "# marker sets": [336, 233],
                                "0_markers": [79, 24], "1_marker": [1068, 429], "2_markers": [26, 7],
                                "3_markers": [0, 0], "4_markers": [0, 0], "5_markers": [0,0],
                                "Completeness": [93.13, 95.38], "Contamination": [3.11, 1.37],
                                "Strain heterogeneity": [38.46, 42.86]})
        checkm_file = Path('test/data/fake_checkm.txt')
        test_checkm = summary_stats.get_checkm_stats(checkm_file)
        assert test_checkm.equals(true_checkm)


class TestDRepStats(unittest.TestCase):
    def test_get_genomes(self):
        drep_mash = Path('test/data/Mdb.csv')
        true_genomes = pd.DataFrame({"genome2": ["Streptomyces_coelicolor_A32_NC_003888_3.fa",
                                                 "Escherichia_coli_str_K-12_substr_MG1655_NC_000913_3.fa",
                                                 "Klebsiella_pneumoniae_subsp_pneumoniae_CP009208_1.fa",
                                             "Salinispora_tropica_CNB-440_NC_009380_1.fa"]})
        true_bins = pd.DataFrame({"genome2": ["test_hiseq_2500.bin.6.fa",
                                                "test_hiseq_2500.bin.4.fa",
                                                "test_hiseq_2500.bin.9.fa",
                                                "test_hiseq_2500.bin.2.fa"]})
        test_genomes, test_bins = summary_stats.get_bins_and_genomes(drep_mash)
        print(true_genomes.head())
        print(test_genomes.head())
        print(true_bins.head())
        print(test_bins.head())
        assert test_genomes.equals(true_genomes)
        assert test_bins.equals(true_bins)

    def test_get_drep(self):
        true_drep = pd.DataFrame({"query": ["Streptomyces_coelicolor_A32_NC_003888_3.fa",
                                             "Salinispora_tropica_CNB-440_NC_009380_1.fa",
                                             "Escherichia_coli_str_K-12_substr_MG1655_NC_000913_3.fa",
                                             "Klebsiella_pneumoniae_subsp_pneumoniae_CP009208_1.fa",
                                             np.nan],
                                  "reference": ["test_hiseq_2500.bin.2.fa",
                                                "test_hiseq_2500.bin.6.fa",
                                                "test_hiseq_2500.bin.9.fa",
                                                np.nan,
                                                "test_hiseq_2500.bin.4.fa"],
                                  "ref_coverage": [0.9987903, 0.94825554, 0.96552, np.nan, np.nan],
                                  "query_coverage": [0.95382077, 0.98048997, 0.9568264, np.nan, np.nan],
                                  "ani": [0.9999186, 0.9999351, 0.9998399, np.nan, np.nan]})
        drep_mummer = Path('test/data/Ndb.csv')
        drep_mash = Path('test/data/Mdb.csv')
        test_drep = summary_stats.get_drep_stats(drep_mummer, drep_mash)
        pd.set_option('display.max_columns', None)
        print(true_drep.fillna(0).head())
        print(test_drep.fillna(0).head())
        compare_drep = pd.concat([test_drep, true_drep]).reset_index(drop=True)
        compare_groupby = compare_drep.groupby(list(compare_drep.columns))
        idx = [x[0] for x in compare_groupby.groups.values() if len(x) == 1]
        print(compare_drep.reindex(idx))
        assert test_drep.fillna(value=0).equals(true_drep.fillna(value=0))