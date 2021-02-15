import unittest

import run_magician


class TestRunMagician(unittest.TestCase):
    def test_bad_file_name(self):
        bad_result = "test_ä.txt"
        with self.assertRaisesRegex(ValueError, r"""Arguments can only consist of alphanumeric characters, quote marks, \., \/, \_, \- and space."""):
            run_magician.run_snakemake(bad_result)

    def test_missing_profile_info(self):
        result = "test.txt"
        profile = "own"
        with self.assertRaisesRegex(ValueError,
                                    "Both name of the custom error profile and read length of the error profile must be given when using own profiles."):
            run_magician.run_snakemake(result, profile)

    def test_too_much_info(self):
        result = "test.txt"
        profile = "mbarc"
        extra_name = "TestR"
        with self.assertRaisesRegex(ValueError,
                                    "Name of the error profile and read length can only be specified when using own profiles."):
            run_magician.run_snakemake(result, profile, extra_name)