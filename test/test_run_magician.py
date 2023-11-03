import pathlib
import unittest

import run_magician


class TestRunMagician(unittest.TestCase):
    snake_path = pathlib.Path(__file__).resolve().parent.parent / "snakefiles" / "Snakefile"
    profile_type = "mbarc"
    profile_base = ""
    readlength = ""
    insert_size = 270
    cluster_cmd = ""
    def test_default_command_success(self):
        """Run Snakemake with default settings."""
        expected_command = ["snakemake", "all_bin_summaries", "-s", self.snake_path,
                            "--config", 'profile_type="mbarc"',
                            'profile_name="False"', 'readlength="False"',
                            'insert_size=270', "-n"]
        snake_flags = ["-n"]
        test_command = run_magician.get_snake_cmd("all_bin_summaries", self.profile_type, self.profile_base,
                                                  self.readlength, self.insert_size, self.cluster_cmd, *snake_flags)
        assert test_command == expected_command


    def test_run_cluster_command(self):
        """Run Snakemake in cluster mode"""
        expected_command = ["snakemake", "all_bin_summaries", "-s", self.snake_path,
                            "--cluster", "qsub -pe threaded {threads}",
                            "--config", 'profile_type="mbarc"',
                            'profile_name="False"', 'readlength="False"',
                            'insert_size=270', "-n"]
        snake_flags = ["-n"]
        cluster_cmd = "qsub -pe threaded {threads}"
        test_command = run_magician.get_snake_cmd("all_bin_summaries", self.profile_type, self.profile_base,
                                                  self.readlength, self.insert_size, cluster_cmd,
                                                  *snake_flags)
        assert test_command == expected_command

    def test_change_settings(self):
        """Use non-default settings for various parameters."""
        expected_command = ["snakemake", "all_bin_summaries", "-s", self.snake_path,
                            "--config", 'profile_type="mbarc"',
                            'profile_name="False"', 'readlength="False"',
                            'insert_size=500', "-n"]
        snake_flags = ["-n"]
        insert_size = 500
        test_command = test_command = run_magician.get_snake_cmd("all_bin_summaries", self.profile_type,
                                                                 self.profile_base,
                                                                 self.readlength, insert_size, self.cluster_cmd,
                                                                 *snake_flags)
        assert test_command == expected_command

    def test_bad_file_name(self):
        bad_result = "test_ä.txt"
        with self.assertRaisesRegex(ValueError,
                                    r"""Arguments can only consist of alphanumeric characters, quote marks, \., \/, \_, \- and space."""):
            run_magician.get_snake_cmd(bad_result)

    def test_missing_profile_info(self):
        result = "test.txt"
        profile = "own"
        with self.assertRaisesRegex(ValueError,
                                    "Both name of the custom error profile and read length of the error profile must be given when using own profiles."):
            run_magician.get_snake_cmd(result, profile)

    def test_too_much_info(self):
        result = "test.txt"
        profile = "mbarc"
        extra_name = "TestR"
        with self.assertRaisesRegex(ValueError,
                                    "Name of the error profile and read length can only be specified when using own profiles."):
            run_magician.get_snake_cmd(result, profile, extra_name)