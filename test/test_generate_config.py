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
