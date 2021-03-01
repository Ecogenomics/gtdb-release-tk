###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import os
import shutil
import subprocess
import tempfile
import unittest

from gtdb_release_tk.common import sha256

DIR_TEST_DATA = os.path.join(os.getcwd(), 'data')


class TestCLI(unittest.TestCase):

    def setUp(self):
        self.out_dir = tempfile.mkdtemp(prefix='gtdb_release_tk_tmp_')
        self.test_metadata = os.path.join(DIR_TEST_DATA, 'test_metadata.tsv')
        self.test_taxonomy = os.path.join(DIR_TEST_DATA, 'test_taxonomy.tsv')
        pass

    def tearDown(self):
        shutil.rmtree(self.out_dir)
        pass

    def test_json_tree_parser(self):
        args = ['python', '-m', 'gtdb_release_tk', 'json_tree_file',
                '--release_number', '95',
                '--output_dir', self.out_dir,
                '--taxonomy_file', self.test_taxonomy,
                '--metadata_file', self.test_metadata]
        p = subprocess.Popen(args, encoding='utf-8')
        p.communicate()
        self.assertEqual(p.returncode, 0)

        path_out = os.path.join(self.out_dir, 'genome_taxonomy_r95_count.json')
        test_hash = sha256(path_out)
        exp_hash = '6699c5f8c24c5fa70062f9a98d0806278be5030cb75d5e44246678bd1f386e45'
        self.assertEqual(test_hash, exp_hash)
