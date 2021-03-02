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
from deepdiff import DeepDiff
import json

from gtdb_release_tk.common import sha256
from gtdb_release_tk.files.taxonomy import ArcTaxonomyFile, BacTaxonomyFile
from gtdb_release_tk.files.website_tree_json import WebsiteTreeJsonFile

DIR_TEST_DATA = os.path.join(os.getcwd(), 'data')


class TestCLI(unittest.TestCase):

    def setUp(self):
        self.out_dir = tempfile.mkdtemp(prefix='gtdb_release_tk_tmp_')
        self.test_metadata = os.path.join(DIR_TEST_DATA, 'test_metadata.tsv')
        self.test_sp_clusters = os.path.join(DIR_TEST_DATA, 'test_sp_clusters.tsv')
        self.test_taxonomy = os.path.join(DIR_TEST_DATA, 'test_taxonomy.tsv')

        self.ar122_taxonomy_r95 = os.path.join(DIR_TEST_DATA, 'ar122_taxonomy_r95.tsv')
        self.bac120_taxonomy_r95 = os.path.join(DIR_TEST_DATA, 'bac120_taxonomy_r95.tsv')
        self.sp_clusters_r95 = os.path.join(DIR_TEST_DATA, 'sp_clusters_r95.tsv')
        self.ar122_metadata_r95 = os.path.join(DIR_TEST_DATA, 'ar122_metadata_r95.tsv')
        self.both_metadata_r95 = os.path.join(DIR_TEST_DATA, 'both_metadata_r95.tsv')
        self.both_taxonomy_r95 = os.path.join(DIR_TEST_DATA, 'both_taxonomy_r95.tsv')

        self.genome_taxonomy_r95_count = os.path.join(DIR_TEST_DATA, 'genome_taxonomy_r95_count.json')
        pass


    def tearDown(self):
        shutil.rmtree(self.out_dir)
        pass

    def test_taxonomy_files(self):
        args = ['python', '-m', 'gtdb_release_tk', 'taxonomy_files',
                '--metadata_file', self.both_metadata_r95,
                '--gtdb_sp_clusters_file', self.sp_clusters_r95,
                '--output_dir', self.out_dir,
                '--release_number', '95']
        p = subprocess.Popen(args, encoding='utf-8')
        p.communicate()
        self.assertEqual(p.returncode, 0)

        test_arc_tax = ArcTaxonomyFile.read(os.path.join(self.out_dir, 'ar122_taxonomy_r95.tsv'))
        exp_arc_tax = ArcTaxonomyFile.read(self.ar122_taxonomy_r95)
        self.assertEqual(exp_arc_tax, test_arc_tax)

        test_bac_tax = BacTaxonomyFile.read(os.path.join(self.out_dir, 'bac120_taxonomy_r95.tsv'))
        exp_bac_tax = BacTaxonomyFile.read(self.bac120_taxonomy_r95)
        self.assertEqual(test_bac_tax, exp_bac_tax)

    def test_json_tree_parser(self):
        args = ['python', '-m', 'gtdb_release_tk', 'json_tree_file',
                '--release_number', '95',
                '--output_dir', self.out_dir,
                '--taxonomy_file', self.both_taxonomy_r95,
                '--metadata_file', self.both_metadata_r95]
        p = subprocess.Popen(args, encoding='utf-8')
        p.communicate()
        self.assertEqual(p.returncode, 0)
