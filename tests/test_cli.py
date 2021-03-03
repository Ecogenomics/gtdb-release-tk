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
from gtdb_release_tk.files.taxonomy import ArcTaxonomyFile, BacTaxonomyFile
from tests import maybe_gunzip

DIR_TEST_DATA = os.path.join(os.getcwd(), 'data')


class TestCLI(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ar122_sp_labels_tree = os.path.join(DIR_TEST_DATA, 'ar122.sp_labels.tree')
        cls.ar122_metadata_r95 = os.path.join(DIR_TEST_DATA, 'ar122_metadata_r95.tsv')
        cls.ar122_r95_tree = os.path.join(DIR_TEST_DATA, 'ar122_r95.tree')
        cls.ar122_taxonomy_r95 = os.path.join(DIR_TEST_DATA, 'ar122_taxonomy_r95.tsv')
        cls.bac120_sp_labels_tree = os.path.join(DIR_TEST_DATA, 'bac120.sp_labels.tree')
        cls.bac120_metadata_r95 = os.path.join(DIR_TEST_DATA, 'bac120_metadata_r95.tsv')
        cls.bac120_taxonomy_r95 = os.path.join(DIR_TEST_DATA, 'bac120_taxonomy_r95.tsv')
        cls.bac120_r95_tree = os.path.join(DIR_TEST_DATA, 'bac120_r95.tree')
        cls.both_metadata_r95 = os.path.join(DIR_TEST_DATA, 'both_metadata_r95.tsv')
        cls.both_taxonomy_r95 = os.path.join(DIR_TEST_DATA, 'both_taxonomy_r95.tsv')
        cls.genome_taxonomy_r95_count = os.path.join(DIR_TEST_DATA, 'genome_taxonomy_r95_count.json')
        cls.sp_clusters_r95 = os.path.join(DIR_TEST_DATA, 'sp_clusters_r95.tsv')
        cls.test_metadata = os.path.join(DIR_TEST_DATA, 'test_metadata.tsv')
        cls.test_sp_clusters = os.path.join(DIR_TEST_DATA, 'test_sp_clusters.tsv')
        cls.test_taxonomy = os.path.join(DIR_TEST_DATA, 'test_taxonomy.tsv')

        [maybe_gunzip(v) for k, v in vars(cls).items() if isinstance(v, str) and os.path.isfile(v + '.gz')]

    def setUp(self):
        self.out_dir = tempfile.mkdtemp(prefix='gtdb_release_tk_tmp_')

    def tearDown(self):
        shutil.rmtree(self.out_dir)

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

    def test_tree_files(self):
        args = ['python', '-m', 'gtdb_release_tk', 'tree_files',
                '--metadata_file', self.both_metadata_r95,
                '--bac_tree', self.ar122_r95_tree,
                '--ar_tree', self.bac120_r95_tree,
                '--release_number', '95',
                '--output_dir', self.out_dir]
        p = subprocess.Popen(args, encoding='utf-8')
        p.communicate()
        self.assertEqual(p.returncode, 0)

    def test_sp_cluster_file(self):
        # TODO: Does this do anything now that the UBAs are gone?
        args = ['python', '-m', 'gtdb_release_tk', 'sp_cluster_file',
                '--metadata_file', self.both_metadata_r95,
                '--gtdb_sp_clusters_file', self.sp_clusters_r95,
                '--release_number', '95',
                '--output_dir', self.out_dir]
        p = subprocess.Popen(args, encoding='utf-8')
        p.communicate()
        self.assertEqual(p.returncode, 0)

    def test_hq_genome_file(self):
        args = ['python', '-m', 'gtdb_release_tk', 'hq_genome_file',
                '--metadata_file', self.both_metadata_r95,
                '--release_number', '95',
                '--output_dir', self.out_dir]
        p = subprocess.Popen(args, encoding='utf-8')
        p.communicate()
        self.assertEqual(p.returncode, 0)

    def test_ssu_files(self):
        # args = ['python', '-m', 'gtdb_release_tk', 'ssu_files',
        #         '--metadata_file', self.both_metadata_r95,
        #         '--gtdb_sp_clusters_file', self.sp_clusters_r95,
        #         '--'
        #         '--output_dir', self.out_dir]
        # p = subprocess.Popen(args, encoding='utf-8')
        # p.communicate()
        # self.assertEqual(p.returncode, 0)
        pass


    def test_marker_files(self):
        pass

    def test_protein_files(self):
        pass

    def test_nucleotide_files(self):
        pass

    def test_msa_files(self):
        pass

    def test_gene_files(self):
        pass

    def test_metadata_files(self):
        pass

    def test_dict_file(self):
        args = ['python', '-m', 'gtdb_release_tk', 'dict_file',
                '--taxonomy_file', self.both_taxonomy_r95,
                '--release_number', '95',
                '--output_dir', self.out_dir]
        p = subprocess.Popen(args, encoding='utf-8')
        p.communicate()
        self.assertEqual(p.returncode, 0)


    def test_arb_files(self):
        pass

    def test_validate(self):
        pass

    def test_add_sp_label(self):
        pass

    def test_tax_comp_files(self):
        pass

    def test_json_tree_parser(self):
        args = ['python', '-m', 'gtdb_release_tk', 'json_tree_file',
                '--release_number', '95',
                '--output_dir', self.out_dir,
                '--taxonomy_file', self.both_taxonomy_r95,
                '--metadata_file', self.both_metadata_r95]
        p = subprocess.Popen(args, encoding='utf-8')
        p.communicate()
        self.assertEqual(p.returncode, 0)

    def test_reps_per_rank(self):
        pass

    def test_genomic_stats(self):
        pass

    def test_genome_quality(self):
        args = ['python', '-m', 'gtdb_release_tk', 'genome_quality',
                '--bac120_metadata_file', self.bac120_metadata_r95,
                '--ar122_metadata_file', self.ar122_metadata_r95,
                '--release_number', '95',
                '--output_dir', self.out_dir]
        p = subprocess.Popen(args, encoding='utf-8')
        p.communicate()
        self.assertEqual(p.returncode, 0)

    def test_genome_category_rank(self):
        args = ['python', '-m', 'gtdb_release_tk', 'genome_category_rank',
                '--bac120_metadata_file', self.bac120_metadata_r95,
                '--ar122_metadata_file', self.ar122_metadata_r95,
                '--release_number', '95',
                '--output_dir', self.out_dir]
        p = subprocess.Popen(args, encoding='utf-8')
        p.communicate()
        self.assertEqual(p.returncode, 0)

    def test_nomenclatural_rank(self):
        args = ['python', '-m', 'gtdb_release_tk', 'nomenclatural_rank',
                '--bac120_metadata_file', self.bac120_metadata_r95,
                '--ar122_metadata_file', self.ar122_metadata_r95,
                '--release_number', '95',
                '--output_dir', self.out_dir,
                '--domain', 'Both']
        p = subprocess.Popen(args, encoding='utf-8')
        p.communicate()
        self.assertEqual(p.returncode, 0)

    def test_ncbi_compare(self):
        pass

    def test_sp_rep_type(self):
        args = ['python', '-m', 'gtdb_release_tk', 'sp_rep_type',
                '--bac120_metadata_file', self.bac120_metadata_r95,
                '--ar122_metadata_file', self.ar122_metadata_r95,
                '--release_number', '95',
                '--output_dir', self.out_dir,
                '--domain', 'Both']
        p = subprocess.Popen(args, encoding='utf-8')
        p.communicate()
        self.assertEqual(p.returncode, 0)

        hashes = (('gtdb_r95_sp_rep_type.tsv', '5d166232f8f00adfbf4a2f1fd096c19a1aad6a9d40065493cbd1d092ce4895e5'),
                  ('gtdb_r95_sp_rep_type.png', '4623d021977db20ebc67a45204e1ec95f2595f62da810ebc983d3cb28c93baf7'),
                  ('gtdb_r95_sp_rep_type.svg', '3ff64648b87d027c0db83f41c4e1f988477766eeebdf8cdedeb7d9fcbd78b697'))
        for f_name, true_hash in hashes:
            test_hash = sha256(os.path.join(self.out_dir, f_name))
            self.assertEqual(test_hash, true_hash)

    def test_itol(self):
        pass

    def test_taxa_count(self):
        pass

    def test_top_taxa(self):
        pass

    def test_nomenclatural_check(self):
        pass
