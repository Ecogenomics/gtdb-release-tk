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
import unittest

from gtdb_release_tk.files.taxonomy import TaxonomyFile, TaxonomyString

DIR_TEST_DATA = os.path.join(os.path.dirname(os.getcwd()), 'data')


def load_tax_file(path):
    out = dict()
    with open(path) as f:
        for line in f.readlines():
            gid, tax = line.strip().split('\t')
            out[gid] = tax
    return out


class TestTaxonomyFile(unittest.TestCase):

    def setUp(self):
        # Load the test taxonomy file
        self.taxonomy_file = os.path.join(DIR_TEST_DATA, 'test_taxonomy.tsv')
        self.test_tax = load_tax_file(self.taxonomy_file)

    def tearDown(self):
        pass

    def test_taxonomy_file(self):
        tf = TaxonomyFile.read(self.taxonomy_file)
        test = {k: str(v) for k, v in tf.data.items()}
        self.assertDictEqual(test, self.test_tax)

    def test_taxonomy_string(self):
        test = 'd__D;p__P;c__C;o__O;f__F;G__g;s__S'
        ts = TaxonomyString(test)
        self.assertEqual(str(ts), test)
