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

import unittest

from gtdb_release_tk.files.gtdb_dict import GTDBDictFile
from gtdb_release_tk.files.taxonomy import TaxonomyFile
from gtdb_release_tk.models.taxonomy_rank import TaxonomyRank
from gtdb_release_tk.models.taxonomy_string import TaxonomyString


class TestMetadataFile(unittest.TestCase):

    def test_create(self):
        tf_data = {'a': TaxonomyString('d__D;p__P;c__C;o__O;f__F;g__G;s__Foo bar')}
        tf = TaxonomyFile(tf_data)
        df = GTDBDictFile.create(tf)
        expected = {'d__D', 'D', 'p__P', 'P', 'c__C', 'C', 'o__O', 'O', 'f__F',
                    'F', 'g__G', 'G', 's__Foo', 'bar', 'F. bar', 'Foo', 's__F. bar',
                    'Foo bar', 's__Foo bar'}
        self.assertSetEqual(df.data, expected)

    def test_parse_rank(self):
        tr = TaxonomyRank('d__Domain')
        self.assertEqual(GTDBDictFile.parse_rank(tr), ('Domain', 'd__Domain'))

    def test_parse_rank_species(self):
        tr = TaxonomyRank('s__Foo bar')
        test = set(GTDBDictFile.parse_rank_species(tr))
        expected = {'s__Foo', 'bar', 'F. bar', 'Foo', 's__F. bar'}
        self.assertSetEqual(test, expected)
