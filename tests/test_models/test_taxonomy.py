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

from gtdb_release_tk.exceptions import InvalidTaxonomy
from gtdb_release_tk.models.taxonomy import Taxonomy


class TestTaxonomy(unittest.TestCase):

    def test_taxonomy(self):
        test = 'd__D;p__P;c__C;o__O;f__F;g__g;s__S'
        ts = Taxonomy(test)
        self.assertEqual(str(ts), test)

    def test_taxonomy_equal(self):
        test = 'd__D;p__P;c__C;o__O;f__F;g__g;s__S'
        ts_1 = Taxonomy(test)
        ts_2 = Taxonomy(test)
        self.assertEqual(ts_1, ts_2)

    def test_taxonomy_empty(self):
        test = 'd__D;p__;c__C;o__O;f__F;g__g;s__S'
        ts = Taxonomy(test)
        self.assertTrue(ts.empty(ts.p))
        self.assertFalse(ts.empty(ts.c))

    def test_ranks(self):
        ranks = ('d__D', 'p__P', 'c__C', 'o__O', 'f__F', 'g__G', 's__')
        ts = Taxonomy('; '.join(ranks))
        for test_rank, exp_rank in zip(ts.ranks, ranks):
            self.assertEqual(test_rank, exp_rank)

    def test_taxonomy_raises_InvalidTaxonomy(self):
        test = 'd__D;p__;c__C;o__O;f__F;g__g;'
        self.assertRaises(InvalidTaxonomy, Taxonomy, test)
