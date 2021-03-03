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

from gtdb_release_tk.models.taxonomy_rank import TaxonomyRank, RankEnum


class TestTaxonomyRank(unittest.TestCase):

    def test_taxonomy_rank(self):
        tr = TaxonomyRank('d__Archaea')
        self.assertEqual(tr.rank, 'd__Archaea')
        self.assertEqual(tr.name, 'Archaea')
        self.assertEqual(tr.type, RankEnum.DOMAIN)
        self.assertFalse(tr.is_empty())
        self.assertEqual(str(tr), 'd__Archaea')

    def test_taxonomy_rank_when_empty(self):
        tr = TaxonomyRank('g__')
        self.assertEqual(tr.rank, 'g__')
        self.assertIsNone(tr.name)
        self.assertTrue(tr.is_empty())
        self.assertEqual(str(tr), 'g__')

    def test_taxonomy_rank_equal(self):
        tr1 = TaxonomyRank('d__Archaea')
        tr2 = TaxonomyRank('d__Archaea')
        tr3 = TaxonomyRank('d__Bacteria')
        self.assertEqual(tr1, tr2)
        self.assertNotEqual(tr1, tr3)
