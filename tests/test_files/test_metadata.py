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
import tempfile
import unittest
from gtdb_release_tk.common import canonical_taxon_name
from gtdb_release_tk.files.metadata import MetadataFile


class TestMetadataFile(unittest.TestCase):

    def setUp(self):
        # self.out_dir = tempfile.mkdtemp(prefix='gtdbtk_tmp_')
        pass

    def tearDown(self):
        # shutil.rmtree(self.out_dir)
        pass

    def test_canonical_taxon_name(self):
        mf = MetadataFile.read('/Users/aaron/tmp/ar122_metadata_r95.tsv')



