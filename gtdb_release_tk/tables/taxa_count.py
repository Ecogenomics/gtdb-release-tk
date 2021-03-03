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

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2019'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import logging
import os
from collections import defaultdict
from pathlib import PurePath

from biolib.taxonomy import Taxonomy

import gtdb_release_tk.HTML as HTML
from gtdb_release_tk.common import canonical_taxon_name, summarise_file


class TaxaCount(object):
    """Create table with number of taxa at each taxonomic rank."""

    def __init__(self, release_number, output_dir):

        self.release_number = release_number
        if output_dir:
            self.output_dir = PurePath(output_dir)

        self.logger = logging.getLogger('timestamp')

    def run(self, bac120_metadata_file, ar120_metadata_file):
        """Create table with number of taxa at each taxonomic rank."""

        bac_taxa = defaultdict(set)
        ar_taxa = defaultdict(set)
        for mf in [bac120_metadata_file, ar120_metadata_file]:
            with open(mf, encoding='utf-8') as f:
                header = f.readline().strip().split('\t')

                gtdb_taxonomy_index = header.index('gtdb_taxonomy')
                gtdb_rep_index = header.index('gtdb_representative')

                for line in f:
                    line_split = line.strip().split('\t')

                    gid = line_split[0]

                    gtdb_rep = line_split[gtdb_rep_index]
                    if gtdb_rep != 't':
                        continue

                    taxonomy = line_split[gtdb_taxonomy_index]
                    gtdb_taxa = [t.strip() for t in taxonomy.split(';')]
                    for rank, taxon in enumerate(gtdb_taxa):
                        if rank == 0:
                            continue  # skip domain

                        if rank < 5:
                            # canonicalize names above genus
                            taxon = canonical_taxon_name(taxon)

                        if gtdb_taxa[0] == 'd__Bacteria':
                            bac_taxa[rank].add(taxon)
                        else:
                            ar_taxa[rank].add(taxon)

        # write out table with species counts
        table_data = []
        table_data.append(['', 'Bacteria', 'Archaea', 'Total'])
        print('Rank\tBacteria\tArchaea\tTotal')
        for idx in range(1, 7):
            row = []
            row.append(Taxonomy.rank_labels[idx].capitalize())
            row.append(f'{len(bac_taxa[idx]):,}')
            row.append(f'{len(ar_taxa[idx]):,}')
            row.append(f'{len(bac_taxa[idx]) + len(ar_taxa[idx]):,}')

            table_data.append(row)
            print('\t'.join(row))

        out_prefix = f'gtdb_r{self.release_number}_taxa_count'
        path_html = os.path.join(self.output_dir, f'{out_prefix}.html')
        with open(path_html, 'w') as fout:
            htmlcode = HTML.table(table_data,
                                  table_class='taxa_count',
                                  col_styles=['font-size: small'] * len(table_data[0]),
                                  col_align=['left'] + ['center'] * (len(table_data[0]) - 1),
                                  cellpadding=6)
            fout.write(htmlcode)
        self.logger.info(summarise_file(path_html))
