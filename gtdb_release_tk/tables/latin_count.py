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
from pathlib import PurePath
from collections import defaultdict, namedtuple

from biolib.taxonomy import Taxonomy

import gtdb_release_tk.HTML as HTML
from gtdb_release_tk.taxon_utils import canonical_taxon, is_latin_taxon


class LatinCount(object):
    """Create table indicating percentage of Latin names at each taxonomic rank."""

    def __init__(self, release_number, output_dir):

        self.release_number = release_number
        if output_dir:
            self.output_dir = PurePath(output_dir)

        self.logger = logging.getLogger('timestamp')

    def run(self, bac120_metadata_file, ar120_metadata_file):
        """Create table indicating percentage of Latin names at each taxonomic rank."""

        processed_taxa = set()
        bac_latin = defaultdict(int)
        bac_placeholder = defaultdict(int)
        ar_latin = defaultdict(int)
        ar_placeholder = defaultdict(int)
        for mf in [bac120_metadata_file, ar120_metadata_file]:
            with open(mf, encoding='utf-8') as f:
                header = f.readline().strip().split('\t')

                gtdb_taxonomy_index = header.index('gtdb_taxonomy')
                gtdb_rep_index = header.index('gtdb_representative')

                for line in f:
                    line_split = line.strip().split('\t')

                    gtdb_rep = line_split[gtdb_rep_index]
                    if gtdb_rep != 't':
                        continue

                    taxonomy = line_split[gtdb_taxonomy_index]
                    gtdb_taxa = [t.strip() for t in taxonomy.split(';')]
                    for rank, taxon in enumerate(gtdb_taxa):
                        # suffixed Latin names above the rank of genus
                        # are not considered to be distinct taxa
                        if rank <= Taxonomy.FAMILY_INDEX:
                            taxon = canonical_taxon(taxon)

                        # process each name exactly once
                        if taxon in processed_taxa:
                            continue

                        processed_taxa.add(taxon)

                        if gtdb_taxa[0] == 'd__Bacteria':
                            if is_latin_taxon(taxon):
                                bac_latin[rank] += 1
                            else:
                                bac_placeholder[rank] += 1
                        else:
                            if is_latin_taxon(taxon):
                                ar_latin[rank] += 1
                            else:
                                ar_placeholder[rank] += 1

        # write out table with species counts
        out_prefix = f'gtdb_r{self.release_number}_latin_count'
        fout = open(self.output_dir / f'{out_prefix}.tsv', 'w')

        table_data = []
        header = ['', 'Bacteria: Latin', 'Bacteria: Placeholder',
                  'Archaea: Latin', 'Archaea: Placeholder',
                  'Latin', 'Placeholder']
        table_data.append(header)
        print('\t'.join(header))
        fout.write('\t'.join(header) + '\n')
        for idx in range(1, 7):
            row = []
            row.append(Taxonomy.rank_labels[idx].capitalize())

            bac_total = bac_latin[idx] + bac_placeholder[idx]
            row.append(
                f'{bac_latin[idx]:,} ({100.0*bac_latin[idx]/bac_total:.2f}%)')
            row.append(
                f'{bac_placeholder[idx]:,} ({100.0*bac_placeholder[idx]/bac_total:.2f}%)')

            ar_total = ar_latin[idx] + ar_placeholder[idx]
            row.append(
                f'{ar_latin[idx]:,} ({100.0*ar_latin[idx]/ar_total:.2f}%)')
            row.append(
                f'{ar_placeholder[idx]:,} ({100.0*ar_placeholder[idx]/ar_total:.2f}%)')

            total = bac_total + ar_total
            row.append(
                f'{bac_latin[idx]+ar_latin[idx]:,} ({100.0*(bac_latin[idx]+ar_latin[idx])/total:.2f}%)')
            row.append(
                f'{bac_placeholder[idx]+ar_placeholder[idx]:,} ({100.0*(bac_placeholder[idx]+ar_placeholder[idx])/total:.2f}%)')

            table_data.append(row)
            print('\t'.join(row))
            fout.write('\t'.join(row) + '\n')

        fout.close()

        # create HTML table
        fout = open(self.output_dir / f'{out_prefix}.html', 'w')
        htmlcode = HTML.table(table_data,
                              table_class='taxa_count',
                              col_styles=['font-size: small'] *
                              len(table_data[0]),
                              col_align=['left'] + ['center'] *
                              (len(table_data[0])-1),
                              cellpadding=6)
        fout.write(htmlcode)
        fout.close()
