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

import logging
import os
from collections import defaultdict

import gtdb_release_tk.HTML as HTML
from gtdb_release_tk.common import summarise_file
from gtdb_release_tk.files.metadata import MetadataFile
from gtdb_release_tk.models.taxonomy_rank import RankEnum, RANK_ORDERING

logger = logging.getLogger('timestamp')


class TaxaCountFile(object):
    __slots__ = ('arc_ranks', 'bac_ranks')

    NAME = 'gtdb_r{release}_taxa_count.html'

    def __init__(self, arc_ranks, bac_ranks):
        self.arc_ranks = {k: frozenset(v) for k, v in arc_ranks.items()}
        self.bac_ranks = {k: frozenset(v) for k, v in bac_ranks.items()}

    @classmethod
    def create(cls, bac_mf: MetadataFile, arc_mf: MetadataFile):
        bac_ranks = defaultdict(set)
        arc_ranks = defaultdict(set)
        for mf in (bac_mf, arc_mf):
            for gid, meta in mf.rows.items():

                gtdb_rep = meta.gtdb_representative
                if gtdb_rep != 't':
                    continue

                tax_string = meta.gtdb_taxonomy
                for rank in tax_string.ranks:
                    if rank.type is RankEnum.DOMAIN:
                        continue

                    if rank.type not in (RankEnum.GENUS, RankEnum.SPECIES):
                        rank_str = rank.canonical
                    else:
                        rank_str = rank.rank

                    if tax_string.d.rank == 'd__Bacteria':
                        bac_ranks[rank.type].add(rank_str)
                    else:
                        arc_ranks[rank.type].add(rank_str)
        return cls(arc_ranks, bac_ranks)

    def as_table(self):
        out = list()
        out.append(['', 'Bacteria', 'Archaea', 'Total'])
        print('Rank\tBacteria\tArchaea\tTotal')
        for rank_type in RANK_ORDERING:
            if rank_type is RankEnum.DOMAIN:
                continue
            row = list()
            row.append(rank_type.name.capitalize())
            row.append(f'{len(self.bac_ranks[rank_type]):,}')
            row.append(f'{len(self.arc_ranks[rank_type]):,}')
            row.append(f'{len(self.bac_ranks[rank_type]) + len(self.arc_ranks[rank_type]):,}')
            out.append(row)
            print('\t'.join(row))
        return out

    def write_html(self, root_dir: str, release: str):
        path = os.path.join(root_dir, self.NAME.format(release=release))
        table_data = self.as_table()
        with open(path, 'w') as f:
            htmlcode = HTML.table(table_data,
                                  table_class='taxa_count',
                                  col_styles=['font-size: small'] * len(table_data[0]),
                                  col_align=['left'] + ['center'] * (len(table_data[0]) - 1),
                                  cellpadding=6)
            f.write(htmlcode)
        logger.info(summarise_file(path))
