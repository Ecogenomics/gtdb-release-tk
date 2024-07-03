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
<<<<<<< HEAD
import os
from collections import defaultdict
from typing import Dict, Tuple, FrozenSet
=======
from pathlib import PurePath
from collections import defaultdict
>>>>>>> 5fa37f4b23bb84a84a2c6bc4869b662b2b7c8e28

import gtdb_release_tk.HTML as HTML
from gtdb_release_tk.common import summarise_file
from gtdb_release_tk.files.sp_clusters import SpClustersFile
from gtdb_release_tk.models.taxonomy_rank import RANK_ORDERING, RankEnum

logger = logging.getLogger('timestamp')


class TopTaxaFile(object):
    NAME = 'gtdb_r{release}_top_taxa.{type}.{ext}'

    def __init__(self, sp_count, species_in_taxon):
        self.sp_count: Dict[str, int] = sp_count
        self.species_in_taxon: Dict[RankEnum, Dict[str, FrozenSet[str]]] = species_in_taxon

    @classmethod
    def create(cls, sp_clusters_file: SpClustersFile):
        # parse taxa from GTDB species cluster file
        sp_count = dict()
        species_in_taxon = defaultdict(lambda: defaultdict(set))
        for rep_genome, sp_row in sp_clusters_file.rows.items():

            taxonomy = sp_row.full_taxonomy
            num_genomes = sp_row.no_clustered_genomes
            sp_count[taxonomy.s.rank] = num_genomes

            for rank in taxonomy.ranks:
                species_in_taxon[rank.type][rank.rank].add(taxonomy.s.rank)
        # Make immutable
        species_in_taxon = {k: {kk: frozenset(vv) for kk, vv in v.items()} for k, v in species_in_taxon.items()}

        return cls(sp_count, species_in_taxon)

    def get_ranked_species_count(self) -> Dict[RankEnum, Tuple[Tuple[str, int]]]:
        out = defaultdict(list)
        for rank_type in RANK_ORDERING:
            for taxon in sorted(self.species_in_taxon[rank_type],
                                key=lambda k: len(self.species_in_taxon[rank_type][k]),
                                reverse=True):
                out[rank_type].append((taxon, len(self.species_in_taxon[rank_type][taxon])))
        out = {k: tuple(v) for k, v in out.items()}
        return out

    def get_ranked_genomes_count(self) -> Dict[RankEnum, Tuple[Tuple[str, int]]]:
        out = defaultdict(list)
        for rank_type in RANK_ORDERING:
            for taxon in sorted(self.species_in_taxon[rank_type],
                                key=lambda k: sum([self.sp_count[sp] for sp in self.species_in_taxon[rank_type][k]]),
                                reverse=True):
                genome_count = sum([self.sp_count[sp] for sp in self.species_in_taxon[rank_type][taxon]])
                out[rank_type].append((taxon, genome_count))
        out = {k: tuple(v) for k, v in out.items()}
        return out

    def write_genomes(self, out_dir: str, release_no: str, num_taxa_to_report: int):

        path_html = os.path.join(out_dir, self.NAME.format(release=release_no, type='genomes', ext='html'))
        path_tsv = os.path.join(out_dir, self.NAME.format(release=release_no, type='genomes', ext='tsv'))

        # count number of genome for each taxon
        ranked_genome_count = self.get_ranked_genomes_count()

        table_data = []
        table_data.append(['Domain', '',
                           'Phylum', '',
                           'Class', '',
                           'Order', '',
                           'Family', '',
                           'Genus', '',
                           'Species', ''])
        for idx in range(num_taxa_to_report):
            row = []
            for rank_type in RANK_ORDERING:
                if len(ranked_genome_count[rank_type]) > idx:
                    taxon, count = ranked_genome_count[rank_type][idx]
                    row.append(taxon[3:])
                    row.append(f'{count:,}')
                else:
                    row.append('')
                    row.append('')

            table_data.append(row)

        # Write HTML
        with open(path_html, 'w') as f:
            html_code = HTML.table(table_data,
                                   table_class='top_taxa',
                                   col_styles=['font-size: small'] * len(table_data[0]),
                                   cellpadding=4)
            f.write(html_code)
        logger.info(summarise_file(path_html))

        # Write TSV
        with open(path_tsv, 'w') as f:
            for row in table_data:
                f.write('\t'.join(row) + '\n')
        logger.info(summarise_file(path_tsv))

    def write_species(self, out_dir: str, release_no: str, num_taxa_to_report: int):

        path_html = os.path.join(out_dir, self.NAME.format(release=release_no, type='species', ext='html'))
        path_tsv = os.path.join(out_dir, self.NAME.format(release=release_no, type='species', ext='tsv'))

        ranked_sp_count = self.get_ranked_species_count()

        # write out table with species counts
        table_data = []
        table_data.append(['Domain', '', 'Phylum', '', 'Class', '', 'Order', '', 'Family', '', 'Genus', ''])
        for idx in range(num_taxa_to_report):
            row = []
            for rank_type in RANK_ORDERING:
                if rank_type is RankEnum.SPECIES:
                    continue

                if len(ranked_sp_count[rank_type]) > idx:
                    taxon, count = ranked_sp_count[rank_type][idx]
                    row.append(taxon[3:])
                    row.append(f'{count:,}')
                else:
                    row.append('')
                    row.append('')

            table_data.append(row)

        with open(path_html, 'w') as f:
            html_code = HTML.table(table_data,
                                   table_class='top_taxa',
                                   col_styles=['font-size: small'] * len(table_data[0]),
                                   cellpadding=4)
            f.write(html_code)
        logger.info(summarise_file(path_html))

        # Write TSV
        with open(path_tsv, 'w') as f:
            for row in table_data:
                f.write('\t'.join(row) + '\n')
        logger.info(summarise_file(path_tsv))
