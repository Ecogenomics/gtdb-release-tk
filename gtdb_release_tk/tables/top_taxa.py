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

import os
import sys
import logging
import argparse
import ntpath
import csv
import random
import operator
from pathlib import PurePath
from collections import defaultdict, namedtuple

from biolib.taxonomy import Taxonomy
from biolib.seq_io import read_seq

import gtdb_release_tk.HTML as HTML
   
class TopTaxa(object):
    """Create table with most abundant taxa at each taxonomic rank."""

    def __init__(self, release_number, output_dir):
        
        self.release_number = release_number
        if output_dir:
            self.output_dir = PurePath(output_dir)

        self.logger = logging.getLogger('timestamp')

    def run(self, gtdb_sp_clusters_file, num_taxa_to_report):
        """Create table with most abundant taxa at each taxonomic rank."""
        
        # parse taxa from GTDB species cluster file
        sp_count = {}
        species_in_taxon = defaultdict(lambda: defaultdict(set))
        with open(gtdb_sp_clusters_file) as f:
            header = f.readline().strip().split('\t')
            
            gtdb_sp_index = header.index('GTDB species')
            gtdb_taxonomy_index = header.index('GTDB taxonomy')
            num_clustered_index = header.index('No. clustered genomes')
            
            for line in f:
                line_split = line.strip().split('\t')
                
                gtdb_sp = line_split[gtdb_sp_index]
                gtdb_taxa = [t.strip() for t in line_split[gtdb_taxonomy_index].split(';')] + [gtdb_sp]
                num_genomes = int(line_split[num_clustered_index])
                
                sp_count[gtdb_sp] = num_genomes
                for rank, taxon in enumerate(gtdb_taxa):
                    species_in_taxon[rank][taxon].add(gtdb_sp)
                
        # count number of species for each taxon
        ranked_sp_count = defaultdict(list)
        for rank in range(6):
            for taxon in sorted(species_in_taxon[rank], key=lambda k: len(species_in_taxon[rank][k]), reverse=True):
                ranked_sp_count[rank].append((taxon, len(species_in_taxon[rank][taxon])))

        # write out table with species counts
        table_data = []
        table_data.append(['Domain', '', 'Phylum', '', 'Class', '', 'Order', '', 'Family', '', 'Genus', ''])
        for idx in range(num_taxa_to_report):
            row = []
            for rank in range(6):
                if len(ranked_sp_count[rank]) > idx:
                    taxon, count = ranked_sp_count[rank][idx]
                    row.append(taxon[3:])
                    row.append(f'{count:,}')
                else:
                    row.append('')
                    row.append('')
                
            table_data.append(row)

        out_prefix = f'gtdb_r{self.release_number}_top_taxa.species'
        fout = open(self.output_dir / f'{out_prefix}.html','w')
        htmlcode = HTML.table(table_data,
                                table_class='top_taxa',
                                col_styles=['font-size: small']*len(table_data[0]),
                                cellpadding=4)
        fout.write(htmlcode)
        fout.close()
        
        # count number of genome for each taxon
        ranked_genome_count = defaultdict(list)
        for rank in range(7):
            for taxon in sorted(species_in_taxon[rank], 
                                key=lambda k: sum([sp_count[sp] for sp in species_in_taxon[rank][k]]), 
                                reverse=True):
                genome_count = sum([sp_count[sp] for sp in species_in_taxon[rank][taxon]])
                ranked_genome_count[rank].append((taxon, genome_count))

        # write out table with species counts
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
            for rank in range(7):
                if len(ranked_genome_count[rank]) > idx:
                    taxon, count = ranked_genome_count[rank][idx]
                    row.append(taxon[3:])
                    row.append(f'{count:,}')
                else:
                    row.append('')
                    row.append('')
                
            table_data.append(row)

        out_prefix = f'gtdb_r{self.release_number}_top_taxa.genomes'
        fout = open(self.output_dir / f'{out_prefix}.html','w')
        htmlcode = HTML.table(table_data,
                                table_class='top_taxa',
                                col_styles=['font-size: small']*len(table_data[0]),
                                cellpadding=4)
        fout.write(htmlcode)
        fout.close()
      