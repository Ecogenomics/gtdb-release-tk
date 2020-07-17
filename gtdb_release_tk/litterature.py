#!/usr/bin/env python3

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

__prog_name__ = 'litterature.py'
__prog_desc__ = 'Parse NCBI taxonomy files to produce simplified summary files.'

__author__ = 'Pierre Chaumeil'
__copyright__ = 'Copyright 2020'
__credits__ = ['Donovan Parks', 'Pierre Chaumeil']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Pierre Chaumeil'
__email__ = 'uqpchaum@uq.edu.au'
__status__ = 'Development'

import os
import sys
import argparse
import traceback
import re
from collections import namedtuple, defaultdict
import numpy as np

from biolib.common import check_file_exists, make_sure_path_exists, check_dir_exists


class TaxonomyNCBI(object):
    """Parse NCBI taxonomy files to produce a simplified summary file."""

    def __init__(self, outfile):
        self.outfile = outfile
        self.NodeRecord = namedtuple(
            'NodeRecord', 'parent_tax_id rank division_id genetic_code_id')
        self.NameRecord = namedtuple('NamesRecord', 'name_txt')

        self.bacterial_division = '0'
        self.unassigned_division = '8'

    def _read_nodes(self, nodes_file):
        """Read NCBI nodes.dmp file.

        Parameters
        ----------
        nodes_file : str
          Path to NCBI nodes.dmp file.

        Returns
        -------
        dict : d[tax_id] -> NodeRecord
          Node record for all nodes.
        """

        d = {}
        for line in open(nodes_file):
            line_split = [t.strip() for t in line.split('|')]

            tax_id = line_split[0]
            parent_tax_id = line_split[1]
            rank = line_split[2]
            division_id = line_split[4]
            genetic_code_id = line_split[6]

            d[tax_id] = self.NodeRecord(
                parent_tax_id, rank, division_id, genetic_code_id)

        return d

    def file_len(self, fname):
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 1

    def _read_names(self, names_file):
        """Read NCBI names.dmp file.

        Parameters
        ----------
        names_file : str
          Path to NCBI names.dmp file.

        Returns
        ------- 
        dict : d[tax_id] -> NameRecord
          Name record of nodes marked as 'scientific name'.
        """

        name_classes = []

        d = {}

        count_total = self.file_len(names_file)

        count = 0
        for line in open(names_file):
            count += 1
            print('{}/{}'.format(count, count_total), end='\r')
            line_split = [t.strip() for t in line.split('|')]

            tax_id = line_split[0]
            name_txt = line_split[1]
            unique_name = line_split[2]
            name_class = line_split[3]

            name_classes.append(name_class)

            # if name_class == 'scientific name' or name_class == 'synonym':
            d[tax_id] = self.NameRecord(name_txt)

        for it in set(name_classes):
            print(f'{it}')

        return d

    def _get_leaves_id(self, records):
        tax_ids = []
        parent_ids = []
        for tax_id, record in records.items():
            tax_ids.append(tax_id)
            parent_ids.append(record.parent_tax_id)

        leaf_list = np.setdiff1d(tax_ids, parent_ids)
        print('Read %d leaves.' % len(leaf_list))
        return leaf_list

    def run(self, ncbi_node_file, ncbi_name_file):
        """Read NCBI taxonomy information and create summary output files."""

        node_records = self._read_nodes(ncbi_node_file)
        print('Read %d node records.' % len(node_records))

        name_records = self._read_names(ncbi_name_file)
        print('Read %d name records.' % len(name_records))

        list_id_leafs = self._get_leaves_id(node_records)
        full_ranks = []

        count_id = 0
        node_records_temp = dict(node_records)

        for tax_id, node_info in node_records.items():
            count_id += 1
            print('{}/{}'.format(count_id, len(node_records)), end='\r')
            tax_id_ranks = []
            current_node = node_info
            to_add = False
            while current_node.parent_tax_id != '1':
                tax_id_ranks.append(current_node.parent_tax_id)
                current_node = node_records_temp.get(
                    current_node.parent_tax_id)
            tax_id_ranks.append(tax_id)
            if '2' in tax_id_ranks or '2157' in tax_id_ranks:
                full_ranks.extend(tax_id_ranks)

        set_full_ranks = set(full_ranks)

        outf = open(self.outfile, 'w')

        list_ranks = []
        for inc, rank_id in enumerate(set_full_ranks):
            print('{}/{}'.format(inc, len(set_full_ranks)), end='\r')
            if rank_id in name_records:
                outf.write('{}\tNCBI\n'.format(
                    name_records.get(rank_id).name_txt))
                list_ranks.append(name_records.get(rank_id).name_txt)
        outf.close()
        return list_ranks


class LitteratureParser(object):
    def __init__(self, output_dir):
        self.outdir = output_dir
        make_sure_path_exists(self.outdir)
        self.outfile = os.path.join(self.outdir, 'existing_names.tsv')

    def parse_nom_file(self, lpsn_file, source):
        outf = open(self.outfile, 'a')
        list_ranks = []
        with open(lpsn_file) as lf:
            lf.readline()
            for line in lf:
                spename = line.strip().split('\t')[0].replace('s__', '')
                list_ranks.append(spename)
                outf.write(f'{spename}\t{source}\n')
        outf.close()
        return list_ranks

    def run(self, ncbi_node_file, ncbi_name_file, lpsn_species_file,
            bacdive_species_file, gtdb_taxonomy):
        tn = TaxonomyNCBI(self.outfile)
        ncbi_ranks = tn.run(ncbi_node_file, ncbi_name_file)
        lpsn_ranks = self.parse_nom_file(lpsn_species_file, 'LPSN')
        bacdive_ranks = self.parse_nom_file(bacdive_species_file, 'DSMZ')

        all_ranks = ncbi_ranks + lpsn_ranks + bacdive_ranks

        gtdb_name_in_ncbi = open(os.path.join(
            self.outdir, 'gtdb_name_in_litterature.txt'), 'w')
        gtdb_name_notin_ncbi = open(os.path.join(
            self.outdir, 'gtdb_name_not_in_litterature.txt'), 'w')

        dict_gtdbranks = {}
        with open(gtdb_taxonomy) as gtdbtax:
            gtdb_ranks = []
            for line in gtdbtax:
                taxline = line.strip().split('\t')[1].split(';')
                for idx, tax in enumerate(taxline):
                    gtdb_ranks.append(tax[3:])
                    dict_gtdbranks[tax[3:]] = ';'.join(taxline[0:idx + 1])
        gtdb_ranks = set(gtdb_ranks)

        count_rank = 0
        for rk in gtdb_ranks:
            count_rank += 1
            print('{}/{} of gtdb ranks'.format(count_rank, len(gtdb_ranks)), end='\r')
            if re.search("^[A-Za-z\s]+$", rk):
                if any(rk in s for s in all_ranks):
                    gtdb_name_in_ncbi.write(
                        '{}\t{}\n'.format(rk, dict_gtdbranks.get(rk)))
                else:
                    gtdb_name_notin_ncbi.write(
                        '{}\t{}\n'.format(rk, dict_gtdbranks.get(rk)))

        gtdb_name_notin_ncbi.close()
        gtdb_name_in_ncbi.close()


if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)
    print('  by ' + __author__ + ' (' + __email__ + ')' + '\n')

    try:
        p = TaxonomyNCBI()
        p.run()
    except SystemExit:
        print("\nControlled exit resulting from an unrecoverable error or warning.")
    except:
        print("\nUnexpected error:", sys.exc_info()[0])
        raise
