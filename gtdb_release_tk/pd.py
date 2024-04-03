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
__copyright__ = 'Copyright 2021'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '0.0.1'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import sys
import logging

import dendropy

from gtdb_release_tk.common import (parse_gtdb_sp_clusters,
                                    parse_genome_category,
                                    sp_cluster_type_category)


class PhylogeneticDiversity(object):
    """Phylogenetic diversity of isolate vs. environmental species clusters."""

    def __init__(self):
        """Initialization."""

        self.logger = logging.getLogger('timestamp')

    def pd(self, tree_file, sp_types):
        """Calculate phylogenetic diversity of isolate vs. environmental species clusters."""

        tree = dendropy.Tree.get_from_path(tree_file,
                                           schema='newick',
                                           rooting='force-rooted',
                                           preserve_underscores=True)
        self.logger.info(' - tree contains {:,} genomes'.format(
            sum([1 for _ in tree.leaf_node_iter()])))

        pd_shared = 0
        pd_isolate = 0
        pd_env = 0
        for node in tree.preorder_node_iter():
            if node == tree.seed_node:
                continue

            node_sp_types = set()
            for leaf in node.leaf_iter():
                rid = leaf.taxon.label
                node_sp_types.add(sp_types[rid])

            if 'BOTH' in node_sp_types or len(node_sp_types) >= 2:
                pd_shared += node.edge_length
            elif list(node_sp_types)[0] == 'ISOLATE' and len(node_sp_types) == 1:
                pd_isolate += node.edge_length
            elif list(node_sp_types)[0] == 'ENV' and len(node_sp_types) == 1:
                pd_env += node.edge_length
            else:
                print(f'Error in pd: {node_sp_types}')
                sys.exit(-1)

        return pd_shared, pd_isolate, pd_env

    def run(self, bac120_metadata_file, ar122_metadata_file, bac120_tree, ar122_tree):
        """Phylogenetic diversity of isolate vs. environmental species clusters."""

        for domain, metadata_file, tree_file in (('Archaea', ar122_metadata_file, ar122_tree), ('Bacteria', bac120_metadata_file, bac120_tree)):
            # determine category (isolate or environmental) for each species cluster
            self.logger.info(
                f'Parsing metadata for {domain} species clusters:')

            sp_clusters = parse_gtdb_sp_clusters(metadata_file)
            self.logger.info(
                f' - identified {len(sp_clusters):,} species clusters')

            genome_cat = parse_genome_category(metadata_file)
            self.logger.info(
                f' - determined genome category for {len(genome_cat):,} genomes')
            self.logger.info(' - identified the following categories: {}'.format(
                set(genome_cat.values())))

            sp_types = {}
            for rid, gids in sp_clusters.items():
                sp_types[rid] = sp_cluster_type_category(gids, genome_cat)

            # calculate phylogenetic diversity
            self.logger.info('Calculating phylogenetic diversity:')
            pd_shared, pd_isolate, pd_env = self.pd(tree_file, sp_types)
            total_pd = pd_shared + pd_isolate + pd_env
            self.logger.info(
                ' - total: {:.2f}'.format(total_pd))
            self.logger.info(
                ' - shared: {:.2f} ({:.2f}%)'.format(pd_shared, 100*pd_shared/total_pd))
            self.logger.info(
                ' - isolate: {:.2f} ({:.2f}%)'.format(pd_isolate, 100*pd_isolate/total_pd))
            self.logger.info(
                ' - environmental: {:.2f} ({:.2f}%)'.format(pd_env, 100*pd_env/total_pd))
