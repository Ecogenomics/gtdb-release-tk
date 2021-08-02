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

from biolib.newick import parse_label

import dendropy

   
class iTOL(object):
    """Generate tree and iTOL files for producing iTOL tree image."""

    def __init__(self, release_number, output_dir):
        
        self.release_number = release_number
        if output_dir:
            self.output_dir = PurePath(output_dir)

        self.logger = logging.getLogger('timestamp')
        
        self.colors = ['#a6cee3',
                        '#1f78b4',
                        '#b2df8a',
                        '#33a02c',
                        '#fb9a99',
                        '#e31a1c',
                        '#fdbf6f',
                        '#ff7f00',
                        '#cab2d6',
                        '#6a3d9a',
                        '#ffff99',
                        '#b15928']

    def run(self, 
                bac120_tree,
                ar122_tree):
        """Generate tree and iTOL files for producing iTOL tree image."""
        
        self.logger.info('Creating trees with iTOL labels.')
        for tree_file, output_tree, itol_colors, itol_labels in [(bac120_tree, 
                                                                    f'bac120_r{self.release_number}.itol.tree', 
                                                                    f'bac120_r{self.release_number}.itol_phyla_colors.txt',
                                                                    f'bac120_r{self.release_number}.itol_phyla_labels.txt'),
                                                                    (ar122_tree, 
                                                                    f'ar122_r{self.release_number}.itol.tree', 
                                                                    f'ar122_r{self.release_number}.itol_phyla_colors.txt',
                                                                    f'ar122_r{self.release_number}.itol_phyla_labels.txt')]:
                                        
            self.logger.info(f'Reading {tree_file} reference tree.')
            domain_tree = dendropy.Tree.get_from_path(tree_file,
                                                   schema='newick',
                                                   rooting='force-rooted',
                                                   preserve_underscores=True)
            self.logger.info(' ...tree contains {:,} genomes.'.format(
                sum([1 for leaf in domain_tree.leaf_node_iter()])))

            phyla = set()
            for node in domain_tree.preorder_node_iter():
                if node.is_leaf():
                    continue
                    
                support, taxon, auxiliary_info = parse_label(node.label)
                if taxon:
                    taxa = taxon.split(';')[0]
                    
                    if taxa.startswith('p__'):
                        node.label = taxa[3:]
                        phyla.add(taxa[3:])
                    else:
                        node.label = None

            domain_tree.write_to_path(self.output_dir / output_tree,
                                   schema='newick',
                                   suppress_rooting=True,
                                   unquoted_underscores=True)
                                   
            self.logger.info('Identified {:,} phyla in {}.'.format(len(phyla), tree_file))
                                   
            # create iTOL metadata for coloring phyla
            self.logger.info(f'Creating iTOL metadata for coloring phyla: {itol_colors}.')
            fout = open(self.output_dir / itol_colors, 'w')
            fout.write('TREE_COLORS\n')
            fout.write('SEPARATOR TAB\n')
            fout.write('DATA\n')
            
            color_index = 0
            for phylum in phyla:
                fout.write('{}\trange\t{}\t{}\n'.format(
                                phylum,
                                self.colors[color_index],
                                phylum))
                                
                color_index += 1
                if color_index >= len(self.colors):
                    color_index = 0
                                
            fout.close()
            
            # create iTOL metadata for phylum labels
            self.logger.info(f'Creating iTOL metadata for phyla labels: {itol_labels}.')
            fout = open(self.output_dir / itol_labels, 'w')
            fout.write('DATASET_TEXT\n')
            fout.write('SEPARATOR TAB\n')
            fout.write('DATASET_LABEL\tPhylum labels\n')
            fout.write('COLOR\t#000000\n')
            fout.write('MARGIN\t0\n')
            fout.write('SHOW_INTERNAL\t1\n')
            fout.write('LABEL_ROTATION\t0\n')
            fout.write('STRAIGHT_LABELS\t0\n')
            fout.write('ALIGN_TO_TREE\t0\n')
            fout.write('SIZE_FACTOR\t1\n')
            fout.write('DATA\n')
            
            color_index = 0
            for phylum in phyla:
                fout.write('{}\t{}\t-1\t{}\tnormal\t1\t0\n'.format(
                                phylum,
                                phylum,
                                self.colors[color_index]))

                color_index += 1
                if color_index >= len(self.colors):
                    color_index = 0
                                
            fout.close()
 